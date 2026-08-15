#!/usr/bin/env python3

"""Transcriptome read rescue for the streaming quantification pass.

The batch rescue in IsoformReadRescue is shaped around having every candidate at once:
it collects a read-name set, writes a transcript FASTA and a reads FASTA, and shells out
to the minimap2 CLI.  A streaming pass has one read in hand at a time and forgets it, so
that shape does not transfer -- and skipping rescue instead is not a smaller version of
the same feature.  On ONT chr20, 14,455 of 120,370 records are rescue candidates and
10,119 of those are reads the extractor discarded for low percent identity, which the
streaming loop never even saw.

So this keeps a resident minimap2 index in-process, via mappy -- the Python binding that
ships inside the minimap2 repo and calls mm_map() directly -- and aligns each candidate
as it arrives.

Deliberately NOT reimplemented here:

  - the candidate population.  The categories are IsoformReadRescue's, and the branch
    that produces them lives in StreamingQuant.stream_assign beside the discard test it
    has to agree with.
  - the transcript targets.  _build_transcript_models() from the batch path, over the
    transcripts of the current splice graph, so a read can rescue only to the targets it
    could rescue to today.  A global transcriptome index would make the comparison
    against the batch path meaningless.
  - acceptance and projection.  mappy hits are converted into pysam records and pushed
    through _parse_rescue_alignments() itself.  A second, simpler acceptance path is how
    the two would silently drift apart, and the rules there are load-bearing enough that
    reproducing them by hand would not survive one edit to either copy.

Index residency and the GIL
---------------------------
Alignment in mappy holds the GIL -- there is no nogil block in mappy.pyx -- so threads
serialize on mm_map() and buy nothing.  That is not a problem here: LRAA already forks a
process per contig-strand, and each worker builds its own index over only that
contig-strand's transcripts.  One small index per process, one mm_map() at a time per
process, no shared aligner and therefore no contention.  The socket-server design the
investigation notes describe is for a single index shared by many clients and is not
needed at this granularity.

Known parity gaps against the batch invocation
----------------------------------------------
The batch path runs `minimap2 -a -N 50 --secondary=yes -f <frac>`.  Reproduced here:
`-N 50` is Aligner(best_n=50); `--secondary=yes` is minimap2's default and mappy does
not set MM_F_NO_PRINT_2ND, so all secondary hits come back; `-a` is implicit, since
mappy always sets MM_F_CIGAR.  Two things are NOT reproducible through the binding:

  - `-f`.  mappy exposes no equivalent, and mm_mapopt_t is a cdef attribute of Aligner
    with no Python setter.  `-f 0` resolves to mid_occ = max_mid_occ = 1e6, i.e. no
    minimizer-occurrence filtering, while the default 2e-4 fraction is floored at
    min_mid_occ = 10.  On an index of one contig-strand's transcripts, whether that
    differs is measurable rather than arguable; see the report.
  - the alignment score.  mappy's Alignment carries mlen, blen and NM but not
    p->dp_score, which is what SAM's AS tag holds.  _alignment_score() therefore takes
    its documented no-AS branch, matched-minus-NM, for streaming hits and the DP score
    for batch hits.  That changes only which hits tie for best, so it changes acceptance
    only where require_unique_path_across_best_hits has to arbitrate.

Everything else about a hit -- reference, coordinates, CIGAR, NM, secondary and
supplementary status -- is reproduced exactly, including minimap2's rule for which of
several primary hits is the SAM primary and which are supplementary.
"""

import logging
import os
import tempfile
import time
from collections import defaultdict

import pysam

import IsoformReadRescue
import LRAA_Globals
import Util_funcs

logger = logging.getLogger(__name__)


# CIGAR operation codes.  mappy emits (length, op) with op in {0:M, 1:I, 2:D, 3:N};
# pysam wants (op, length).  Named rather than inlined because the two libraries order
# the pair differently and a silent swap produces a plausible-looking wrong alignment.
BAM_CSOFT_CLIP = 4

# SAM flag bits set on the converted records.
SAM_FLAG_REVERSE = 16
SAM_FLAG_SECONDARY = 256
SAM_FLAG_SUPPLEMENTARY = 2048


def mappy_available():
    """Whether the mappy binding can be imported, without raising if it cannot.

    Separate from the import inside StreamingRescuer so the CLI can refuse a run up
    front rather than after the splice graph is built, per contig-strand, in a worker.
    """
    try:
        import mappy  # noqa: F401
    except ImportError:
        return False
    return True


class RescueCandidateCensus:
    """Records what the streaming pass offers rescue, aligning nothing.

    The reproduction control for the candidate set.  It rides the same offer() hook the
    real rescuer does, so what it reports is the population production gating actually
    produces rather than a re-derivation of that gating -- which is the only version of
    this measurement worth having, since a re-derivation agreeing with itself proves
    nothing.
    """

    def __init__(self):
        self.candidate_read_names = defaultdict(set)

    def offer(self, read, category):
        self.candidate_read_names[category].add(
            Util_funcs.get_read_name_include_sc_encoding(read)
        )
        return ()

    def close(self):
        return

    def stats(self):
        return {
            category: len(read_names)
            for category, read_names in self.candidate_read_names.items()
        }


def write_candidate_dump(out_path, candidate_read_names_by_category):
    """read_name<TAB>category, sorted, one line per (read, category).

    Sorted so two dumps of the same population are byte-identical regardless of the
    order reads arrived in, which is what makes `diff` a usable comparison and what
    keeps the artifact from depending on worker completion order.
    """
    rows = sorted(
        (read_name, category)
        for category, read_names in candidate_read_names_by_category.items()
        for read_name in read_names
    )
    with open(out_path, "wt") as ofh:
        for read_name, category in rows:
            ofh.write(f"{read_name}\t{category}\n")
    logger.info("wrote %d rescue candidate rows to %s", len(rows), out_path)
    return len(rows)


class StreamingRescuer:
    """A resident minimap2 index over one contig-strand's transcripts.

    Construct once per contig-strand worker, call offer() per candidate read, close()
    when the stream ends.  offer() returns the splice-graph paths the read rescued to,
    which is zero or one path under the rescue acceptance rules -- the caller looks each
    up in the assignment table exactly as it does a genomically derived path, so nothing
    about how a rescued read is quantified is special-cased.
    """

    def __init__(
        self,
        splice_graph,
        transcripts,
        contig_seq_str,
        read_path_mapper,
        include_monoexonic=False,
        tmp_dir=None,
        record_candidate_names=False,
    ):
        import mappy

        self._splice_graph = splice_graph
        self._read_path_mapper = read_path_mapper
        self._transcript_models = IsoformReadRescue._build_transcript_models(
            splice_graph,
            transcripts,
            contig_seq_str,
            include_monoexonic=include_monoexonic,
        )
        self._aligner = None
        self._buf = None
        self._header = None
        self._target_id_to_tid = {}
        self.reads_offered = defaultdict(int)
        self.reads_rescued = defaultdict(int)
        self.reads_without_sequence = 0
        self.hits_examined = 0
        self.align_seconds = 0.0
        self.accept_seconds = 0.0
        # Off by default: at a billion reads a candidate name set is exactly the
        # unbounded per-read state this whole pass exists not to hold. Turned on only by
        # the candidate-parity measurement, where the population is 14k names.
        self._record_candidate_names = record_candidate_names
        self.candidate_read_names = defaultdict(set)

        if not self._transcript_models:
            # Not an error: a contig-strand can hold no multi-exonic target, and the
            # batch path returns early on the same condition.  offer() then declines
            # everything, which is what having no target to rescue against means.
            logger.info(
                "[%s%s] streaming rescue: no transcript models to align against",
                splice_graph.get_contig_acc(),
                splice_graph.get_contig_strand(),
            )
            return

        transcript_fa = self._write_transcript_fasta(tmp_dir)
        try:
            preset = IsoformReadRescue._resolve_rescue_minimap2_preset()
            self._aligner = mappy.Aligner(
                fn_idx_in=transcript_fa,
                preset=preset,
                # -N 50, the batch invocation's secondary-hit budget.  mappy's default
                # is 5, which would silently discard most of what the batch path ranks.
                best_n=50,
                n_threads=max(1, int(LRAA_Globals.config.get("tool_threads", 1) or 1)),
            )
            # Load failure is falsy rather than an exception: Aligner.__bool__ reports
            # whether the index pointer is non-NULL, so `except` catches nothing here.
            if not self._aligner:
                raise RuntimeError(
                    "mappy could not load a transcript index from "
                    f"{transcript_fa}; streaming transcriptome rescue cannot run"
                )
            # An index built with --idx-no-seq answers every query with zero hits and no
            # error at all, because mappy always requests base-level alignment and
            # map() returns immediately when the index holds no sequences.  Aligner.seq()
            # is the only way to tell that index from a working one.
            probe_name = self._aligner.seq_names[0]
            if not self._aligner.seq(probe_name, 0, 1):
                raise RuntimeError(
                    "mappy loaded a transcript index that retains no sequences, which "
                    "yields zero hits with no error; streaming transcriptome rescue "
                    "cannot run against it"
                )
            # A multi-part index is silently truncated to its first part -- mappy.pyx
            # reads only that one, with no warning -- so a read would be offered fewer
            # targets than it is entitled to and simply fail to rescue.  Indexes built
            # from a FASTA at the default -I 8G are single-part, and this is what
            # establishes that rather than assuming it.
            if self._aligner.n_seq != len(self._transcript_models):
                raise RuntimeError(
                    "mappy sees {} of {} transcript targets, so the index is truncated "
                    "(a multi-part index exposes only its first part); streaming "
                    "transcriptome rescue would offer reads an incomplete target "
                    "set".format(self._aligner.n_seq, len(self._transcript_models))
                )
        finally:
            # The index is built and resident by now; mappy does not read the file
            # again, and leaving it would accumulate one FASTA per contig-strand.
            try:
                os.unlink(transcript_fa)
                os.rmdir(os.path.dirname(transcript_fa))
            except OSError:
                pass

        # One scratch buffer for the life of the worker.  Letting map() allocate its own
        # pays an mm_tbuf_init/destroy pair per read, which at rescue volumes is the
        # difference between a fixed setup cost and a per-read one.
        self._buf = mappy.ThreadBuffer()

        target_ids = list(self._transcript_models.keys())
        self._header = pysam.AlignmentHeader.from_references(
            target_ids,
            [len(self._transcript_models[t]["sequence"]) for t in target_ids],
        )
        self._target_id_to_tid = {t: i for i, t in enumerate(target_ids)}

        logger.info(
            "[%s%s] streaming rescue: resident mappy index over %d transcript targets",
            splice_graph.get_contig_acc(),
            splice_graph.get_contig_strand(),
            self._aligner.n_seq,
        )

    def _write_transcript_fasta(self, tmp_dir):
        # mappy indexes from a file: Aligner(seq=...) takes a single sequence and names
        # it "N/A", so it cannot express a multi-target index.  Written into the run's
        # per-contig-strand temp dir like the batch path's scratch, and deleted as soon
        # as the index is resident.
        parent = tmp_dir or os.environ.get("LRAA_TMP_DIR") or os.getcwd()
        if not os.path.isdir(parent):
            parent = os.getcwd()
        holder = tempfile.mkdtemp(
            prefix="stream_rescue_idx.{}.{}.".format(
                self._splice_graph.get_contig_acc(),
                self._splice_graph.get_contig_strand(),
            ),
            dir=parent,
        )
        transcript_fa = os.path.join(holder, "transcripts.fa")
        IsoformReadRescue._write_transcript_fasta(
            transcript_fa, self._transcript_models
        )
        return transcript_fa

    def offer(self, read, category):
        """Rescue paths for one candidate bam record, or () if it does not rescue.

        `read` is the genome-aligned record the streaming pass has in hand.  Its
        sequence is taken and reverse-complemented for a reverse-strand record, exactly
        as _collect_read_sequences does, so what reaches the aligner is the cDNA-oriented
        read the batch path would have written into its reads FASTA.
        """
        self.reads_offered[category] += 1
        read_name = Util_funcs.get_read_name_include_sc_encoding(read)
        if self._record_candidate_names:
            self.candidate_read_names[category].add(read_name)
        if self._aligner is None:
            return ()

        seq = read.query_sequence
        if not seq:
            # Only a primary record carries the full sequence, and the streaming pass
            # only ever holds primary records, so this is a bam without SEQ rather than
            # a filtering artifact.  The batch path drops these too.
            self.reads_without_sequence += 1
            return ()
        if read.is_reverse:
            seq = IsoformReadRescue._reverse_complement(seq)

        # Timed separately from acceptance so the report can say what a read costs to
        # ALIGN, which is the number that decides whether this scales, rather than what a
        # read costs to align plus project plus refine boundaries.
        started = time.perf_counter()
        alignments = self._pysam_records_for(read_name, seq)
        self.align_seconds += time.perf_counter() - started
        if not alignments:
            return ()
        self.hits_examined += len(alignments)
        accept_started = time.perf_counter()

        # The genome baseline this rescue has to beat, computed from the very record in
        # hand rather than re-fetched: an alignment explaining less of the read than the
        # genome already does cannot correct anything.
        genome_explained = IsoformReadRescue._explained_read_bases(read)
        genome_gap_id = IsoformReadRescue._gap_aware_identity(read)

        rescued_mps, _details = IsoformReadRescue._parse_rescue_alignments(
            alignments,
            self._splice_graph,
            self._transcript_models,
            read_path_mapper=self._read_path_mapper,
            require_unique_path_across_best_hits=True,
            split_multipaths_by_gene=False,
            read_name_to_genome_explained=(
                {} if genome_explained is None else {read_name: genome_explained}
            ),
            read_name_to_genome_gap_id=(
                {} if genome_gap_id is None else {read_name: genome_gap_id}
            ),
        )
        self.accept_seconds += time.perf_counter() - accept_started
        if not rescued_mps:
            return ()
        self.reads_rescued[category] += 1
        return tuple(tuple(mp.get_simple_path()) for mp in rescued_mps)

    def _pysam_records_for(self, read_name, seq):
        """One pysam record per mappy hit, flagged as minimap2's SAM writer would flag it.

        minimap2 labels hits in the order mm_map returns them: a hit whose id differs
        from its parent is secondary; among the rest the FIRST is the SAM primary and
        every later one is supplementary (mm_set_sam_pri).  mappy preserves that order
        but exposes only is_primary, so the primary/supplementary split has to be
        rebuilt here -- and it matters, because the batch acceptance drops supplementary
        hits and keeps secondary ones.

        The query name is NOT passed to map().  mappy 2.28's map() has no `name`
        parameter -- that arrived in a later release -- and it would be inert here in any
        case: mm_map uses the query name only for all-vs-all self-hit suppression, which
        none of the map-* presets enable.
        """
        records = []
        primaries_seen = 0
        for hit in self._aligner.map(seq, buf=self._buf):
            if hit.ctg not in self._target_id_to_tid:
                continue
            flag = 0
            if hit.strand < 0:
                flag |= SAM_FLAG_REVERSE
            if hit.is_primary:
                primaries_seen += 1
                if primaries_seen > 1:
                    flag |= SAM_FLAG_SUPPLEMENTARY
            else:
                flag |= SAM_FLAG_SECONDARY

            # Soft clips in SAM orientation.  For a reverse hit the stored query is the
            # reverse complement of what mappy was given, while q_st/q_en stay on the
            # original strand (PAF convention), so the two clip lengths swap.
            left_clip, right_clip = hit.q_st, len(seq) - hit.q_en
            if hit.strand < 0:
                left_clip, right_clip = right_clip, left_clip

            cigar = []
            if left_clip:
                cigar.append((BAM_CSOFT_CLIP, left_clip))
            cigar.extend((op, length) for length, op in hit.cigar)
            if right_clip:
                cigar.append((BAM_CSOFT_CLIP, right_clip))

            record = pysam.AlignedSegment(self._header)
            record.query_name = read_name
            record.flag = flag
            record.reference_id = self._target_id_to_tid[hit.ctg]
            record.reference_start = hit.r_st
            record.mapping_quality = hit.mapq
            record.cigartuples = cigar
            # No SEQ: acceptance reads read length off the CIGAR (infer_read_length), so
            # carrying the sequence on every hit would only copy it once per hit.
            record.set_tag("NM", int(hit.NM), value_type="i")
            records.append(record)
        return records

    def close(self):
        self._aligner = None
        self._buf = None

    def stats(self):
        return {
            "reads_offered": dict(self.reads_offered),
            "reads_rescued": dict(self.reads_rescued),
            "reads_without_sequence": self.reads_without_sequence,
            "hits_examined": self.hits_examined,
            "n_targets": 0 if self._aligner is None else self._aligner.n_seq,
        }
