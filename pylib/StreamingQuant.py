#!/usr/bin/env python3

"""Assign every read in a bam without holding any of them.

The default quantification path materializes each shard's alignments, collapses them
into multipaths, runs EM, and emits tracking from the retained read associations. That
is fine at millions of reads and impossible at billions: the alignment objects alone
ran to 1,266 MB for 173k reads on a chr20 HiFi shard, and the read names needed for
tracking are themselves ~90 bytes apiece in Python.

This module implements the second half of a two-pass alternative. A first pass
quantifies normally against a coverage-normalized bam -- whose read count is bounded by
the depth target far more tightly than by library size -- and hands over the abundances
it settled on. This pass then streams the full bam and, for each read, looks up the
answer that pass already computed, writes the tracking row, and forgets the read.

Two properties measured on chr20 make it work:

  - a read's path through the splice graph is a function of the read and the graph
    alone, identical across runs (86,655/86,655 reads), so this pass can recompute it
    without anything carried over;
  - 99.97% of full-bam reads land on a path the first pass already resolved, provided
    the table also records the paths that matched *nothing*. Recording only successes
    drops that to 89%, because ~28% of paths have no compatible isoform and would be
    re-tested for every read carrying them.

Resident memory here is the table plus per-isoform accumulators -- thousands of
entries -- and is independent of how many reads stream past.
"""

import gzip
import logging
import os
from collections import defaultdict

import pysam

import LRAA_Globals
import MultiPath
import Util_funcs
from LRAA_Globals import SPACER
from Pretty_alignment import Pretty_alignment

logger = logging.getLogger(__name__)


class AssignmentTable:
    """Canonical path -> the per-isoform split a completed quantification settled on.

    A lookup returning None means the path was never seen, and the caller has to resolve
    it. A lookup returning an empty list means it was seen and matched no isoform, which
    is a real answer and must be cached as one.
    """

    __slots__ = ("_rows",)

    def __init__(self):
        self._rows = dict()

    def __len__(self):
        return len(self._rows)

    @property
    def num_negative(self):
        return sum(1 for v in self._rows.values() if not v)

    def lookup(self, canonical_path):
        return self._rows.get(canonical_path)

    def add(self, canonical_path, rows):
        self._rows[canonical_path] = rows

    @classmethod
    def build(cls, splice_graph, mp_to_transcripts, theta, em_was_run,
              unassigned_mp_count_pairs, transcript_id_to_component_id):
        """Collapse a finished quantification into the rows a streaming pass emits.

        Fractions are recomputed from the converged theta rather than reused from the ones EM
        returned. Those two now AGREE: EM runs a final E-step at the theta it returns
        (EM.py, after the iteration loop), so its fractions and counts describe the same
        abundances it hands back. Recomputing here is therefore equivalent for a path the
        first pass saw, and necessary for one it did not -- the resolver has no returned
        fraction to reuse and must derive it from theta.

        Before that fix the two differed, unboundedly on genes that exhausted max_iter, and
        this recomputation existed to keep one semantics in the file rather than two. Keeping
        it means seen and resolved paths still run through a single expression, which is what
        makes them comparable at all.

        The denominator is per read-sharing COMPONENT, because that is the unit EM runs on
        and therefore the unit theta is normalized within. A path compatible with isoforms
        of two genes in one component gets ONE split across all of them; per-gene groups
        would each sum to 1 and count the read twice.

        With no theta -- run_EM False -- the split is equal across a gene's compatible
        isoforms, which is what _estimate_isoform_read_support does in that mode. Falling
        through to a theta-weighted formula would divide by zero and silently assign
        nothing.

        Everything a tracking row needs except the read name is fixed per (path, isoform),
        so it is computed once here rather than per read.
        """
        table = cls()
        use_3p = LRAA_Globals.config["weight_reads_by_3prime_agreement"]
        rows = defaultdict(list)

        # Told, not inferred, in both directions. An empty theta is ambiguous on its own:
        # run_EM=False, or EM ran and the state was reset between the quantification and this
        # call. Inferring from emptiness takes the equal-split branch in the second case and
        # emits a table that is complete, internally consistent, and built from no abundance
        # information at all.
        if em_was_run and not theta:
            raise RuntimeError(
                "EM ran but no theta was supplied: the quantification state this table "
                "would be built from does not belong to that EM run"
            )
        # The converse is equally wrong and quieter: theta present with EM off would take the
        # weighted branch while claiming to reproduce the no-EM equal split.
        if theta and not em_was_run:
            raise RuntimeError(
                "theta was supplied but EM did not run: the no-EM path assigns an equal "
                "split, so weighting by these abundances would not reproduce it"
            )
        have_theta = bool(theta)

        # No default when a transcript is absent. Every transcript this quantification
        # assigned has a component id, so a missing one means the component map does not
        # belong to the run that produced these assignments -- and a fabricated id would
        # silently give that transcript a denominator of its own.
        def _component_of(transcript):
            tid = transcript.get_transcript_id()
            try:
                return transcript_id_to_component_id[tid]
            except KeyError:
                raise RuntimeError(
                    "no component id for {}: this component map does not belong to the "
                    "quantification that produced these assignments".format(tid)
                )

        for mp, transcripts_assigned in mp_to_transcripts.items():
            canon = splice_graph.canonical_simple_path(mp.get_simple_path())
            # One row set per path. Two multipaths reaching one canonical path would append
            # both sets under this key, and a streamed read on that path would then be
            # written once per set -- its fractions summing past 1.0 with nothing out of
            # place in the file. Raising rather than merging because the correct merge is
            # undefined: the two multipaths carry different weights and different mp ids.
            if canon in rows:
                raise RuntimeError(
                    "two multipaths map to canonical path {}: cannot build one row set "
                    "for it".format(canon)
                )
            rows[canon] = rows_for_multipath(
                mp,
                transcripts_assigned,
                theta,
                em_was_run,
                use_3p,
                lambda t, _mp=mp: float(t.get_multipath_weight(_mp)),
                _component_of,
            )

        for canon, entries in rows.items():
            table.add(canon, entries)

        # Paths that were evaluated and matched nothing. Caching the negative is what keeps
        # a streaming pass from re-resolving them read after read.
        negatives = 0
        for mp_count_pair in unassigned_mp_count_pairs:
            mp, _count = mp_count_pair.get_multipath_and_count()
            canon = splice_graph.canonical_simple_path(mp.get_simple_path())
            # An unassigned multipath sharing a canonical path with an assigned one is the
            # same collision the positive side refuses, and skipping it quietly would be a
            # decision rather than a no-op: default mode dropped this multipath's reads,
            # while a table keeping the positive entry assigns them its split. Refused in
            # both places or the two sides disagree about whether canonicalization is
            # injective.
            if canon in table._rows:
                raise RuntimeError(
                    "unassigned multipath shares canonical path {} with an assigned one: "
                    "cannot record it as matching nothing".format(canon)
                )
            table.add(canon, [])
            negatives += 1

        logger.info(
            "streaming assignment table: %d paths (%d with isoforms, %d recorded as "
            "matching none); theta %s",
            len(table), len(table) - table.num_negative, negatives,
            "from EM" if have_theta else "absent, equal splits",
        )
        return table


def rows_for_multipath(mp, transcripts_assigned, theta, em_was_run, use_3p, weight_of,
                       component_of):
    """The tracking rows one multipath contributes: (gene, tx, hash, exons, mp, frac, w).

    The single definition of the split, called both when precomputing the table from a
    finished quantification and when a streaming pass resolves a path that quantification
    never saw. Two copies would agree the day they were written and then drift, and the
    drift would appear as reads assigned differently depending on whether their path
    happened to survive coverage normalization -- indistinguishable, in the output, from
    real biology.

    `weight_of(transcript)` supplies the 3'-agreement weight, because its two callers hold
    it differently: the table reads the value stored on the transcript during assignment,
    while the resolver computes it for a path that has none. Both values come from
    Quantify._assign_read_weights_based_on_read_end_agreement, which is a function of path
    geometry alone, so they agree by construction rather than by coincidence.

    `component_of(transcript)` supplies the read-sharing component the transcript was
    quantified in. The denominator is per COMPONENT, not per gene: EM runs once per component
    of genes joined by shared reads, so theta is normalized within a component and a per-gene
    denominator would renormalize a subset of it. For a read compatible with isoforms of two
    genes in one component that mistake is not subtle -- each gene's group sums to 1 and the
    read is counted twice.

    With no theta -- run_EM False -- the split is equal across the component's compatible
    isoforms, which is what _estimate_isoform_read_support does in that mode, called once per
    component.
    """
    by_component = defaultdict(list)
    gene_of_tid = {}
    component_of_tid = {}
    for t in transcripts_assigned:
        tid = t.get_transcript_id()
        gene_id = t.get_gene_id()
        # Tripwire for the stated invariant that isoforms are unique per gene. Independent of
        # grouping now that grouping is by component, but still worth catching: Quantify
        # validates this at indexing time, and this catches a caller that assembled the list
        # some other way.
        prior = gene_of_tid.setdefault(tid, gene_id)
        if prior != gene_id:
            raise RuntimeError(
                "transcript {} appears under two gene ids ({}, {}): isoforms must be "
                "unique per gene".format(tid, prior, gene_id)
            )
        component_id = component_of(t)
        prior_component = component_of_tid.setdefault(tid, component_id)
        if prior_component != component_id:
            raise RuntimeError(
                "transcript {} appears under two component ids ({}, {}): the split "
                "denominator would depend on which copy was read".format(
                    tid, prior_component, component_id
                )
            )
        by_component[component_id].append(t)

    # A path compatible with isoforms from more than one component cannot be split at all.
    #
    # Pass 1 cannot produce this: a multipath compatible with two genes' isoforms is exactly
    # what joins those genes into one component, so anything it saw is single-component by
    # construction. It arises when the streaming pass resolves a path the first pass never
    # saw, because coverage normalization sampled away the only read bridging two genes.
    # Those components then carry independently normalized thetas: pooling them compares
    # numbers from two different normalizations, and splitting each separately sums the read's
    # fractions to the number of components.
    #
    # Neither is an approximation worth making, so the run is refused. The normalized bam is
    # not a sufficient foundation for this locus, and the unseen-path-rate guard cannot speak
    # to it -- that bounds how much resolution happened, not whether a split was meaningful.
    if len(by_component) > 1:
        raise RuntimeError(
            "path {} is compatible with isoforms from {} separately quantified components "
            "({}). Their abundances were normalized independently, so this read cannot be "
            "split across them: pooling mixes two normalizations, and splitting per "
            "component assigns the read once per component. Coverage normalization removed "
            "the read that would have joined these genes into one component, so the "
            "normalized bam is an insufficient foundation here -- raise "
            "--normalize_max_cov_level, or run without --stream_reads.".format(
                mp.get_simple_path(),
                len(by_component),
                ",".join(sorted(str(c) for c in by_component)),
            )
        )

    out = []
    for group in by_component.values():
        # Everything below derives from this mapping rather than the list, so a repeat
        # cannot emit a row twice or skew a denominator by list length.
        by_tid = {}
        for t in group:
            by_tid[t.get_transcript_id()] = t
        if len(by_tid) != len(group):
            raise RuntimeError(
                "duplicate transcript ids in one component's candidate group: {}".format(
                    sorted(t.get_transcript_id() for t in group)
                )
            )

        weights = {
            tid: (float(weight_of(t)) if use_3p else 1.0) for tid, t in by_tid.items()
        }
        if theta:
            # No default. When EM ran, every transcript it quantified has a theta, so a
            # missing key means this transcript was not part of the run that produced it.
            # Defaulting to 0.0 would emit a well-formed row assigning the read nothing,
            # which no consumer can tell from a real zero-abundance isoform.
            missing = [tid for tid in by_tid if tid not in theta]
            if missing:
                raise RuntimeError(
                    "no EM theta for {}: assignment state does not belong to the "
                    "quantification that produced this theta".format(
                        ",".join(sorted(missing))
                    )
                )
            th = {tid: float(theta[tid]) for tid in by_tid}
            denom = sum(weights[k] * th[k] for k in weights)
            # Zero denominator yields zero fractions, so the read's fractions sum to 0
            # rather than 1. Deliberate: it is what the default path's E-step does
            # (EM.py:309-313), and reproducing that is the requirement. An equal split here
            # would be better arithmetic and a divergence. Unreached on ONT and PacBio
            # chr20 -- 0 zero-valued fractions in 719k rows.
            fracs = {
                k: (weights[k] * th[k] / denom if denom > 0 else 0.0) for k in weights
            }
        else:
            if em_was_run:
                raise RuntimeError(
                    "EM ran but no theta was supplied for this multipath's transcripts"
                )
            # Equal split over the gene's isoforms, matching the no-EM branch of
            # _estimate_isoform_read_support, which divides by the number of compatible
            # transcripts within the gene it is called for (Quantify.py:162 runs it per
            # gene) and applies no 3' weighting to the fraction.
            fracs = {k: 1.0 / len(by_tid) for k in weights}

        for tid, t in by_tid.items():
            num_exons = t.get_num_exon_segments()
            out.append(
                (t.get_output_gene_id(), tid, _splice_hash_code(t, num_exons),
                 num_exons, mp.get_id(), fracs[tid], weights[tid])
            )

    # deterministic row order: tracking output must not follow dict order
    out.sort(key=lambda r: (r[1], r[0]))
    return out


def _splice_hash_code(transcript, num_exons):
    """The tracking file's transcript_splice_hash_code, as report_quant_results makes it."""
    if num_exons > 1:
        return Util_funcs.get_hash_code(transcript.get_introns_string())
    return transcript.get_transcript_id()


def make_path_resolver(splice_graph, quantify_obj, theta, em_was_run,
                       transcript_id_to_component_id,
                       fraction_read_align_overlap=None):
    """Resolve a splice-graph path the first pass never saw, as that pass would have.

    Reads surviving coverage normalization do not cover every path in the full bam, so a
    streaming pass meets paths absent from the precomputed table. Dropping them would lose
    reads; guessing would invent numbers. This runs the same anchoring, the same
    compatibility cascade and the same 3'-weighting the first pass ran, against the same
    splice graph and the same transcript set, then splits by the same theta.

    Every step is a call into Quantify rather than a reimplementation. The cascade in
    particular is ten tests whose ORDER decides which isoform a read lands on, so a second
    copy of it would not fail loudly -- it would quietly assign some reads differently.

    Returns the same 7-tuple rows as the table holds, or [] for a path that anchors to no
    gene or matches no isoform. [] is a real answer and is cached as one: it is what stops
    a streaming pass from re-resolving an unmatchable path once per read.
    """
    if fraction_read_align_overlap is None:
        fraction_read_align_overlap = LRAA_Globals.config["fraction_read_align_overlap"]
    use_3p = LRAA_Globals.config["weight_reads_by_3prime_agreement"]

    # Same no-default rule as the table's: an unseen path's candidate isoforms all come from
    # the first pass, so each has a component id. A fabricated one would hand a transcript a
    # denominator of its own and hide the cross-component case this must refuse.
    def _component_of(transcript):
        tid = transcript.get_transcript_id()
        try:
            return transcript_id_to_component_id[tid]
        except KeyError:
            raise RuntimeError(
                "no component id for {}: this component map does not belong to the "
                "quantification whose transcripts this resolver anchors against".format(tid)
            )

    def resolve(simple_path):
        # A single-read multipath standing for this path. Read count 1: the resolver
        # answers "how does one read on this path split", and the caller applies it per
        # read, so a count here would double-apply.
        mp = MultiPath.MultiPath(splice_graph, [simple_path], read_count=1)

        top_genes = quantify_obj._get_all_genes_with_node_matches_to_simplepath(
            simple_path
        )
        if top_genes is None:
            return []

        gene_isoforms = quantify_obj.candidate_isoforms_for_genes(top_genes)
        if not gene_isoforms:
            return []

        transcripts_assigned = quantify_obj.resolve_mp_to_transcripts(
            splice_graph, mp, gene_isoforms, fraction_read_align_overlap
        )
        if transcripts_assigned is None:
            return []

        # Computed rather than read off the transcripts: these objects carry weights for the
        # multipaths of the first pass, and this path is not one of them. The function is a
        # function of path geometry alone, so the value it returns here is the value the
        # first pass would have stored had it seen this path.
        weight_by_tid = quantify_obj._assign_read_weights_based_on_read_end_agreement(
            splice_graph, mp, transcripts_assigned
        )

        return rows_for_multipath(
            mp,
            transcripts_assigned,
            theta,
            em_was_run,
            use_3p,
            lambda t: weight_by_tid[t.get_transcript_id()],
            _component_of,
        )

    return resolve


class StreamingTotals:
    """Per-isoform and per-gene running sums. Sized by the annotation, not the library."""

    def __init__(self):
        self.frac_sum = defaultdict(float)
        self.uniq_reads = defaultdict(int)
        self.gene_frac_sum = defaultdict(float)
        self.reads_streamed = 0
        self.reads_filtered = 0
        self.reads_no_path = 0
        self.reads_spacer_path = 0
        self.reads_unassignable = 0
        self.reads_assigned = 0
        # paths absent from the table and resolved here instead. Counted rather than
        # dropped: every one of their reads is still assigned, and this is the measure
        # of how much work the first pass failed to precompute.
        self.paths_resolved_in_stream = 0
        # reads landing on those paths. The guard divides by reads, so it has to count
        # reads: incrementing only per distinct path made the threshold unreachable --
        # 10k novel paths over 100M reads reads as 0.01%, and a first pass that
        # precomputed nothing would never trip it.
        self.reads_on_stream_resolved_paths = 0
        # Reads matched to isoforms whose fractions all came out zero, which happens when
        # every compatible isoform of the gene has theta 0. Those reads are counted as
        # assigned and contribute nothing, so without this they are invisible: the tracking
        # file holds well-formed rows reading 0.000000 and the read's mass is simply absent
        # from every count. Zero on ONT and PacBio chr20, but reachable wherever
        # normalization thins a rare isoform's evidence to nothing while the full bam still
        # carries reads that only it explains -- which is the regime this mode is for.
        self.reads_zero_fraction = 0

    def record(self, rows):
        report_min = LRAA_Globals.config["unique_read_report_min_frac"]
        any_positive = False
        for gene_id, transcript_id, _hash, _nex, _mp, frac, _w in rows:
            self.frac_sum[transcript_id] += frac
            self.gene_frac_sum[gene_id] += frac
            if frac > 0.0:
                any_positive = True
            if frac >= report_min:
                self.uniq_reads[transcript_id] += 1
        if rows and not any_positive:
            self.reads_zero_fraction += 1


def streaming_read_path(
    read, lraa_obj, splice_graph, contig_seq, try_correct_alignments=True
):
    """The splice-graph path one bam record maps to, computed without retaining it.

    Must agree with what the default path computes for the same record. The default flow
    is Bam_alignment_extractor -> Pretty_alignment_manager (try_correct_alignments, then
    prune_long_terminal_introns, then lighten) -> LRAA._map_read_to_graph, so this applies
    the same operations in the same order. Omitting either would diverge on exactly the
    reads those steps touch, which is a small set and therefore easy to miss.

    Pruning uses the single-alignment primitive, and correction is invoked with quiet=True,
    because both batch helpers log once per call and construct a progress bar per call: at a
    billion reads with ~1.5% soft-clip candidates that is tens of millions of log lines and
    progress objects. The correction arithmetic itself is per-alignment and unaffected.

    A named function rather than an inline block so the parity check exercises the code
    that actually runs, not a copy of it. Verified at 104,402/104,402 reads identical to
    the default path on ONT chr20.
    """
    pretty_alignment = Pretty_alignment.get_pretty_alignment(read)
    if try_correct_alignments and pretty_alignment.is_softclip_realign_candidate():
        Pretty_alignment.try_correct_alignments(
            [pretty_alignment], splice_graph, contig_seq, quiet=True
        )
    pretty_alignment.prune_long_terminal_introns_single()
    pretty_alignment.lighten()
    return lraa_obj._map_read_to_graph(
        pretty_alignment.get_pretty_alignment_segments(),
        snap_nearby_boundary_features=True,
        left_soft_clipping=pretty_alignment.left_soft_clipping,
        right_soft_clipping=pretty_alignment.right_soft_clipping,
    )


def dump_streaming_paths(
    bam_file,
    contig_acc,
    contig_strand,
    contig_seq,
    lraa_obj,
    out_path,
    try_correct_alignments=True,
    region_lend=None,
    region_rend=None,
):
    """Write read_name -> canonical path as the streaming pass computes it.

    Diagnostic. Run immediately after the default pass over the same bam, against the very
    same splice graph and LRAA object, so a diff against the default pass's own dump
    isolates path computation from graph-construction variance. Comparing dumps from two
    separate invocations would not: node ids are process-global counters, and two
    independently built graphs may differ.

    This exercises streaming_read_path, the function stream_assign itself calls, so the
    parity number describes production code rather than a copy of it.
    """
    splice_graph = lraa_obj._splice_graph
    n = 0
    with pysam.AlignmentFile(bam_file, "rb") as reader:
        if region_lend is not None and region_rend is not None:
            fetcher = reader.fetch(contig_acc, region_lend, region_rend)
        else:
            fetcher = reader.fetch(contig_acc)
        with open(out_path, "wt") as fh:
            for read in fetcher:
                # Runs inside run_quant_only's mapping-quality swap, so config already
                # holds min_mapping_quality_for_final_quant -- the value the first pass
                # extracted under. Explicit thresholds must not be hoisted here: they
                # would pin the discovery value and diverge from the pass this
                # reproduces.
                if Util_funcs.quant_discard_reason(read, contig_strand) is not None:
                    continue
                path = streaming_read_path(
                    read, lraa_obj, splice_graph, contig_seq, try_correct_alignments
                )
                if not path or path == SPACER or SPACER in path:
                    continue
                fh.write(
                    "{}\t{}\n".format(
                        Util_funcs.get_read_name_include_sc_encoding(read),
                        splice_graph.canonical_simple_path(path),
                    )
                )
                n += 1
    logger.info(
        "[%s%s] wrote %d streaming read paths to %s", contig_acc, contig_strand, n, out_path
    )
    return n


def stream_assign(
    bam_file,
    contig_acc,
    contig_strand,
    contig_seq,
    lraa_obj,
    table,
    tracking_fh,
    resolver=None,
    try_correct_alignments=True,
    region_lend=None,
    region_rend=None,
    max_fallback_path_frac=None,
):
    """Stream the bam, emit one tracking row per (read, compatible isoform).

    Holds no per-read state: each read is mapped, looked up, written, and dropped.

    A path the table has never seen is resolved once by `resolver` and the verdict is
    cached, so the read at hand and every later one on that path are still assigned. No
    read is dropped for being unfamiliar -- completeness of the tracking file is why this
    pass exists. `max_fallback_path_frac` bounds how much fallback resolution is
    tolerable before the run is refused, since a table built from too thin a normalized
    bam would degrade into re-resolving everything.

    One record per read is a structural guarantee rather than a check:
    Util_funcs.quant_discard_reason rejects secondary and supplementary alignments,
    so there is no grouping for a one-record-at-a-time pass to reconstruct.
    """
    totals = StreamingTotals()
    splice_graph = lraa_obj._splice_graph

    # Paths resolved during this stream. Bounded by the number of distinct novel paths,
    # not by reads, and needed so a path's later reads are attributed to the fallback
    # rather than counted as table hits -- the guard divides by reads, so it has to see
    # every read that depended on in-stream resolution, not just the first.
    stream_resolved = set()

    with pysam.AlignmentFile(bam_file, "rb") as reader:
        if region_lend is not None and region_rend is not None:
            fetcher = reader.fetch(contig_acc, region_lend, region_rend)
        else:
            fetcher = reader.fetch(contig_acc)

        for read in fetcher:
            # Exactly the filter the first pass applied, shared rather than copied.
            # Config is read rather than passed: this runs inside run_quant_only's
            # mapping-quality swap, so it already holds
            # min_mapping_quality_for_final_quant. Pinning explicit thresholds here
            # would reintroduce the discovery value and drop or admit reads the pass
            # this stream extends did not.
            if Util_funcs.quant_discard_reason(read, contig_strand) is not None:
                totals.reads_filtered += 1
                continue

            totals.reads_streamed += 1

            path = streaming_read_path(
                read, lraa_obj, splice_graph, contig_seq, try_correct_alignments
            )
            if not path or path == SPACER:
                totals.reads_no_path += 1
                continue
            if SPACER in path:
                totals.reads_spacer_path += 1
                continue

            canon = splice_graph.canonical_simple_path(path)
            rows = table.lookup(canon)
            if canon in stream_resolved:
                totals.reads_on_stream_resolved_paths += 1
            if rows is None:
                # Never seen. Resolve once with the ordinary compatibility cascade and
                # cache the verdict, positive or negative, so this read and every later
                # one on this path is assigned rather than silently omitted.
                totals.paths_resolved_in_stream += 1
                totals.reads_on_stream_resolved_paths += 1
                stream_resolved.add(canon)
                rows = resolver(path) if resolver is not None else None
                if rows is None:
                    raise RuntimeError(
                        "no resolver supplied for a path absent from the assignment "
                        "table; streaming can neither assign nor drop it: " + str(canon)
                    )
                table.add(canon, rows)
            if not rows:
                totals.reads_unassignable += 1
                continue

            read_name = Util_funcs.get_read_name_include_sc_encoding(read)
            for gene_id, transcript_id, splice_hash, num_exons, mp_id, frac, weight in rows:
                tracking_fh.write(
                    "{}\t{}\t{}\t{}\t{}\t{}\t{:.6f}\t{:.6f}\n".format(
                        gene_id, transcript_id, splice_hash, num_exons, mp_id,
                        read_name, frac, weight,
                    )
                )
            totals.record(rows)
            totals.reads_assigned += 1

    _report(contig_acc, contig_strand, totals, max_fallback_path_frac)
    return totals


def write_expr(transcripts, totals, ofh, quant_only=True):
    """Emit quant.expr from the streamed totals rather than from EM's own counts.

    The counts here are sums of per-read fractional assignments over the whole bam, so
    they are consistent with the tracking file this pass wrote. They are NOT the counts
    the default path reports: this is one expectation step at the abundances the first
    pass estimated, where the default path would go on to re-estimate them. How close
    the two land is a measured quantity, not an identity. Measured on ONT chr20: gene
    totals preserved, and 22 of 1785 expressed transcripts more than 5% apart at the
    isoform level, holding 8% of read mass.

    TPM is normalized over the transcripts passed in, which is this contig/strand's set.
    That matches the default path, whose per-shard quant.expr also sums to 1e6 and which
    the merge renormalizes genome-wide; measured identical on both modes.
    """
    total_reported = sum(totals.frac_sum.get(t.get_transcript_id(), 0.0) for t in transcripts)

    # Checked once, before any row is written, rather than per transcript: an empty
    # transcript list would never reach a per-row check, and raising midway would leave a
    # complete tracking file beside a truncated expr file, which a resuming pipeline could
    # read as success.
    num_total_reads = LRAA_Globals.config.get("num_total_reads")
    if not num_total_reads:
        raise RuntimeError(
            "num_total_reads is unset, so RPM cannot be computed; a zero-filled column is "
            "indistinguishable from a real measurement"
        )

    for transcript in sorted(transcripts, key=lambda t: t.get_transcript_id()):
        transcript_id = transcript.get_transcript_id()
        gene_id = transcript.get_output_gene_id()
        num_exons = transcript.get_num_exon_segments()
        counts = totals.frac_sum.get(transcript_id, 0.0)
        uniq = totals.uniq_reads.get(transcript_id, 0)
        gene_total = totals.gene_frac_sum.get(gene_id, 0.0)
        iso_frac = counts / gene_total if gene_total > 0 else 0.0
        uniq_gene_frac = uniq / gene_total if gene_total > 0 else 0.0
        tpm = counts / total_reported * 1e6 if total_reported > 0 else 0.0
        rpm = counts / num_total_reads * 1e6

        vals = [
            gene_id,
            transcript_id,
            f"{uniq}",
            f"{counts:.1f}",
            f"{iso_frac:.3f}",
            f"{uniq_gene_frac:.3f}",
            f"{tpm:.3f}",
            transcript.get_exons_string(),
            (transcript.get_introns_string() if num_exons > 1 else ""),
            _splice_hash_code(transcript, num_exons),
        ]
        if not quant_only:
            vals += ["", ""]
        vals.append(f"{rpm:.3f}")
        print("\t".join(vals), file=ofh)


def _report(contig_acc, contig_strand, totals, max_fallback_path_frac):
    tag = f"[{contig_acc}{contig_strand}]"
    logger.info(
        "%s streamed %d reads: %d assigned, %d matched no isoform, %d no graph path, "
        "%d spacer path; %d paths (%d reads) resolved here rather than by the first pass",
        tag, totals.reads_streamed, totals.reads_assigned, totals.reads_unassignable,
        totals.reads_no_path, totals.reads_spacer_path, totals.paths_resolved_in_stream,
        totals.reads_on_stream_resolved_paths,
    )
    if totals.reads_zero_fraction:
        # Logged at WARNING because nothing in the output shows it: these reads have
        # tracking rows reading 0.000000 and contribute no mass to any count, so a reader
        # comparing counts against read totals would find the shortfall and have no way to
        # attribute it. Means the first pass gave every isoform compatible with these reads
        # an abundance of zero -- normalization thinned the evidence away while the full bam
        # still carries reads only those isoforms explain.
        logger.warning(
            "%s %d reads matched isoforms whose abundances were all zero, so they were "
            "assigned no mass; raise --normalize_max_cov_level if this is a large share",
            tag, totals.reads_zero_fraction,
        )
    accounted = totals.reads_assigned + totals.reads_unassignable
    if not accounted:
        return
    # Every read is assigned either way, so this bounds wasted work, not lost data: a table
    # built from too thin a normalized bam degrades into re-resolving everything, which is
    # slow rather than wrong.
    #
    # Checked after the stream, so it reports rather than prevents: the shard has already
    # done the work by the time this runs. Deliberate -- the rate is only trustworthy over
    # the whole shard, since coordinate order means early reads sample one locus, not the
    # library. It fails the run so a pipeline does not silently keep paying the cost.
    #
    # Measured over READS, not distinct paths. Counting paths made this unreachable: ten
    # thousand novel paths across a hundred million reads reads as 0.01%, so a first pass
    # that precomputed nothing would never have tripped it.
    frac_fallback = totals.reads_on_stream_resolved_paths / accounted
    logger.info(
        "%s reads needing in-stream path resolution: %.4f%% of accounted reads",
        tag, 100.0 * frac_fallback,
    )
    # `is not None`, so an explicit 0.0 means "refuse any fallback" rather than "disabled".
    # A falsy test silently turns the strictest possible setting into no setting at all.
    if max_fallback_path_frac is not None and frac_fallback > max_fallback_path_frac:
        raise RuntimeError(
            "{} {:.2%} of accounted reads landed on paths the first pass never saw, above "
            "the {:.2%} limit. Assignments are still complete -- those paths were resolved "
            "here -- but the first pass precomputed too little to be worth the second, so "
            "the compatibility work is being done twice. Raise "
            "--normalize_max_cov_level, or drop --stream_reads.".format(
                tag, frac_fallback, max_fallback_path_frac
            )
        )


def open_tracking_writer(path):
    """Tracking output, gzipped when the path says so.

    At a billion reads this file carries ~2.1 billion rows at ~140 bytes each, so
    whether it is compressed is the difference between roughly 0.3 TB and 0.05 TB.
    """
    if str(path).endswith(".gz"):
        return gzip.open(path, "wt")
    return open(path, "wt")
