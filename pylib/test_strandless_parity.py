#!/usr/bin/env python3

"""The equivalence gate for strandless chunking, held to its own contract.

Two things are asserted here, and they are different in kind.

The STRAND-ASSIGNMENT INVARIANT is a claim about the pipeline: a read's strand is
decided from that read alone, so a split run inside a chunk must agree with a
split run over the whole contig.  The fixture makes that claim non-trivial on
purpose -- reads whose splice dinucleotides CONTRADICT the aligner's flag, so the
splitter rewrites the flag and the invariant has something to be wrong about.  A
fixture with no flipped read would pass while testing nothing, so the flip count
is asserted rather than assumed.

The GATE ITSELF is a measuring instrument, and an instrument that cannot register
a fault is worthless.  So each check is also run against a deliberately broken
arrangement: a chunk whose annotation disagrees with the genome's, a
strand-filtered extraction from a raw BAM (the ordering mistake), a tracking table
whose rows were reordered against one whose values were changed.  The gate has to
tell those apart, and name what it found rather than count it.
"""

import gzip
import json
import os

import pysam
import pytest

import StrandlessParity as gate


CONTIG = "chrT"
LENGTH = 30000
BOUNDARY = 15000

# 1-based inclusive intron coordinates -> the dinucleotide pair to plant there.
# "GT..AG" reads as forward, "CT..AC" as its reverse complement, per
# separate_bam_by_strand.majority_vote_intron_orient.
FORWARD = ("GT", "AG")
REVERSE = ("CT", "AC")


class Corpus:
    """A contig whose splice sites and annotation are chosen, not incidental."""

    def __init__(self, tmp_path):
        self.dir = tmp_path
        self.sequence = list(("ACGT" * (LENGTH // 4 + 1))[:LENGTH])
        self.reads = []  # (name, flag, blocks)
        self.genes = []  # (gene_id, strand, exons)

    # -- construction -----------------------------------------------------------

    def spliced_read(self, name, flag, left, right, dinucs):
        """One read with exactly one intron, whose dinucleotides are planted."""

        intron_lend = left[1] + 1
        intron_rend = right[0] - 1
        self.sequence[intron_lend - 1] = dinucs[0][0]
        self.sequence[intron_lend] = dinucs[0][1]
        self.sequence[intron_rend - 2] = dinucs[1][0]
        self.sequence[intron_rend - 1] = dinucs[1][1]
        self.reads.append((name, flag, [left, right]))
        return self

    def unspliced_read(self, name, flag, block):
        self.reads.append((name, flag, [block]))
        return self

    def gene(self, gene_id, strand, exons):
        self.genes.append((gene_id, strand, exons))
        return self

    def build(self):
        self.fasta = str(self.dir / "genome.fa")
        sequence = "".join(self.sequence)
        with open(self.fasta, "wt") as ofh:
            print(">{}".format(CONTIG), file=ofh)
            for i in range(0, LENGTH, 60):
                print(sequence[i : i + 60], file=ofh)
        pysam.faidx(self.fasta)

        self.bam = str(self.dir / "reads.bam")
        header = {
            "HD": {"VN": "1.6", "SO": "coordinate"},
            "SQ": [{"SN": CONTIG, "LN": LENGTH}],
        }
        records = [self._alignment(*read) for read in self.reads]
        records.sort(key=lambda aln: aln.reference_start)
        with pysam.AlignmentFile(self.bam, "wb", header=header) as ofh:
            for record in records:
                ofh.write(record)
        pysam.index(self.bam)

        self.gtf = str(self.dir / "annot.gtf")
        with open(self.gtf, "wt") as ofh:
            for gene_id, strand, exons in self.genes:
                lend = min(e[0] for e in exons)
                rend = max(e[1] for e in exons)
                transcript_id = gene_id + ".1"
                attrs = 'gene_id "{}"; transcript_id "{}";'.format(
                    gene_id, transcript_id
                )
                for feature, flend, frend in (
                    ("gene", lend, rend),
                    ("transcript", lend, rend),
                ) + tuple(("exon",) + e for e in sorted(exons)):
                    print(
                        "\t".join(
                            (
                                CONTIG,
                                "test",
                                feature,
                                str(flend),
                                str(frend),
                                ".",
                                strand,
                                ".",
                                attrs,
                            )
                        ),
                        file=ofh,
                    )
        return self

    def _alignment(self, name, flag, blocks):
        aln = pysam.AlignedSegment()
        aln.query_name = name
        aln.flag = flag
        aln.reference_id = 0
        aln.reference_start = blocks[0][0] - 1
        aln.mapping_quality = 60
        cigar = []
        for i, (lend, rend) in enumerate(blocks):
            if i:
                cigar.append((3, lend - blocks[i - 1][1] - 1))
            cigar.append((0, rend - lend + 1))
        aln.cigar = cigar
        length = sum(rend - lend + 1 for lend, rend in blocks)
        aln.query_sequence = "A" * length
        aln.query_qualities = pysam.qualitystring_to_array("I" * length)
        aln.set_tag("NM", 0)
        return aln


@pytest.fixture
def corpus(tmp_path):
    """Six reads chosen so every branch of the strand decision is exercised.

    ``flip_to_plus`` and ``flip_to_minus`` are the ones that matter: their splice
    dinucleotides contradict the aligner, so the splitter REWRITES the flag.  They
    are what makes the invariant falsifiable, and what the ordering mistake
    misassigns.
    """

    return (
        Corpus(tmp_path)
        # aligner says reverse, the splice sites say forward
        .spliced_read("flip_to_plus", 16, (2000, 2200), (2600, 3000), FORWARD)
        # aligner says forward, the splice sites say reverse
        .spliced_read("flip_to_minus", 0, (5000, 5200), (5600, 6000), REVERSE)
        # aligner and splice sites agree; lands in the second chunk
        .spliced_read("agree_plus", 0, (20000, 20200), (20600, 21000), FORWARD)
        # unspliced: no splice evidence, decided by the annotation it overlaps
        .unspliced_read("annot_to_plus", 16, (2100, 2250))
        # unspliced and unannotated: nothing to infer from, keeps the aligner flag
        .unspliced_read("undecidable", 16, (10000, 10200))
        # straddles the boundary: dropped by both chunks, kept by the whole scope
        .unspliced_read("spans_boundary", 0, (BOUNDARY - 100, BOUNDARY + 100))
        .gene("geneP", "+", [(2000, 2300), (2700, 3000)])
        .gene("geneM", "-", [(5000, 5300), (5700, 6000)])
        .gene("geneQ", "+", [(20000, 20300), (20700, 21000)])
        .build()
    )


def intervals():
    return [
        gate.Interval(CONTIG, 1, BOUNDARY),
        gate.Interval(CONTIG, BOUNDARY + 1, LENGTH),
    ]


def run_invariant(corpus, work_dir, **kwargs):
    return gate.check_strand_invariant(
        corpus.bam,
        corpus.fasta,
        str(work_dir),
        intervals=kwargs.pop("intervals", intervals()),
        gtf=corpus.gtf,
        contig=CONTIG,
        **kwargs,
    )


def by_name(records):
    return {record["read_name"]: record for record in records}


# ------------------------------------------------- the fixture is what it claims


def test_the_splitter_decides_the_strands_the_fixture_was_built_for(corpus, tmp_path):
    """Pins the fixture. Every later assertion is meaningless if this drifts."""

    assignments, counters, _ = gate.split_assignments(
        corpus.bam,
        str(tmp_path / "whole"),
        200000,
        infer_read_orient=True,
        genome_fa=corpus.fasta,
        gtf=corpus.gtf,
    )
    strands = {name: strand for (name, _), strand in assignments.items()}
    assert strands["flip_to_plus"] == "+"  # dinucleotides beat the reverse flag
    assert strands["flip_to_minus"] == "-"  # dinucleotides beat the forward flag
    assert strands["agree_plus"] == "+"
    assert strands["annot_to_plus"] == "+"  # annotation beats the reverse flag
    assert strands["undecidable"] == "-"  # nothing to infer from, flag stands
    assert counters["num_records_strand_flipped"] == 3
    assert counters["num_inferred_by_splice_dinucs"] == 3
    assert counters["num_inferred_by_annot_overlap"] == 1
    assert counters["num_records_strand_uncertain"] == 2  # undecidable, spans_boundary


# --------------------------------------------------------------- the invariant


def test_invariant_holds_on_the_aligner_flag(corpus, tmp_path):
    """The mode the pipeline runs today: no inference, so no read can move."""

    report = run_invariant(corpus, tmp_path / "flag")

    assert report["mode"] == "aligner_flag"
    assert report["verdict"] == "INVARIANT HOLDS"
    assert report["agreement"] == 1.0
    assert report["records_disagreed"] == 0
    assert report["whole_scope_records_strand_flipped"] == 0


def test_invariant_holds_for_inferred_orientation_despite_flips(corpus, tmp_path):
    """The forward-looking mode, and the only one where a flip exists at all.

    The flip assertions are the point: they prove the check was asked a question
    it could have failed.  Three reads have their flag rewritten genome-wide and
    the same three inside their chunks, and every compared read still agrees.
    """

    report = run_invariant(corpus, tmp_path / "infer", infer_read_orient=True)

    assert report["whole_scope_records_strand_flipped"] == 3
    assert report["in_chunk_records_strand_flipped"] == 3
    assert report["records_compared"] == 5  # six reads, one severed
    assert report["records_agreed"] == 5
    assert report["agreement"] == 1.0
    assert report["verdict"] == "INVARIANT HOLDS"


def test_a_read_the_chunks_never_saw_is_named_with_its_reason(corpus, tmp_path):
    """Not a count. The severed accounting differs between the two designs, so a
    count cannot say whether a difference is correct; a name and a reason can."""

    report = run_invariant(corpus, tmp_path / "severed", infer_read_orient=True)

    only_whole = by_name(report["records_only_in_whole_scope"])
    assert set(only_whole) == {"spans_boundary"}
    record = only_whole["spans_boundary"]
    assert record["reason"].startswith("overhang_dropped_at_")
    assert "chrT_00" in record["reason"] and "chrT_01" in record["reason"]
    assert record["whole_scope_strand"] == "+"
    # and the same read is named in the per-chunk drop accounting
    assert set(report["overhang_dropped_read_names"]) == {"spans_boundary"}


def test_a_chunk_whose_annotation_disagrees_is_caught_and_the_read_is_named(
    corpus, tmp_path
):
    """The instrument must register a fault.

    A chunk's GTF is the region's annotation rebased.  If per-chunk scoping ever
    dropped or altered a locus an in-chunk read overlaps, the annotation fallback
    would answer differently inside the chunk from outside it -- silently, because
    both answers are well-formed.  Flipping one gene's strand in one chunk's GTF
    stages exactly that, and the gate must name the read it moved.
    """

    work = tmp_path / "corrupt"
    run_invariant(corpus, work, infer_read_orient=True)

    chunk_gtf = os.path.join(str(work), "extracted", "chrT_00", "chunk.gtf")
    with open(chunk_gtf, "rt") as fh:
        lines = fh.readlines()
    with open(chunk_gtf, "wt") as fh:
        for line in lines:
            fields = line.rstrip("\n").split("\t")
            if "geneP" in fields[-1]:
                fields[6] = "-"
            fh.write("\t".join(fields) + "\n")

    report = gate.check_strand_invariant(
        corpus.bam,
        corpus.fasta,
        str(tmp_path / "corrupt_check"),
        chunks_root=os.path.join(str(work), "extracted"),
        gtf=corpus.gtf,
        contig=CONTIG,
        infer_read_orient=True,
    )

    assert report["verdict"] == "INVARIANT VIOLATED"
    assert report["records_disagreed"] == 1
    disagreement = report["disagreements"][0]
    assert disagreement["read_name"] == "annot_to_plus"
    assert disagreement["whole_scope_strand"] == "+"
    assert disagreement["in_chunk_strand"] == "-"
    # the raw flag is reported too: without it a reader cannot tell which side
    # rewrote the record
    assert disagreement["aligner_flag_strand"] == "-"
    assert report["agreement"] == 4 / 5


def test_the_invariant_refuses_a_stranded_region():
    """A strand-filtered extraction from a raw BAM filters on the flag the
    splitter has not rewritten yet. The check will not be pointed at one."""

    with pytest.raises(gate.ParityError, match="names a strand"):
        gate.intervals_from_region_strings(["chrT+:1-1000"])


def test_the_invariant_refuses_a_stranded_chunk_manifest(corpus, tmp_path):
    """Same prohibition, reached through a real run's chunk directory."""

    work = tmp_path / "stranded"
    run_invariant(corpus, work)
    manifest_path = os.path.join(
        str(work), "extracted", "chrT_00", "chunk.partition.json"
    )
    with open(manifest_path, "rt") as fh:
        manifest = json.load(fh)
    manifest["strand"] = "+"
    with open(manifest_path, "wt") as fh:
        json.dump(manifest, fh)

    with pytest.raises(gate.ParityError, match="already committed the ordering"):
        gate.check_strand_invariant(
            corpus.bam,
            corpus.fasta,
            str(tmp_path / "stranded_check"),
            chunks_root=os.path.join(str(work), "extracted"),
        )


def test_overlapping_intervals_are_refused(corpus, tmp_path):
    """Two intervals holding one read would compare it twice and could hide a
    disagreement behind an agreement."""

    with pytest.raises(gate.ParityError, match="not disjoint"):
        gate.check_strand_invariant(
            corpus.bam,
            corpus.fasta,
            str(tmp_path / "overlap"),
            intervals=[
                gate.Interval(CONTIG, 1, 20000),
                gate.Interval(CONTIG, 15000, LENGTH),
            ],
        )


# ------------------------------------------------------- the ordering mistake


def test_the_ordering_mistake_is_measured_and_its_reads_named(corpus, tmp_path):
    """Extracting a STRANDED region from a raw BAM, on purpose.

    ``retained_for_extraction`` filters on the raw flag; the splitter rewrites
    that flag.  So a ``chrT+`` extraction keeps ``flip_to_minus`` -- forward by
    the aligner, reverse by its splice sites -- and every later stage treats it as
    forward.  Nothing downstream revisits it, which is why the mistake has to be
    caught here.
    """

    report = gate.measure_ordering_violation(
        corpus.bam,
        corpus.fasta,
        "{}+:1-{}".format(CONTIG, BOUNDARY),
        str(tmp_path / "wrong"),
        gtf=corpus.gtf,
    )

    assert report["misassigned_read_names"] == ["flip_to_minus"]
    assert report["records_the_split_would_move"] == 1
    assert "ORDERING MISTAKE IS DETECTABLE" in report["verdict"]


def test_the_strandless_extraction_keeps_both_strands(corpus, tmp_path):
    """The reason the ordering can be fixed at all: a strandless region emits
    both orientations, so the split still has both to choose between."""

    report = run_invariant(corpus, tmp_path / "both", infer_read_orient=True)
    first = next(c for c in report["per_chunk"] if c["chunk_id"] == "chrT_00")

    # flip_to_plus (flag reverse) and flip_to_minus (flag forward) both survive
    # extraction, together with annot_to_plus and undecidable
    assert first["records_split_in_chunk"] == 4
    assert first["disagreed"] == 0


# ----------------------------------------------------------------- arm compare


EXPR_HEADER = [
    "gene_id",
    "transcript_id",
    "uniq_reads",
    "all_reads",
    "isoform_fraction",
    "unique_gene_read_fraction",
    "TPM",
    "exons",
    "introns",
    "splice_hash_code",
    "RPM_total_reads",
]
TRACKING_HEADER = [
    "gene_id",
    "transcript_id",
    "transcript_splice_hash_code",
    "num_exons",
    "mp_id",
    "read_name",
    "frac_assigned",
    "read_weight",
]


def expr_row(transcript_id, all_reads, tpm="1000.000"):
    return [
        "G^1",
        transcript_id,
        str(int(float(all_reads))),
        "{:.1f}".format(float(all_reads)),
        "1.000",
        "1.000",
        tpm,
        "chrT:(+)[[1, 100]]",
        "",
        "abc123",
        "100.000",
    ]


def tracking_row(read_name, transcript_id, frac="1.000000", mp_id="MP1"):
    return [
        "G^1",
        transcript_id,
        "abc123",
        "1",
        mp_id,
        read_name,
        frac,
        "1.000000",
    ]


def make_arm(
    root,
    expr_rows,
    tracking_rows,
    num_total_reads=100,
    severed_at_cuts=(),
    overhang=None,
    chunk_reads=None,
    strandless=False,
):
    """A minimal pipeline output directory: exactly what ``compare-arms`` reads."""

    root = str(root)
    os.makedirs(os.path.join(root, "merged"), exist_ok=True)
    os.makedirs(os.path.join(root, "cuts"), exist_ok=True)
    expr_path = os.path.join(root, "merged", "chunked.quant.expr")
    with open(expr_path, "wt") as fh:
        print("\t".join(EXPR_HEADER), file=fh)
        for row in expr_rows:
            print("\t".join(row), file=fh)
    track_path = os.path.join(root, "merged", "chunked.quant.tracking.gz")
    with gzip.open(track_path, "wt") as fh:
        print("\t".join(TRACKING_HEADER), file=fh)
        for row in tracking_rows:
            print("\t".join(row), file=fh)

    if severed_at_cuts:
        with open(os.path.join(root, "cuts", "plus.dropped_reads.txt"), "wt") as fh:
            for name in severed_at_cuts:
                print(name, file=fh)
    for chunk_id, names in (overhang or {}).items():
        cdir = os.path.join(root, "chunks", chunk_id)
        os.makedirs(cdir, exist_ok=True)
        with open(os.path.join(cdir, "chunk.dropped_reads.txt"), "wt") as fh:
            for name in names:
                print(name, file=fh)
    for chunk_id, names in (chunk_reads or {}).items():
        cdir = os.path.join(root, "chunks", chunk_id)
        os.makedirs(cdir, exist_ok=True)
        header = {"HD": {"VN": "1.6"}, "SQ": [{"SN": CONTIG, "LN": LENGTH}]}
        with pysam.AlignmentFile(
            os.path.join(cdir, "chunk.bam"), "wb", header=header
        ) as fh:
            for name in names:
                aln = pysam.AlignedSegment()
                aln.query_name = name
                aln.flag = 0
                aln.reference_id = 0
                aln.reference_start = 100
                aln.cigar = [(0, 50)]
                aln.query_sequence = "A" * 50
                aln.mapping_quality = 60
                fh.write(aln)

    with open(os.path.join(root, "outputs.json"), "wt") as fh:
        json.dump(
            {
                "num_total_reads": num_total_reads,
                "chunked": {
                    "quant_expr": expr_path,
                    "quant_tracking": track_path,
                    "expr_rows": len(expr_rows),
                    "tracking_rows": len(tracking_rows),
                },
            },
            fh,
        )
    with open(os.path.join(root, "timing.json"), "wt") as fh:
        json.dump({"strandless_chunks": strandless, "num_total_reads": num_total_reads}, fh)
    return root


BASE_EXPR = [expr_row("T1", 2), expr_row("T2", 3)]
BASE_TRACKING = [
    tracking_row("r1", "T1"),
    tracking_row("r2", "T1"),
    tracking_row("r3", "T2"),
]


def test_identical_arms_are_equivalent(tmp_path):
    a = make_arm(tmp_path / "a", BASE_EXPR, BASE_TRACKING)
    b = make_arm(tmp_path / "b", BASE_EXPR, BASE_TRACKING, strandless=True)

    report = gate.compare_arms(a, b)

    assert report["expr"]["byte_identical"] is True
    assert report["tracking"]["byte_identical"] is True
    assert report["verdict"].startswith("EQUIVALENT:")
    assert report["passed"] is True


def test_a_tracking_order_difference_is_reported_as_ordering_and_proven(tmp_path):
    """Ordering is claimed only when the row multisets are equal, which no
    content difference can satisfy. That is the proof, not the claim."""

    a = make_arm(tmp_path / "a", BASE_EXPR, BASE_TRACKING)
    b = make_arm(
        tmp_path / "b", BASE_EXPR, list(reversed(BASE_TRACKING)), strandless=True
    )

    report = gate.compare_arms(a, b)

    assert report["tracking"]["byte_identical"] is False
    assert report["tracking"]["rows_identical_in_order"] is False
    assert report["tracking"]["row_multisets_equal"] is True
    assert report["tracking"]["difference_is_ordering_only"] is True
    assert report["verdict"].startswith("EQUIVALENT UP TO ROW ORDER")
    assert report["passed"] is True


def test_a_transcript_emitted_twice_is_named_and_not_keyed_over(tmp_path):
    """The defect this gate found on its first real run.

    A strandless chunk's GTF carries BOTH orientations. Handed unchanged to both
    of the chunk's per-orientation quant runs, it made every transcript appear
    twice in the merged table -- each copy quantified against the other strand's
    reads, so neither is a copy of the other. Keying on transcript_id would pick
    one arbitrarily and report a plausible number; the gate must refuse to key and
    name the duplicates instead.
    """

    a = make_arm(tmp_path / "a", BASE_EXPR, BASE_TRACKING)
    b = make_arm(
        tmp_path / "b",
        BASE_EXPR + [expr_row("T1", 0), expr_row("T2", 0)],
        BASE_TRACKING,
        strandless=True,
    )

    report = gate.compare_arms(a, b)

    assert report["expr"]["duplicate_transcript_ids_b"] == {"T1": 2, "T2": 2}
    assert report["expr"]["duplicate_transcript_ids_a"] == {}
    assert report["expr"]["content_identical"] is False
    assert "not_comparable" in report["expr"]
    assert report["passed"] is False


def test_a_content_difference_is_never_called_ordering(tmp_path):
    changed = [
        tracking_row("r1", "T1", frac="0.500000"),
        tracking_row("r2", "T1"),
        tracking_row("r3", "T2"),
    ]
    a = make_arm(tmp_path / "a", BASE_EXPR, BASE_TRACKING)
    b = make_arm(tmp_path / "b", BASE_EXPR, changed, strandless=True)

    report = gate.compare_arms(a, b)

    assert report["tracking"]["row_multisets_equal"] is False
    assert report["tracking"]["difference_is_ordering_only"] is False
    assert report["tracking"]["frac_assigned_keys_differing"] == 1
    assert report["verdict"].startswith("DIVERGENT AND NOT EXPLAINED")
    assert report["passed"] is False


def test_mp_id_alone_does_not_make_the_arms_differ(tmp_path):
    """mp_id is a per-process MultiPath counter, not data: two runs of the same
    input number their multipaths differently."""

    renumbered = [
        tracking_row("r1", "T1", mp_id="MP9"),
        tracking_row("r2", "T1", mp_id="MP8"),
        tracking_row("r3", "T2", mp_id="MP7"),
    ]
    a = make_arm(tmp_path / "a", BASE_EXPR, BASE_TRACKING)
    b = make_arm(tmp_path / "b", BASE_EXPR, renumbered, strandless=True)

    report = gate.compare_arms(a, b)

    assert report["tracking"]["byte_identical"] is False
    assert report["tracking"]["rows_identical_in_order"] is True
    assert report["tracking"]["columns_excluded_as_incomparable"] == ["mp_id"]
    assert report["passed"] is True
    # and it is described as what it is, not as an ordering difference: the rows
    # are in the same order, so calling it ordering would be a false statement
    assert report["verdict"] == (
        "EQUIVALENT UP TO mp_id: every expr value agrees and every tracking row "
        "agrees in the same order, differing only in mp_id, which is a "
        "per-process counter rather than data"
    )


def test_a_differing_denominator_is_a_precondition_violation(tmp_path):
    """-N is RPM_total_reads' denominator and is passed identically to every
    chunk. While it differs no per-transcript comparison means anything."""

    a = make_arm(tmp_path / "a", BASE_EXPR, BASE_TRACKING, num_total_reads=100)
    b = make_arm(
        tmp_path / "b", BASE_EXPR, BASE_TRACKING, num_total_reads=99, strandless=True
    )

    report = gate.compare_arms(a, b)

    assert len(report["preconditions_violated"]) == 1
    assert "num_total_reads differs" in report["preconditions_violated"][0]


def test_divergence_confined_to_the_severed_difference_is_explained(tmp_path):
    """The arms cut in different places, so they sever different reads. A row
    that moved because of one of those reads is an input difference, not a
    defect -- and the read is named."""

    a = make_arm(
        tmp_path / "a",
        [expr_row("T1", 2), expr_row("T2", 3)],
        BASE_TRACKING,
        overhang={"chrT_plus_00": ["r9"]},
    )
    b = make_arm(
        tmp_path / "b",
        [expr_row("T1", 3), expr_row("T2", 3)],
        BASE_TRACKING + [tracking_row("r9", "T1")],
        strandless=True,
        chunk_reads={"chrT_00": ["r1", "r2", "r3", "r9"]},
    )

    report = gate.compare_arms(a, b)

    assert report["severed"]["symmetric_difference"] == ["r9"]
    assert report["attribution"]["fully_explained_by_severing"] is True
    assert report["attribution"]["transcripts_touched_by_severed_reads"] == ["T1"]
    assert report["verdict"].startswith("DIVERGENT, FULLY EXPLAINED BY SEVERING")
    assert report["passed"] is True

    named = report["reads_absent_by_name"]["absent_from_strand_first"]
    assert [record["read_name"] for record in named] == ["r9"]
    assert named[0]["reason"] == "overhang_dropped_at_chrT_plus_00"


def test_divergence_with_no_severed_read_behind_it_is_unexplained(tmp_path):
    """The case that must never be absorbed: identical inputs, different answer.
    The transcript is named, because a count would not let anyone find it."""

    a = make_arm(tmp_path / "a", [expr_row("T1", 2), expr_row("T2", 3)], BASE_TRACKING)
    b = make_arm(
        tmp_path / "b",
        [expr_row("T1", 7), expr_row("T2", 3)],
        BASE_TRACKING,
        strandless=True,
    )

    report = gate.compare_arms(a, b)

    assert report["severed"]["symmetric_difference"] == []
    assert report["attribution"]["fully_explained_by_severing"] is False
    assert report["attribution"]["unexplained_expr_columns"]["all_reads"] == ["T1"]
    assert report["expr"]["numeric"]["all_reads"]["max_abs_delta"] == 5.0
    assert report["verdict"].startswith("DIVERGENT AND NOT EXPLAINED")
    assert report["passed"] is False


def test_tpm_is_reported_separately_because_it_is_globally_coupled(tmp_path):
    """TPM's denominator is the merged table's total all_reads, so one differing
    read rescales every row. Counting it as unexplained would condemn the whole
    table on the strength of a single legitimate severing."""

    a = make_arm(
        tmp_path / "a",
        [expr_row("T1", 2, tpm="400000.000"), expr_row("T2", 3, tpm="600000.000")],
        BASE_TRACKING,
        overhang={"chrT_plus_00": ["r9"]},
    )
    b = make_arm(
        tmp_path / "b",
        [expr_row("T1", 3, tpm="500000.000"), expr_row("T2", 3, tpm="500000.000")],
        BASE_TRACKING + [tracking_row("r9", "T1")],
        strandless=True,
        chunk_reads={"chrT_00": ["r9"]},
    )

    report = gate.compare_arms(a, b)

    assert report["attribution"]["tpm_transcripts_differing"] == 2
    assert report["attribution"]["tpm_max_abs_delta"] == 100000.0
    # T2's TPM moved without T2 having a severed read, and that is not a defect
    assert "TPM" not in report["attribution"]["unexplained_expr_columns"]
    assert report["attribution"]["fully_explained_by_severing"] is True


def test_a_read_in_neither_severed_set_is_still_named_and_explained(tmp_path):
    """A read absent from one arm's output that no cut lost: the gate resolves
    the reason from that arm's own artifacts rather than guessing."""

    a = make_arm(tmp_path / "a", BASE_EXPR, BASE_TRACKING + [tracking_row("r9", "T1")])
    b = make_arm(
        tmp_path / "b",
        BASE_EXPR,
        BASE_TRACKING,
        strandless=True,
        chunk_reads={"chrT_00": ["r1", "r2", "r3", "r9"]},
    )

    report = gate.compare_arms(a, b)

    named = report["reads_absent_by_name"]["absent_from_strandless"]
    assert [record["read_name"] for record in named] == ["r9"]
    assert named[0]["reason"] == (
        "in_chunk_input_chrT_00_but_assigned_no_transcript"
    )
    assert report["passed"] is False


def test_the_chunk_inputs_are_compared_by_name_not_by_count(tmp_path):
    """The arms' inputs, diffed before their outputs are trusted.

    Chunk BAMs carry rebased coordinates and the two arms cut in different places,
    so a position is not comparable between them and a name is. A count would say
    "one fewer read" and leave the reader with no way to find which.
    """

    a = make_arm(
        tmp_path / "a",
        BASE_EXPR,
        BASE_TRACKING,
        chunk_reads={"chrT_plus_00": ["r1", "r2", "r3", "r9"]},
    )
    b = make_arm(
        tmp_path / "b",
        BASE_EXPR,
        BASE_TRACKING,
        strandless=True,
        overhang={"chrT_00": ["r9"]},
        chunk_reads={"chrT_00": ["r1", "r2", "r3"]},
    )

    report = gate.compare_arms(a, b)
    inputs = report["chunk_inputs"]

    assert inputs["names_identical"] is False
    assert inputs["distinct_read_names_strand_first"] == 4
    assert inputs["distinct_read_names_strandless"] == 3
    assert [r["read_name"] for r in inputs["only_in_strand_first"]] == ["r9"]
    assert inputs["only_in_strand_first"][0]["reason"] == "overhang_dropped_at_chrT_00"
    assert inputs["only_in_strand_first"][0]["strand_first_chunks"] == ["chrT_plus_00"]
    assert inputs["only_in_strandless"] == []
