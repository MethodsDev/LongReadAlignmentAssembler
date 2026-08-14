#!/usr/bin/env python3

"""Tests for util/misc/extract_contig_region_inputs.py.

The utility partitions a contig-strand into strictly disjoint chunks: no halo, no
widening. A boundary must clear every same-strand gene locus, with a margin, and
that rule is absolute -- a locus cannot be dropped and re-derived. An alignment
crossing a boundary does NOT forbid it: the read is dropped, counted and named.
These tests hold it to that contract -- every retained record in exactly one chunk
or in the dropped set, never both and never neither -- and to the mechanical
contract that the extracted inputs are self-consistent.
"""

import importlib.util
import json
import re
from importlib.machinery import SourceFileLoader
from pathlib import Path

import pysam
import pytest


def _load_extractor():
    path = (
        Path(__file__).resolve().parents[1]
        / "util"
        / "misc"
        / "extract_contig_region_inputs.py"
    )
    loader = SourceFileLoader("extract_contig_region_inputs_under_test", str(path))
    spec = importlib.util.spec_from_loader(loader.name, loader)
    module = importlib.util.module_from_spec(spec)
    loader.exec_module(module)
    return module


extractor = _load_extractor()


class Fixture:
    """Builds a synthetic contig with an indexed FASTA, BAM and GTF."""

    def __init__(self, tmp_path, contig="chrT", length=30000):
        self.tmp_path = tmp_path
        self.contig = contig
        self.length = length
        self.reads = []  # (name, strand, blocks, flag)
        self.transcripts = []  # (gene_id, transcript_id, strand, exons)

    def add_read(self, name, strand, blocks, flag=None):
        self.reads.append((name, strand, [tuple(b) for b in blocks], flag))
        return self

    def add_transcript(self, gene_id, transcript_id, strand, exons):
        self.transcripts.append(
            (gene_id, transcript_id, strand, [tuple(e) for e in exons])
        )
        return self

    # -- expected values, in source coordinates ---------------------------------

    def read_blocks(self, name):
        for read_name, _, blocks, _ in self.reads:
            if read_name == name:
                return blocks
        raise KeyError(name)

    def transcript_exons(self, transcript_id):
        for _, tid, _, exons in self.transcripts:
            if tid == transcript_id:
                return exons
        raise KeyError(transcript_id)

    # -- materialisation --------------------------------------------------------

    def build(self):
        self.fasta = str(self.tmp_path / "genome.fa")
        # deterministic, unremarkable sequence: these tests are about coordinates
        sequence = ("ACGT" * (self.length // 4 + 1))[: self.length]
        with open(self.fasta, "wt") as ofh:
            print(">{}".format(self.contig), file=ofh)
            for i in range(0, self.length, 60):
                print(sequence[i : i + 60], file=ofh)
        pysam.faidx(self.fasta)

        self.bam = str(self.tmp_path / "reads.bam")
        header = {
            "HD": {"VN": "1.6", "SO": "coordinate"},
            "SQ": [{"SN": self.contig, "LN": self.length}],
        }
        records = [self._alignment(*read) for read in self.reads]
        records.sort(key=lambda aln: aln.reference_start)
        with pysam.AlignmentFile(self.bam, "wb", header=header) as ofh:
            for record in records:
                ofh.write(record)
        pysam.index(self.bam)

        self.gtf = str(self.tmp_path / "annot.gtf")
        with open(self.gtf, "wt") as ofh:
            for gene_id, transcript_id, strand, exons in self.transcripts:
                lend = min(e[0] for e in exons)
                rend = max(e[1] for e in exons)
                attrs = 'gene_id "{}";'.format(gene_id)
                print(
                    "\t".join(
                        (
                            self.contig,
                            "test",
                            "gene",
                            str(lend),
                            str(rend),
                            ".",
                            strand,
                            ".",
                            attrs,
                        )
                    ),
                    file=ofh,
                )
                tx_attrs = 'gene_id "{}"; transcript_id "{}";'.format(
                    gene_id, transcript_id
                )
                print(
                    "\t".join(
                        (
                            self.contig,
                            "test",
                            "transcript",
                            str(lend),
                            str(rend),
                            ".",
                            strand,
                            ".",
                            tx_attrs,
                        )
                    ),
                    file=ofh,
                )
                for number, (exon_lend, exon_rend) in enumerate(sorted(exons), 1):
                    print(
                        "\t".join(
                            (
                                self.contig,
                                "test",
                                "exon",
                                str(exon_lend),
                                str(exon_rend),
                                ".",
                                strand,
                                ".",
                                tx_attrs + ' exon_number {};'.format(number),
                            )
                        ),
                        file=ofh,
                    )
        return self

    def _alignment(self, name, strand, blocks, flag):
        aln = pysam.AlignedSegment()
        aln.query_name = name
        aln.flag = (16 if strand == "-" else 0) if flag is None else flag
        aln.reference_id = 0
        aln.reference_start = blocks[0][0] - 1
        aln.mapping_quality = 60
        cigar = []
        for i, (lend, rend) in enumerate(blocks):
            if i:
                cigar.append((3, lend - blocks[i - 1][1] - 1))
            cigar.append((0, rend - lend + 1))
        aln.cigar = cigar
        query_length = sum(rend - lend + 1 for lend, rend in blocks)
        aln.query_sequence = "A" * query_length
        aln.query_qualities = pysam.qualitystring_to_array("I" * query_length)
        aln.set_tag("NM", 0)
        return aln


def _emitted_alignments(bam_filename):
    """name -> (start, end, cigarstring) in mini-contig coordinates."""
    out = {}
    with pysam.AlignmentFile(bam_filename, "rb") as bam:
        for aln in bam:
            out[aln.query_name] = (
                aln.reference_start + 1,
                aln.reference_end,
                aln.cigarstring,
            )
    return out


def _emitted_gtf(gtf_filename):
    """transcript_id -> {'exons': [...], 'bounds': (lend, rend), 'features': [...]}"""
    models = {}
    with open(gtf_filename, "rt") as fh:
        for line in fh:
            vals = line.rstrip("\n").split("\t")
            transcript_id = extractor._attribute(vals[8], "transcript_id")
            if transcript_id is None:
                continue
            model = models.setdefault(
                transcript_id, {"exons": [], "bounds": None, "features": []}
            )
            lend, rend = int(vals[3]), int(vals[4])
            model["features"].append(vals[2])
            if vals[2] == "exon":
                model["exons"].append((lend, rend))
            elif vals[2] == "transcript":
                model["bounds"] = (lend, rend)
    for model in models.values():
        model["exons"].sort()
    return models


def _single_gene_fixture(tmp_path, strand):
    """Loci and reads flanking, touching and inside a partition [10001, 20000]."""

    fixture = Fixture(tmp_path, length=30000)
    # left neighbour, comfortably clear of the boundary
    fixture.add_transcript("gLeft", "tLeft", strand, [(5000, 5200), (5500, 5700)])
    # a locus whose first exon starts exactly on the left boundary
    fixture.add_transcript(
        "gEdgeL", "tEdgeL", strand, [(10001, 10200), (10500, 10800)]
    )
    fixture.add_transcript(
        "gInner", "tInner", strand, [(14000, 14300), (15000, 15400)]
    )
    # a locus whose last exon ends exactly on the right boundary
    fixture.add_transcript(
        "gEdgeR", "tEdgeR", strand, [(19000, 19300), (19700, 20000)]
    )
    fixture.add_transcript("gRight", "tRight", strand, [(25000, 25400)])

    # reads touching each edge exactly, and one spliced read inside
    fixture.add_read("rEdgeL", strand, [(10001, 10200), (10500, 10800)])
    fixture.add_read("rInner", strand, [(14000, 14300), (15000, 15400)])
    fixture.add_read("rEdgeR", strand, [(19000, 19300), (19700, 20000)])
    fixture.add_read("rOutsideLeft", strand, [(5000, 5200)])
    fixture.add_read("rOutsideRight", strand, [(25000, 25400)])
    # opposite strand, inside the region: must be filtered out entirely
    other = "-" if strand == "+" else "+"
    fixture.add_read("rOtherStrand", other, [(14100, 14400)])
    return fixture.build()


@pytest.mark.parametrize("strand", ["+", "-"])
def test_coordinate_round_trip_recovers_source_coordinates(tmp_path, strand):
    fixture = _single_gene_fixture(tmp_path, strand)
    region = "{}{}:10001-20000".format(fixture.contig, strand)

    manifest = extractor.extract_partition(
        genome_fa=fixture.fasta,
        bam=fixture.bam,
        region=region,
        output_prefix=str(tmp_path / "part"),
        gtf=fixture.gtf,
        margin=0,
    )

    offset = manifest["offset"]
    assert offset == 10000
    assert manifest["mini_length"] == 10000
    assert (manifest["partition_lend"], manifest["partition_rend"]) == (10001, 20000)

    emitted = _emitted_alignments(manifest["files"]["bam"])
    # exactly the same-strand reads inside the partition, both edge-touching ones
    assert set(emitted) == {"rEdgeL", "rInner", "rEdgeR"}
    for name, (start, end, cigar) in emitted.items():
        blocks = fixture.read_blocks(name)
        assert start + offset == blocks[0][0], name
        assert end + offset == blocks[-1][1], name
        assert 1 <= start and end <= manifest["mini_length"], name
        # block structure survives translation intact
        translated = [(lend + offset, rend + offset) for lend, rend in _cigar_blocks(start, cigar)]
        assert translated == blocks, name

    models = _emitted_gtf(manifest["files"]["gtf"])
    assert set(models) == {"tEdgeL", "tInner", "tEdgeR"}
    for transcript_id, model in models.items():
        source_exons = fixture.transcript_exons(transcript_id)
        assert [
            (lend + offset, rend + offset) for lend, rend in model["exons"]
        ] == sorted(source_exons), transcript_id
        assert model["bounds"] is not None
        assert (
            model["bounds"][0] + offset,
            model["bounds"][1] + offset,
        ) == (
            min(e[0] for e in source_exons),
            max(e[1] for e in source_exons),
        ), transcript_id
        assert all(
            1 <= lend and rend <= manifest["mini_length"]
            for lend, rend in model["exons"]
        )


def _cigar_blocks(start, cigarstring):
    """1-based blocks from a CIGAR string, splitting only at N."""
    blocks = []
    position = start
    block_start = start
    for length, op in re.findall(r"(\d+)([A-Z=])", cigarstring):
        length = int(length)
        if op in ("M", "D", "=", "X"):
            position += length
        elif op == "N":
            blocks.append((block_start, position - 1))
            position += length
            block_start = position
    blocks.append((block_start, position - 1))
    return blocks


def test_output_bam_header_declares_the_mini_contig_length(tmp_path):
    fixture = _single_gene_fixture(tmp_path, "+")
    manifest = extractor.extract_partition(
        genome_fa=fixture.fasta,
        bam=fixture.bam,
        region="{}+:10001-20000".format(fixture.contig),
        output_prefix=str(tmp_path / "part"),
        gtf=fixture.gtf,
        margin=0,
    )
    with pysam.AlignmentFile(manifest["files"]["bam"], "rb") as bam:
        sq = bam.header.to_dict()["SQ"]
    assert sq == [{"SN": fixture.contig, "LN": 10000}]
    # the source contig is longer, so inheriting its header would be inconsistent
    assert fixture.length == 30000
    with pysam.FastaFile(manifest["files"]["fasta"]) as fasta:
        assert fasta.get_reference_length(fixture.contig) == 10000


def test_manifest_records_partition_bounds_margin_and_offset(tmp_path):
    fixture = _single_gene_fixture(tmp_path, "+")
    manifest = extractor.extract_partition(
        genome_fa=fixture.fasta,
        bam=fixture.bam,
        region="{}+:10001-20000".format(fixture.contig),
        output_prefix=str(tmp_path / "part"),
        gtf=fixture.gtf,
        margin=0,
    )
    assert manifest["disjoint_hard_cut"] is True
    assert manifest["margin"] == 0
    assert manifest["contig_length"] == 30000
    assert manifest["offset"] == manifest["partition_lend"] - 1
    assert manifest["counts"]["alignments_emitted"] == 3
    assert manifest["counts"]["gtf_transcripts_emitted"] == 3
    assert manifest["emitted_transcript_ids"] == ["tEdgeL", "tEdgeR", "tInner"]


def _overhang_fixture(tmp_path, side):
    fixture = Fixture(tmp_path, length=30000)
    fixture.add_transcript("gIn", "tIn", "+", [(12000, 12300), (13000, 13400)])
    if side == "left":
        # first block and intron left of the boundary, last block inside: the
        # exact shape of the read lost at 10 genes per segment
        fixture.add_read("rSpanLeft", "+", [(9500, 9700), (10200, 11000)])
    else:
        fixture.add_read("rSpanRight", "+", [(19500, 19800), (20200, 21000)])
    return fixture.build()


@pytest.mark.parametrize("side", ["left", "right"])
def test_read_spanning_a_boundary_is_dropped_by_name_not_truncated(tmp_path, side):
    """Policy: the read is dropped, the boundary stands.

    Previously this position was refused outright. Refusing made granularity a
    function of where alignments happen to leave gaps; dropping makes it a
    function of the requested span, and the price is counted instead of hidden.
    What must NOT happen either way is a truncated alignment or an out-of-range
    record, so the emitted BAM is checked for both.
    """

    fixture = _overhang_fixture(tmp_path, side)
    prefix = str(tmp_path / "part")
    read_name = "rSpanLeft" if side == "left" else "rSpanRight"

    manifest = extractor.extract_partition(
        genome_fa=fixture.fasta,
        bam=fixture.bam,
        region="{}+:10001-20000".format(fixture.contig),
        output_prefix=prefix,
        gtf=fixture.gtf,
        margin=0,
    )

    # counted, named in the manifest, and named in a file a filter can consume
    assert manifest["counts"]["alignments_dropped_overhang"] == 1
    assert manifest["dropped_read_names"] == [read_name]
    assert Path(manifest["files"]["dropped_reads"]).read_text().split() == [read_name]

    # not truncated, and not present in any form
    emitted = _emitted_alignments(manifest["files"]["bam"])
    assert read_name not in emitted
    for start, end, _ in emitted.values():
        assert start >= 1
        assert end <= manifest["mini_length"]

    # the chunk was written, unlike under the old refusal
    for suffix in (".fa", ".bam", ".gtf", ".partition.json", ".dropped_reads.txt"):
        assert Path(prefix + suffix).exists()


@pytest.mark.parametrize("side", ["left", "right"])
def test_a_spanning_alignment_no_longer_makes_a_position_inadmissible(tmp_path, side):
    """The alignment half of the old admissibility rule is gone.

    ``admissibility_offenders`` reports annotated loci only, so a boundary crossed
    by nothing but reads is legal. The reads that cross it are dropped by both
    neighbours -- each contains neither -- and everything else lands exactly once.
    """

    fixture = _overhang_fixture(tmp_path, side)
    annotation = extractor.load_gtf(fixture.gtf, fixture.contig, "+")
    boundary = 10000 if side == "left" else 20000
    read_name = "rSpanLeft" if side == "left" else "rSpanRight"

    region = extractor.parse_region(
        "{}+:1-{}".format(fixture.contig, boundary)
    )
    assert (
        extractor.admissibility_offenders(annotation, region, fixture.length, margin=0)
        == []
    )

    # islands are annotated loci only, so the spanning read forms none
    islands = extractor.find_islands(annotation, fixture.contig, "+", margin=0)
    zones = extractor.cut_zones(islands, 1, fixture.length)
    assert any(lo <= boundary <= hi for lo, hi in zones)

    partitions = [(1, boundary), (boundary + 1, fixture.length)]
    seen = {}
    dropped = []
    for index, (lend, rend) in enumerate(partitions):
        manifest = extractor.extract_partition(
            genome_fa=fixture.fasta,
            bam=fixture.bam,
            region="{}+:{}-{}".format(fixture.contig, lend, rend),
            output_prefix=str(tmp_path / "zone{}".format(index)),
            gtf=fixture.gtf,
            margin=0,
        )
        offset = manifest["offset"]
        dropped.extend(manifest["dropped_read_names"])
        for name, (start, end, cigar) in _emitted_alignments(
            manifest["files"]["bam"]
        ).items():
            assert name not in seen
            seen[name] = index
            # complete: every block recovers its source coordinates
            assert [
                (lo + offset, hi + offset) for lo, hi in _cigar_blocks(start, cigar)
            ] == fixture.read_blocks(name), name

    # the severed read is reported by both neighbours and emitted by neither
    assert dropped == [read_name, read_name]
    assert read_name not in seen
    source = {name for name, _, _, _ in fixture.reads}
    assert set(seen) | {read_name} == source


def _straddle_fixture(tmp_path):
    """Loci straddling both intended boundaries, plus one wholly inside."""

    fixture = Fixture(tmp_path, length=40000)
    fixture.add_transcript(
        "gStraddleL", "tStraddleL", "+", [(9800, 9900), (10100, 10400)]
    )
    fixture.add_transcript("gInner", "tInner", "+", [(14000, 14300), (15000, 15400)])
    fixture.add_transcript(
        "gStraddleR", "tStraddleR", "+", [(19900, 19950), (20100, 20500)]
    )
    fixture.add_read("rInner", "+", [(14000, 14300), (15000, 15400)])
    fixture.add_read("rStraddleL", "+", [(9800, 9900), (10100, 10400)])
    fixture.add_read("rStraddleR", "+", [(19900, 19950), (20100, 20500)])
    return fixture.build()


def _edge_shapes_fixture(tmp_path):
    """Loci adjacent to, and nested inside, the boundaries at 10000 and 20000."""

    fixture = Fixture(tmp_path, length=40000)
    # ends exactly ON the left cut, so it belongs to the left neighbour
    fixture.add_transcript("gAdjOutL", "tAdjOutL", "+", [(9000, 9600), (9800, 10000)])
    # starts exactly on the first base of the partition
    fixture.add_transcript("gAdjInL", "tAdjInL", "+", [(10001, 10200), (10400, 10600)])
    fixture.add_transcript("gNested", "tNested", "+", [(12000, 12400), (12900, 13000)])
    # ends exactly on the last base of the partition
    fixture.add_transcript("gAdjInR", "tAdjInR", "+", [(19000, 19200), (19800, 20000)])
    # starts exactly on the first base after the partition
    fixture.add_transcript("gAdjOutR", "tAdjOutR", "+", [(20001, 20400)])
    fixture.add_read("rAdjInL", "+", [(10001, 10200), (10400, 10600)])
    fixture.add_read("rAdjInR", "+", [(19000, 19200), (19800, 20000)])
    fixture.add_read("rAdjOutL", "+", [(9000, 9600)])
    fixture.add_read("rAdjOutR", "+", [(20001, 20400)])
    return fixture.build()


def test_locus_straddling_either_boundary_is_rejected_by_name(tmp_path):
    fixture = _straddle_fixture(tmp_path)
    annotation = extractor.load_gtf(fixture.gtf, fixture.contig, "+")
    prefix = str(tmp_path / "bad")

    with pytest.raises(extractor.ExtractionError) as excinfo:
        extractor.extract_partition(
            genome_fa=fixture.fasta,
            bam=fixture.bam,
            region="{}+:10001-20000".format(fixture.contig),
            output_prefix=prefix,
            gtf=fixture.gtf,
            margin=0,
        )
    message = str(excinfo.value)
    assert "gStraddleL" in message and "gStraddleR" in message
    assert "straddles" in message
    for suffix in (".fa", ".bam", ".gtf", ".partition.json"):
        assert not Path(prefix + suffix).exists()

    # each edge is independently reported, so a one-sided fault is nameable too
    left_only = extractor.admissibility_offenders(
        annotation,
        extractor.parse_region("{}+:10001-40000".format(fixture.contig)),
        fixture.length,
        margin=0,
    )
    assert [o for o in left_only if "gStraddleL" in o]
    assert not [o for o in left_only if "gStraddleR" in o]
    right_only = extractor.admissibility_offenders(
        annotation,
        extractor.parse_region("{}+:1-20000".format(fixture.contig)),
        fixture.length,
        margin=0,
    )
    assert [o for o in right_only if "gStraddleR" in o]
    assert not [o for o in right_only if "gStraddleL" in o]


def test_boundary_touching_loci_are_emitted_whole_at_both_edges(tmp_path):
    fixture = _edge_shapes_fixture(tmp_path)
    manifest = extractor.extract_partition(
        genome_fa=fixture.fasta,
        bam=fixture.bam,
        region="{}+:10001-20000".format(fixture.contig),
        output_prefix=str(tmp_path / "edges"),
        gtf=fixture.gtf,
        margin=0,
    )
    offset = manifest["offset"]
    models = _emitted_gtf(manifest["files"]["gtf"])
    # adjacent-but-outside belongs to the neighbour; touching-inside and nested
    # are emitted, and emitted whole
    assert set(models) == {"tAdjInL", "tNested", "tAdjInR"}
    for transcript_id, model in models.items():
        assert [
            (lo + offset, hi + offset) for lo, hi in model["exons"]
        ] == sorted(fixture.transcript_exons(transcript_id)), transcript_id
    # the exon touching each edge lands on the first and last base of the contig
    assert models["tAdjInL"]["exons"][0][0] == 1
    assert models["tAdjInR"]["exons"][-1][1] == manifest["mini_length"]
    assert set(_emitted_alignments(manifest["files"]["bam"])) == {
        "rAdjInL",
        "rAdjInR",
    }


# Boundaries inside the annotated cut zones of ``_straddle_fixture`` at margin 0:
# loci occupy 9800-10400, 14000-15400 and 19900-20500, so 12000 and 17000 are both
# clear of every locus on either side.
_STRADDLE_PARTITIONS = [(1, 12000), (12001, 17000), (17001, 40000)]


def test_no_transcript_is_ever_emitted_partially(tmp_path):
    """Even straddling loci land whole, in exactly one partition."""

    fixture = _straddle_fixture(tmp_path)
    annotation = extractor.load_gtf(fixture.gtf, fixture.contig, "+")

    # the chosen boundaries really are clear of every locus
    for lend, rend in _STRADDLE_PARTITIONS:
        region = extractor.parse_region(
            "{}+:{}-{}".format(fixture.contig, lend, rend)
        )
        assert (
            extractor.admissibility_offenders(
                annotation, region, fixture.length, margin=0
            )
            == []
        )

    seen = {}
    for index, (lend, rend) in enumerate(_STRADDLE_PARTITIONS):
        manifest = extractor.extract_partition(
            genome_fa=fixture.fasta,
            bam=fixture.bam,
            region="{}+:{}-{}".format(fixture.contig, lend, rend),
            output_prefix=str(tmp_path / "p{}".format(index)),
            gtf=fixture.gtf,
            margin=0,
        )
        offset = manifest["offset"]
        for transcript_id, model in _emitted_gtf(manifest["files"]["gtf"]).items():
            assert transcript_id not in seen, transcript_id
            seen[transcript_id] = True
            assert [
                (lo + offset, hi + offset) for lo, hi in model["exons"]
            ] == sorted(fixture.transcript_exons(transcript_id)), transcript_id
            assert model["features"].count("transcript") == 1

    assert set(seen) == {tid for _, tid, _, _ in fixture.transcripts}


def _multi_locus_fixture(tmp_path):
    """Six well-separated loci with reads, plus two overlapping loci."""

    fixture = Fixture(tmp_path, length=60000)
    layout = [
        ("g1", "t1", "+", [(2000, 2300), (2800, 3100)]),
        ("g2", "t2", "+", [(9000, 9400), (9900, 10200)]),
        # g3 and g3b overlap: one island, never separable
        ("g3", "t3", "+", [(17000, 17300), (17800, 18200)]),
        ("g3b", "t3b", "+", [(18100, 18400), (18900, 19100)]),
        ("g4", "t4", "+", [(27000, 27500)]),
        ("g5", "t5", "+", [(36000, 36300), (36900, 37200)]),
        ("g6", "t6", "+", [(45000, 45400), (46000, 46300)]),
        ("m1", "mt1", "-", [(6000, 6300), (6900, 7200)]),
        ("m2", "mt2", "-", [(23000, 23400)]),
        ("m3", "mt3", "-", [(41000, 41300), (41900, 42200)]),
    ]
    for gene_id, transcript_id, strand, exons in layout:
        fixture.add_transcript(gene_id, transcript_id, strand, exons)
        # two reads per locus, one matching the exon chain, one single-block
        fixture.add_read(
            "r_{}_spliced".format(transcript_id),
            strand,
            exons if len(exons) > 1 else [exons[0]],
        )
        fixture.add_read(
            "r_{}_single".format(transcript_id),
            strand,
            [(exons[0][0] + 20, exons[0][1] - 20)],
        )
    # intergenic reads, which form islands of their own
    fixture.add_read("r_intergenic_a", "+", [(31000, 31500)])
    fixture.add_read("r_intergenic_b", "-", [(50000, 50600)])
    return fixture.build()


# Boundaries inside the annotated cut zones of ``_multi_locus_fixture`` at margin
# 200, on BOTH strands. Plus-strand loci leave the zones [3300, 8799],
# [10400, 16799], [19300, 26799], [27700, 35799], [37400, 44799], [46500, 59999];
# minus-strand loci leave [1, 5799], [7400, 22799], [23600, 40799], [42400, 59999].
# 12000, 24000 and 40000 sit inside a zone on either strand.
_MULTI_LOCUS_PARTITIONS = [
    (1, 12000),
    (12001, 24000),
    (24001, 40000),
    (40001, 60000),
]


@pytest.mark.parametrize("strand", ["+", "-"])
def test_every_record_lands_in_exactly_one_partition(tmp_path, strand):
    fixture = _multi_locus_fixture(tmp_path)
    annotation = extractor.load_gtf(fixture.gtf, fixture.contig, strand)

    # the partitions tile the contig-strand without gap or overlap
    assert _MULTI_LOCUS_PARTITIONS[0][0] == 1
    assert _MULTI_LOCUS_PARTITIONS[-1][1] == fixture.length
    for (_, rend), (lend, _) in zip(
        _MULTI_LOCUS_PARTITIONS, _MULTI_LOCUS_PARTITIONS[1:]
    ):
        assert lend == rend + 1

    source_transcripts = {
        tid for _, tid, tx_strand, _ in fixture.transcripts if tx_strand == strand
    }
    source_reads = {
        name for name, read_strand, _, _ in fixture.reads if read_strand == strand
    }

    transcript_homes = {}
    read_homes = {}
    dropped = set()
    for index, (lend, rend) in enumerate(_MULTI_LOCUS_PARTITIONS):
        region = extractor.parse_region(
            "{}{}:{}-{}".format(fixture.contig, strand, lend, rend)
        )
        assert (
            extractor.admissibility_offenders(
                annotation, region, fixture.length, margin=200
            )
            == []
        ), region
        manifest = extractor.extract_partition(
            genome_fa=fixture.fasta,
            bam=fixture.bam,
            region=region,
            output_prefix=str(tmp_path / "p{}_{}".format(strand, index)),
            gtf=fixture.gtf,
            margin=200,
        )
        offset = manifest["offset"]
        dropped.update(manifest["dropped_read_names"])
        for transcript_id, model in _emitted_gtf(manifest["files"]["gtf"]).items():
            assert transcript_id not in transcript_homes, transcript_id
            transcript_homes[transcript_id] = index
            assert [
                (lo + offset, hi + offset) for lo, hi in model["exons"]
            ] == sorted(fixture.transcript_exons(transcript_id))
        for name, (start, end, cigar) in _emitted_alignments(
            manifest["files"]["bam"]
        ).items():
            assert name not in read_homes, name
            read_homes[name] = index
            blocks = fixture.read_blocks(name)
            assert [
                (lo + offset, hi + offset) for lo, hi in _cigar_blocks(start, cigar)
            ] == blocks, name

    assert set(transcript_homes) == source_transcripts
    # these boundaries sever nothing, so the accounting closes with an empty
    # dropped set rather than by absorbing the difference
    assert dropped == set()
    assert set(read_homes) == source_reads


def test_a_straddling_read_is_dropped_and_every_other_record_lands_once(tmp_path):
    """The accounting identity: emitted + dropped == retained, disjointly."""

    fixture = _multi_locus_fixture(tmp_path)
    # a read deliberately laid across the 24000 boundary
    fixture.add_read("rAcross24000", "+", [(23800, 24200)])
    fixture.build()

    read_homes = {}
    dropped = []
    for index, (lend, rend) in enumerate(_MULTI_LOCUS_PARTITIONS):
        manifest = extractor.extract_partition(
            genome_fa=fixture.fasta,
            bam=fixture.bam,
            region="{}+:{}-{}".format(fixture.contig, lend, rend),
            output_prefix=str(tmp_path / "acc{}".format(index)),
            gtf=fixture.gtf,
            margin=200,
        )
        dropped.extend(manifest["dropped_read_names"])
        for name in _emitted_alignments(manifest["files"]["bam"]):
            assert name not in read_homes, name
            read_homes[name] = index

    assert dropped == ["rAcross24000", "rAcross24000"]
    source = {name for name, strand, _, _ in fixture.reads if strand == "+"}
    assert set(read_homes) | set(dropped) == source
    assert not (set(read_homes) & set(dropped))


@pytest.mark.parametrize("strand", ["+", "-"])
def test_no_read_overlaps_a_transcript_outside_its_own_partition(tmp_path, strand):
    """Disjointness gives candidate-set closure, checked against OVERLAP.

    Overlap, not splice compatibility: LRAA's majority-voting fallback assigns on
    a shared graph node rather than a matching intron chain, so overlap is the
    weaker and safer notion to verify. Reads dropped for overhang are out of
    scope: they are emitted nowhere, so they have no candidate set to be short.
    """

    fixture = _multi_locus_fixture(tmp_path)
    annotation = extractor.load_gtf(fixture.gtf, fixture.contig, strand)

    assert (
        extractor.overlap_offenders(
            fixture.bam,
            annotation,
            fixture.contig,
            strand,
            _MULTI_LOCUS_PARTITIONS,
        )
        == []
    )

    # and the check is not vacuous: a boundary drawn through a locus offends
    locus = next(
        gene for gene in annotation.genes.values() if len(gene.transcript_ids) > 0
    )
    bad_cut = (locus.lend + locus.rend) // 2
    bad = [(1, bad_cut), (bad_cut + 1, fixture.length)]
    offenders = extractor.overlap_offenders(
        fixture.bam, annotation, fixture.contig, strand, bad
    )
    assert offenders, "cutting through a locus must be reported"


def test_islands_bound_granularity_and_zones_are_exactly_the_legal_cuts(tmp_path):
    fixture = _multi_locus_fixture(tmp_path)
    annotation = extractor.load_gtf(fixture.gtf, fixture.contig, "+")
    islands = extractor.find_islands(annotation, fixture.contig, "+", margin=200)

    # every locus is inside exactly one island. Alignments no longer form islands:
    # under the drop-overhang policy an alignment is a cost, not a prohibition.
    assert sum(len(island.gene_ids) for island in islands) == len(annotation.genes)
    assert all(not hasattr(island, "n_alignments") for island in islands)

    # the overlapping pair g3/g3b is welded into one island
    welded = [island for island in islands if len(island.gene_ids) > 1]
    assert [sorted(island.gene_ids) for island in welded] == [["g3", "g3b"]]

    zones = extractor.cut_zones(islands, 1, fixture.length)
    for lo, hi in zones:
        for cut in (lo, (lo + hi) // 2, hi):
            region = extractor.parse_region(
                "{}+:1-{}".format(fixture.contig, cut)
            )
            assert (
                extractor.admissibility_offenders(
                    annotation, region, fixture.length, margin=200
                )
                == []
            ), cut
    for island in islands:
        inside = (island.lend + island.rend) // 2
        region = extractor.parse_region("{}+:1-{}".format(fixture.contig, inside))
        assert extractor.admissibility_offenders(
            annotation, region, fixture.length, margin=200
        ), inside


def test_margin_demands_clearance_beyond_mere_non_crossing(tmp_path):
    fixture = Fixture(tmp_path, length=30000)
    fixture.add_transcript("gA", "tA", "+", [(9000, 9900)])
    fixture.add_transcript("gB", "tB", "+", [(10050, 11000)])
    fixture.add_read("rA", "+", [(9000, 9900)])
    fixture.add_read("rB", "+", [(10050, 11000)])
    fixture.build()

    region = "{}+:1-9950".format(fixture.contig)
    # nothing crosses 9950, so margin 0 accepts it
    manifest = extractor.extract_partition(
        genome_fa=fixture.fasta,
        bam=fixture.bam,
        region=region,
        output_prefix=str(tmp_path / "m0"),
        gtf=fixture.gtf,
        margin=0,
    )
    assert manifest["counts"]["gtf_transcripts_emitted"] == 1

    with pytest.raises(extractor.ExtractionError) as excinfo:
        extractor.extract_partition(
            genome_fa=fixture.fasta,
            bam=fixture.bam,
            region=region,
            output_prefix=str(tmp_path / "m200"),
            gtf=fixture.gtf,
            margin=200,
        )
    message = str(excinfo.value)
    assert "within the 200 bp margin" in message
    assert "gA" in message or "gB" in message


def test_secondary_alignments_are_excluded_or_rejected_never_mishandled(tmp_path):
    fixture = Fixture(tmp_path, length=30000)
    fixture.add_transcript("gIn", "tIn", "+", [(12000, 12300), (13000, 13400)])
    fixture.add_read("rPrimary", "+", [(12000, 12300), (13000, 13400)])
    # a secondary alignment that straddles the right boundary: it must neither
    # block the cut nor reach the output
    fixture.add_read("rSecondary", "+", [(19900, 20100)], flag=256)
    # and a supplementary one, likewise
    fixture.add_read("rSupplementary", "+", [(19950, 20050)], flag=2048)
    fixture.build()

    region = "{}+:10001-20000".format(fixture.contig)
    manifest = extractor.extract_partition(
        genome_fa=fixture.fasta,
        bam=fixture.bam,
        region=region,
        output_prefix=str(tmp_path / "excl"),
        gtf=fixture.gtf,
        margin=0,
        secondary_alignments="exclude",
    )
    assert manifest["counts"]["nonprimary_excluded"] == 2
    assert set(_emitted_alignments(manifest["files"]["bam"])) == {"rPrimary"}

    with pytest.raises(extractor.ExtractionError) as excinfo:
        extractor.extract_partition(
            genome_fa=fixture.fasta,
            bam=fixture.bam,
            region=region,
            output_prefix=str(tmp_path / "rej"),
            gtf=fixture.gtf,
            margin=0,
            secondary_alignments="reject",
        )
    assert "secondary/supplementary" in str(excinfo.value)
    assert "rSecondary" in str(excinfo.value)


def test_gtf_line_without_identifying_attributes_is_rejected(tmp_path):
    fixture = Fixture(tmp_path, length=5000)
    fixture.add_read("r", "+", [(1000, 1200)])
    fixture.build()
    with open(fixture.gtf, "wt") as ofh:
        print(
            "\t".join(
                (fixture.contig, "test", "exon", "1000", "1200", ".", "+", ".", "note")
            ),
            file=ofh,
        )
    with pytest.raises(extractor.ExtractionError) as excinfo:
        extractor.load_gtf(fixture.gtf, fixture.contig, "+")
    assert "neither gene_id nor transcript_id" in str(excinfo.value)


def test_paired_alignment_is_rejected_rather_than_silently_rebased(tmp_path):
    fixture = Fixture(tmp_path, length=20000)
    fixture.add_read("rPaired", "+", [(12000, 12300)], flag=1)
    fixture.build()
    with pytest.raises(extractor.ExtractionError) as excinfo:
        extractor.extract_partition(
            genome_fa=fixture.fasta,
            bam=fixture.bam,
            region="{}+:10001-15000".format(fixture.contig),
            output_prefix=str(tmp_path / "paired"),
            margin=0,
        )
    assert "flagged paired" in str(excinfo.value)


def test_cli_extracts_a_chunk_and_names_the_reads_it_dropped(tmp_path, capsys):
    fixture = _overhang_fixture(tmp_path, "right")
    prefix = str(tmp_path / "cli")

    rc = extractor.main(
        [
            "--genome_fa",
            fixture.fasta,
            "--bam",
            fixture.bam,
            "--gtf",
            fixture.gtf,
            "--region",
            "{}+:10001-20000".format(fixture.contig),
            "--output_prefix",
            prefix,
            "--margin",
            "0",
        ]
    )
    assert rc == 0

    manifest = json.loads(Path(prefix + ".partition.json").read_text())
    assert manifest["counts"]["alignments_dropped_overhang"] == 1
    assert manifest["dropped_read_names"] == ["rSpanRight"]
    assert manifest["counts"]["alignments_emitted"] == 0
    assert Path(prefix + ".dropped_reads.txt").read_text().split() == ["rSpanRight"]

    # the drop is announced with its reason, not left to the manifest alone
    stderr = capsys.readouterr().err
    assert "rSpanRight" in stderr
    assert "dropped" in stderr
    assert "neither truncates nor widens" in stderr


def test_containment_assertions_catch_an_escaped_record_with_a_named_message(
    tmp_path,
):
    """Last-resort guards, exercised directly.

    Under the drop-overhang policy emission filters an overhanging alignment out
    before rebasing it, so these assertions are no longer reachable through
    ``extract_partition``. That is the point of keeping them: they are the
    backstop that makes "no out-of-range record is ever written" true by
    construction rather than by trusting the filter above them. Called directly,
    because a guard that cannot be reached cannot be tested through the front door.
    """

    fixture = _overhang_fixture(tmp_path, "right")
    region = extractor.parse_region("{}+:10001-20000".format(fixture.contig))
    offset = region.lend - 1
    mini_length = region.rend - region.lend + 1

    with pysam.AlignmentFile(fixture.bam, "rb") as bam:
        overhanging = next(
            aln for aln in bam.fetch(fixture.contig) if aln.query_name == "rSpanRight"
        )
        mini_header = extractor._mini_header(bam.header, fixture.contig, mini_length)

    with pytest.raises(extractor.ExtractionError) as excinfo:
        extractor._rebase_alignment(
            overhanging, offset, mini_header, fixture.contig, mini_length
        )
    message = str(excinfo.value)
    assert "rSpanRight" in message
    assert "mini contig of 1-10000" in message

    # The GTF guard cannot be reached through extract_partition: a gene's span is
    # the min/max over every line carrying its gene_id, so a contained gene has
    # only contained lines, and a straddling gene is not selected at all -- which
    # is exactly the whole-or-nothing property. Exercise the guard directly.
    line = "chrT\ttest\texon\t19900\t20500\t.\t+\t.\tgene_id \"gX\";\n"
    with pytest.raises(extractor.ExtractionError) as excinfo:
        extractor._rebase_gtf_line(line, 0, "chrT", 20000, "transcript tX")
    message = str(excinfo.value)
    assert "transcript tX" in message
    assert "mini contig of 1-20000" in message
    # and it accepts a feature that does fit, rebasing both coordinates
    assert extractor._rebase_gtf_line(
        line, 19899, "chrT", 20000, "transcript tX"
    ).split("\t")[3:5] == ["1", "601"]


@pytest.mark.parametrize(
    "region_str",
    ["chr20:100", "chr20+:200-100", "chr20+:0-100", "not a region", "chr20+:a-b"],
)
def test_malformed_region_strings_are_rejected(region_str):
    with pytest.raises(extractor.ExtractionError):
        extractor.parse_region(region_str)


def test_region_parsing_round_trips_strand_and_commas():
    region = extractor.parse_region("chr20-:1,000-2,000")
    assert region == extractor.Region("chr20", "-", 1000, 2000)
    assert extractor.parse_region("chrT:5-9").strand == ""


@pytest.mark.parametrize("strand", ["+", "-"])
def test_margin_guarantees_flanking_context_for_every_emitted_gtf_feature(
    tmp_path, strand
):
    """Every emitted annotated model sits at least ``margin`` bases inside both edges.

    That is what makes a rebased mini contig safe for the sequence-context windows
    LRAA evaluates around a model's ends -- the 40 bp polyA-signal window
    (`polyA_signal_window` [-40, -10]) and the 20 bp internal-priming window
    (`Util_funcs.check_internal_priming` check_length=20) -- and for the 50 bp
    read-enrichment window (`TSS_window_read_enrich_len`). It follows from the
    locus rule: a locus spanning [s, e] blocks every cut in
    [s - margin, e + margin - 1], so if both edges are admissible then
    s >= lend + margin and e <= rend - margin.

    ALIGNMENTS are deliberately excluded from this guarantee; see
    ``test_a_read_flush_against_a_boundary_is_still_emitted_at_a_positive_margin``.
    """

    margin = 200
    fixture = _multi_locus_fixture(tmp_path)

    checked = 0
    for index, (lend, rend) in enumerate(_MULTI_LOCUS_PARTITIONS):
        manifest = extractor.extract_partition(
            genome_fa=fixture.fasta,
            bam=fixture.bam,
            region="{}{}:{}-{}".format(fixture.contig, strand, lend, rend),
            output_prefix=str(tmp_path / "clear{}_{}".format(strand, index)),
            gtf=fixture.gtf,
            margin=margin,
        )
        mini_length = manifest["mini_length"]
        first_base_free = lend > 1
        last_base_free = rend < fixture.length
        for model in _emitted_gtf(manifest["files"]["gtf"]).values():
            if first_base_free:
                assert model["exons"][0][0] > margin
            if last_base_free:
                assert model["exons"][-1][1] <= mini_length - margin
            checked += 1
    assert checked > 0


def test_a_read_flush_against_a_boundary_is_still_emitted_at_a_positive_margin(
    tmp_path,
):
    """The margin governs loci, not alignments, and that is deliberate.

    Extending it to reads would mean dropping wholly contained alignments, making
    the extractor's dropped count exceed the spanning count the selector computed
    in advance. An accounting identity that holds is worth more here than a
    flanking guarantee for reads, whose ends are evidence rather than a model
    whose sequence context is being judged.
    """

    fixture = Fixture(tmp_path, length=30000)
    fixture.add_transcript("gFar", "tFar", "+", [(3000, 3400)])
    # ends on the chunk's very last base, and starts on its very first
    fixture.add_read("rFlushRight", "+", [(19000, 20000)])
    fixture.add_read("rFlushLeft", "+", [(10001, 11000)])
    fixture.build()

    manifest = extractor.extract_partition(
        genome_fa=fixture.fasta,
        bam=fixture.bam,
        region="{}+:10001-20000".format(fixture.contig),
        output_prefix=str(tmp_path / "flush"),
        gtf=fixture.gtf,
        margin=200,
    )

    emitted = _emitted_alignments(manifest["files"]["bam"])
    assert set(emitted) == {"rFlushLeft", "rFlushRight"}
    assert emitted["rFlushLeft"][0] == 1
    assert emitted["rFlushRight"][1] == manifest["mini_length"]
    assert manifest["counts"]["alignments_dropped_overhang"] == 0
    assert manifest["dropped_read_names"] == []
