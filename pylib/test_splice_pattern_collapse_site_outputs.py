"""The collapse must report WHICH site each count and PAS motif belongs to.

A collapsed model spans every terminal variant its members had, so `PolyA_sites
"1700,1900,2100"` alone says three sites exist and nothing about which carries 60 reads
or which is an oligo-dT artifact. The per-site lists are only meaningful if index i of
every list describes the same site, and the site beds are only honest if support is
deduplicated rather than summed -- every transcript ending at a site copies the support
of one splice-graph node, so summing multiplies a site by however many isoforms share it.

PAS metadata belongs to the 3' end alone. Attaching it to a TSS looks plausible and is
wrong: MEASURED on testing/single_contig, sixteen transcripts share the TSS at
minigenome:2031 with identical support 98 and split 10/6 on PAS_offset (-24 vs -18),
purely because they have different PolyA ends.
"""
import os
import subprocess
import sys

import pytest

sys.path.insert(0, os.path.dirname(os.path.realpath(__file__)))

import SplicePatternCollapse
from SplicePatternCollapse import GeneConflictError, SiteSupportConflictError

SCRIPT = os.path.join(
    os.path.dirname(os.path.realpath(__file__)),
    "..",
    "util",
    "collapse_LRAA_GTF_by_splice_pattern.py",
)


def _transcript(
    transcript_id,
    exons,
    strand="+",
    gene_id="MYGENE^ENSG001",
    TSS=None,
    PolyA=None,
    PAS=None,
    PAS_offset=None,
    InternalPriming=None,
    contig="chr1",
):
    attrs = 'gene_id "{}"; transcript_id "{}";'.format(gene_id, transcript_id)
    for key, value in (
        ("TSS", "True" if TSS is not None else None),
        ("PolyA", "True" if PolyA is not None else None),
        ("TSS_read_count", TSS),
        ("PolyA_read_count", PolyA),
        ("PAS", PAS),
        ("PAS_offset", PAS_offset),
        ("InternalPriming", InternalPriming),
    ):
        if value is not None:
            attrs += ' {} "{}";'.format(key, value)

    rows = [
        "\t".join(
            [
                contig,
                "LRAA",
                "transcript",
                str(exons[0][0]),
                str(exons[-1][1]),
                ".",
                strand,
                ".",
                attrs,
            ]
        )
    ]
    exon_attrs = 'gene_id "{}"; transcript_id "{}";'.format(gene_id, transcript_id)
    for lend, rend in exons:
        rows.append(
            "\t".join(
                [contig, "LRAA", "exon", str(lend), str(rend), ".", strand, ".", exon_attrs]
            )
        )
    return rows


def _write_gtf(path, transcript_rows):
    with open(path, "wt") as fh:
        fh.write("# LRAA version v-test\n")
        for rows in transcript_rows:
            fh.write("\n".join(rows) + "\n")
    return str(path)


# Three isoforms of one intron chain, differing only at their termini: two TSS and three
# PolyA sites, each with its own support and PAS annotation.
def _alt_termini_gtf(path):
    return _write_gtf(
        path,
        [
            _transcript(
                "t:iso-1",
                [(1000, 1200), (1500, 1700)],
                TSS=40,
                PolyA=25,
                PAS="AATAAA",
                PAS_offset=-21,
                InternalPriming="False",
            ),
            _transcript(
                "t:iso-2",
                [(900, 1200), (1500, 1900)],
                TSS=12,
                PolyA=60,
                PAS="ATTAAA",
                PAS_offset=-16,
                InternalPriming="False",
            ),
            _transcript(
                "t:iso-3",
                [(1000, 1200), (1500, 2100)],
                TSS=40,
                PolyA=8,
                PAS="none",
                InternalPriming="True",
            ),
        ],
    )


def _attrs(gtf_path, predicate=lambda attrs: True):
    found = []
    for line in open(gtf_path):
        if line.startswith("#"):
            continue
        cols = line.rstrip("\n").split("\t")
        if len(cols) < 9 or cols[2] != "transcript":
            continue
        attrs = dict(
            pair.strip().split(' "', 1)
            for pair in cols[8].rstrip(";").split(";")
            if ' "' in pair
        )
        attrs = {k: v.rstrip('"') for k, v in attrs.items()}
        if predicate(attrs):
            found.append(attrs)
    return found


def _bed_rows(path):
    rows = []
    for line in open(path):
        if line.startswith("#"):
            continue
        rows.append(line.rstrip("\n").split("\t"))
    return rows


def test_per_site_lists_stay_aligned_with_their_coordinates(tmp_path):
    out = str(tmp_path / "collapsed.gtf")
    SplicePatternCollapse.collapse_gtf(_alt_termini_gtf(tmp_path / "in.gtf"), out)

    merged = _attrs(out, lambda a: "merged_isoforms_shared_splice_pattern" in a)
    assert len(merged) == 1
    merged = merged[0]

    # Index i of every list describes site i.  Note the TSS order is NOT input order:
    # iso-1 comes first in the file with support 40 at coordinate 1000, but 900 sorts
    # ahead of it, so a list built from a separate traversal would invert these.
    assert merged["TSS_sites"] == "900,1000"
    assert merged["TSS_site_support"] == "12,40"

    assert merged["PolyA_sites"] == "1700,1900,2100"
    assert merged["PolyA_site_support"] == "25,60,8"
    assert merged["PolyA_site_PAS"] == "AATAAA,ATTAAA,none"
    # "." holds position 2 rather than the list being short, which would silently shift
    # every later field by one.
    assert merged["PolyA_site_PAS_offset"] == "-21,-16,."
    assert merged["PolyA_site_internal_priming"] == "False,False,True"


def test_pas_annotation_is_not_attached_to_TSS_sites(tmp_path):
    out = str(tmp_path / "collapsed.gtf")
    written = SplicePatternCollapse.collapse_gtf(
        _alt_termini_gtf(tmp_path / "in.gtf"), out
    )

    merged = _attrs(out, lambda a: "merged_isoforms_shared_splice_pattern" in a)[0]
    assert not [key for key in merged if key.startswith("TSS_site_PAS")]
    assert "internal_priming" not in " ".join(
        key for key in merged if key.startswith("TSS")
    )

    # and the TSS bed carries no motif columns either
    header = [
        line for line in open(written["TSS_bed"]) if line.startswith("#chrom")
    ][0]
    assert "pas" not in header
    assert header.rstrip("\n").split("\t")[-1] == "transcript_ids"


def test_site_support_is_deduplicated_rather_than_summed(tmp_path):
    """Two transcripts sharing a TSS report one site at its own support, not twice it."""
    out = str(tmp_path / "collapsed.gtf")
    written = SplicePatternCollapse.collapse_gtf(
        _alt_termini_gtf(tmp_path / "in.gtf"), out
    )

    tss = _bed_rows(written["TSS_bed"])
    assert len(tss) == 2

    shared = [row for row in tss if row[2] == "1000"][0]
    assert shared[6] == "40"  # not 80: iso-1 and iso-3 copy one node's support
    assert shared[7] == "2"
    assert shared[8] == "t:iso-1,t:iso-3"

    # one base, half open
    assert (shared[1], shared[2]) == ("999", "1000")
    assert shared[3] == "TSS:chr1:1000:+"


def test_disagreeing_support_at_one_site_is_refused(tmp_path):
    """Boundary support comes from one splice-graph node, so it cannot disagree."""
    gtf = _write_gtf(
        tmp_path / "in.gtf",
        [
            _transcript("t:iso-1", [(1000, 1200), (1500, 1700)], TSS=40),
            # same TSS coordinate, different count
            _transcript("t:iso-2", [(1000, 1200), (1500, 1900)], TSS=17),
        ],
    )
    with pytest.raises(SiteSupportConflictError) as raised:
        SplicePatternCollapse.collapse_gtf(gtf, str(tmp_path / "collapsed.gtf"))

    message = str(raised.value)
    assert "chr1:1000" in message
    assert "t:iso-1" in message and "t:iso-2" in message


def test_minus_strand_sites_are_ordered_and_attributed_by_coordinate(tmp_path):
    """On '-' the TSS is the high coordinate, and lists stay genomic ascending."""
    gtf = _write_gtf(
        tmp_path / "in.gtf",
        [
            _transcript(
                "t:iso-1",
                [(1000, 1200), (1500, 1700)],
                strand="-",
                TSS=30,
                PolyA=11,
                PAS="AATAAA",
                PAS_offset=-19,
                InternalPriming="False",
            ),
            _transcript(
                "t:iso-2",
                [(900, 1200), (1500, 1900)],
                strand="-",
                TSS=5,
                PolyA=44,
                PAS="none",
                InternalPriming="True",
            ),
        ],
    )
    out = str(tmp_path / "collapsed.gtf")
    written = SplicePatternCollapse.collapse_gtf(gtf, out)

    merged = _attrs(out, lambda a: "merged_isoforms_shared_splice_pattern" in a)[0]
    # 1700 is iso-1's 5' end, 1900 is iso-2's; ascending puts the 3'-most TSS first.
    assert merged["TSS_sites"] == "1700,1900"
    assert merged["TSS_site_support"] == "30,5"
    assert merged["PolyA_sites"] == "900,1000"
    assert merged["PolyA_site_support"] == "44,11"
    assert merged["PolyA_site_PAS"] == "none,AATAAA"
    assert merged["PolyA_site_internal_priming"] == "True,False"

    polyA = {row[2]: row for row in _bed_rows(written["PolyA_bed"])}
    assert polyA["900"][11] == "True"
    assert polyA["1000"][9:12] == ["AATAAA", "-19", "False"]
    assert all(row[5] == "-" for row in polyA.values())


def test_single_member_groups_keep_their_scalar_attributes(tmp_path):
    """Nothing to attribute when a group has one member, so nothing is restated."""
    gtf = _write_gtf(
        tmp_path / "in.gtf",
        [
            _transcript(
                "t:solo",
                [(1000, 1200), (1500, 1700)],
                TSS=7,
                PolyA=9,
                PAS="AATAAA",
                PAS_offset=-30,
                InternalPriming="False",
            )
        ],
    )
    out = str(tmp_path / "collapsed.gtf")
    SplicePatternCollapse.collapse_gtf(gtf, out)

    solo = _attrs(out)[0]
    assert solo["TSS_read_count"] == "7"
    assert solo["PolyA_read_count"] == "9"
    assert solo["PAS"] == "AATAAA" and solo["PAS_offset"] == "-30"
    assert "TSS_site_support" not in solo


def test_beds_exist_and_are_empty_when_no_boundaries_were_called(tmp_path):
    """--no_infer_TSS/--no_infer_PolyA must not turn a declared output into a gap."""
    gtf = _write_gtf(
        tmp_path / "in.gtf",
        [_transcript("t:iso-1", [(1000, 1200), (1500, 1700)])],
    )
    written = SplicePatternCollapse.collapse_gtf(gtf, str(tmp_path / "collapsed.gtf"))

    for key in ("TSS_bed", "PolyA_bed"):
        assert os.path.exists(written[key])
        assert _bed_rows(written[key]) == []
        assert open(written[key]).readline().startswith("#")


# One intron chain under two gene_ids: collapsing would mint one transcript_id twice.
CONFLICTING = [
    _transcript("t:a", [(1000, 1200), (1500, 1700)], gene_id="GENE_A", TSS=5),
    _transcript("t:b", [(1000, 1200), (1500, 1700)], gene_id="GENE_B", PolyA=6),
]


def test_fatal_conflict_writes_only_the_conflicts_report(tmp_path):
    """The standalone contract: refuse, leaving nothing to mistake for a collapse."""
    gtf = _write_gtf(tmp_path / "in.gtf", CONFLICTING)
    out = str(tmp_path / "collapsed.gtf")

    with pytest.raises(GeneConflictError):
        SplicePatternCollapse.collapse_gtf(gtf, out)

    assert os.path.exists(out + ".gene_conflicts.tsv")
    assert not os.path.exists(out)
    assert not os.path.exists(out + ".isoform_merge_report.tsv")
    assert not os.path.exists(out + ".TSS.bed")
    assert not os.path.exists(out + ".PolyA.bed")


def test_nonfatal_conflict_still_writes_the_beds(tmp_path):
    """The integrated contract: the beds describe the input, not the collapse.

    A WDL task that failed here would publish neither them nor the conflicts report.
    """
    gtf = _write_gtf(tmp_path / "in.gtf", CONFLICTING)
    out = str(tmp_path / "collapsed.gtf")

    written = SplicePatternCollapse.collapse_gtf(
        gtf, out, nonfatal_gene_conflicts=True
    )

    assert set(written) == {"gene_conflicts_report", "TSS_bed", "PolyA_bed"}
    assert not os.path.exists(out)
    assert len(_bed_rows(written["TSS_bed"])) == 1
    assert len(_bed_rows(written["PolyA_bed"])) == 1
    assert len(open(written["gene_conflicts_report"]).readlines()) == 3  # header + 2


def test_cli_nonfatal_flag_exits_zero_and_reports(tmp_path):
    gtf = _write_gtf(tmp_path / "in.gtf", CONFLICTING)
    out = str(tmp_path / "collapsed.gtf")

    result = subprocess.run(
        [
            sys.executable,
            SCRIPT,
            "--gtf",
            gtf,
            "--output_gtf",
            out,
            "--nonfatal_gene_conflicts",
        ],
        capture_output=True,
        text=True,
    )
    assert result.returncode == 0
    assert "cross-gene splice patterns" in result.stderr
    assert os.path.exists(out + ".TSS.bed")
    assert not os.path.exists(out)


def test_collapsed_gtf_keeps_the_source_provenance(tmp_path):
    """Reserializing transcripts drops the header; the version must survive it."""
    out = str(tmp_path / "collapsed.gtf")
    SplicePatternCollapse.collapse_gtf(_alt_termini_gtf(tmp_path / "in.gtf"), out)

    comments = [line.rstrip("\n") for line in open(out) if line.startswith("#")]
    assert "# LRAA version v-test" in comments
    assert any("splice-pattern collapsed from in.gtf" in line for line in comments)
