#!/usr/bin/env python3
"""GeneBodyAssigner: gene choice and S/U/A class for reads isoform quant leaves unassigned."""

import io
import os
import sys

import pytest

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import LRAA_Globals
import GeneBodyAssigner as GBA


class Tx:
    """The four Transcript getters GeneBodyIndex reads."""

    def __init__(self, gene, exons, strand="+"):
        self._gene = gene
        self._exons = sorted(exons)
        self._strand = strand

    def get_output_gene_id(self):
        return self._gene

    def get_exon_segments(self):
        return list(self._exons)

    def get_strand(self):
        return self._strand


@pytest.fixture(autouse=True)
def defaults(monkeypatch):
    monkeypatch.setitem(LRAA_Globals.config, "gene_body_min_read_frac_in_span", 0.5)
    monkeypatch.setitem(LRAA_Globals.config, "gene_body_min_feature_bp", 20)
    monkeypatch.setitem(LRAA_Globals.config, "gene_body_ambiguous_policy", "drop")


def only_gene(tied, weights=None):
    assert tied is not None
    resolved, rule = GBA.GeneBodyIndex.resolve(tied, weights)
    assert rule != GBA.RULE_AMBIGUOUS and len(resolved) == 1, resolved
    t, frac = resolved[0]
    assert frac == 1.0
    return t.gene_id, t.splice_class


def ids_fracs(resolved):
    return [(t.gene_id, f) for t, f in resolved]


# G: exons 1000-1999, 5000-5999, 9000-9999 (introns 2000-4999, 6000-8999)
G = Tx("G", [(1000, 1999), (5000, 5999), (9000, 9999)])


def test_intronic_read_goes_to_its_gene_as_unspliced():
    idx = GBA.GeneBodyIndex([G], "+")
    assert only_gene(idx.assign([(3000, 3800)], [])) == ("G", "U")


def test_spliced_exonic_read_is_spliced():
    idx = GBA.GeneBodyIndex([G], "+")
    gene, cls = only_gene(idx.assign([(1500, 1999), (5000, 5400)], [(2000, 4999)]))
    assert (gene, cls) == ("G", "S")


def test_unspliced_exonic_read_is_ambiguous_class():
    idx = GBA.GeneBodyIndex([G], "+")
    assert only_gene(idx.assign([(5100, 5900)], [])) == ("G", "A")


def test_partially_spliced_read_with_retained_intron_is_unspliced():
    idx = GBA.GeneBodyIndex([G], "+")
    # spliced intron 1 out, intron 2 retained
    call = idx.assign([(1800, 1999), (5000, 9200)], [(2000, 4999)])
    assert only_gene(call) == ("G", "U")


def test_exon_end_wobble_below_min_feature_bp_stays_ambiguous_class():
    idx = GBA.GeneBodyIndex([G], "+")
    assert only_gene(idx.assign([(5100, 6010)], [])) == ("G", "A")


def test_antisense_gene_is_never_a_candidate():
    minus = Tx("M", [(2500, 4500)], strand="-")
    idx = GBA.GeneBodyIndex([G, minus], "+")
    assert set(idx.genes) == {"G"}
    # exonic for the antisense gene, intronic for G: G it is
    assert only_gene(idx.assign([(3000, 3800)], [])) == ("G", "U")


def test_mostly_intergenic_read_is_not_counted():
    idx = GBA.GeneBodyIndex([G], "+")
    # 300 of 1000 bases inside the span
    assert idx.assign([(9700, 10699)], []) is None


def test_downstream_overhang_within_threshold_is_counted():
    idx = GBA.GeneBodyIndex([G], "+")
    # 700 of 1000 aligned bases in the span: last exon + 300 bp unannotated 3' extension
    assert only_gene(idx.assign([(9300, 10299)], [])) == ("G", "A")


# Nested: N sits inside G's second intron, same strand.
N = Tx("N", [(6500, 6799), (7500, 7799)])


def test_nested_gene_exon_beats_host_intron():
    idx = GBA.GeneBodyIndex([G, N], "+")
    assert only_gene(idx.assign([(6550, 6790)], [])) == ("N", "A")


def test_read_in_nested_gene_intron_goes_to_innermost_gene():
    idx = GBA.GeneBodyIndex([G, N], "+")
    # inside N's intron and inside G's intron: both contain it, N is smaller
    assert only_gene(idx.assign([(6900, 7400)], [])) == ("N", "U")


def test_host_intron_read_outside_nested_gene_goes_to_host():
    idx = GBA.GeneBodyIndex([G, N], "+")
    assert only_gene(idx.assign([(8000, 8800)], [])) == ("G", "U")


def test_shared_annotated_junction_outranks_exonic_bases():
    # H has a long exon over N's whole first-intron region; N's junction is annotated
    H = Tx("H", [(6000, 8000)])
    idx = GBA.GeneBodyIndex([H, N], "+")
    call = idx.assign([(6600, 6799), (7500, 7600)], [(6800, 7499)])
    assert only_gene(call) == ("N", "S")


def test_readthrough_gene_loses_constituent_intronic_reads():
    A = Tx("A", [(1000, 1999), (3000, 3999)])
    B = Tx("B", [(6000, 6999), (8000, 8999)])
    RT = Tx("A-B", [(1000, 1999), (3000, 3999), (6000, 6999), (8000, 8999)])
    idx = GBA.GeneBodyIndex([A, B, RT], "+")
    assert only_gene(idx.assign([(2200, 2800)], [])) == ("A", "U")
    assert only_gene(idx.assign([(7200, 7800)], [])) == ("B", "U")
    # in the A-B intergenic gap only the readthrough gene spans it
    assert only_gene(idx.assign([(4500, 5500)], [])) == ("A-B", "U")


def test_tandem_read_into_next_gene_goes_to_gene_holding_it():
    A = Tx("A", [(1000, 1999), (3000, 3999)])
    B = Tx("B", [(4100, 4999), (6000, 6999)])
    idx = GBA.GeneBodyIndex([A, B], "+")
    # 600 bp in A's last exon, 100 bp gap, 100 bp in B's first exon
    assert only_gene(idx.assign([(3400, 4199)], []))[0] == "A"


def test_exonic_graze_below_min_does_not_beat_containing_intron():
    # G2 starts inside G's first intron; the read grazes G2's first exon by 10 bp
    G2 = Tx("G2", [(3790, 4200), (4800, 4900)])
    idx = GBA.GeneBodyIndex([G, G2], "+")
    call = idx.assign([(3000, 3799)], [])
    assert call is not None
    # G2 span holds only 10 bp of the read: not a candidate at all
    assert only_gene(call) == ("G", "U")


def test_innermost_span_breaks_containment_tie():
    A = Tx("A", [(1000, 1099), (2900, 2999)])
    B = Tx("B", [(2000, 2099), (5000, 5099)])
    idx = GBA.GeneBodyIndex([A, B], "+")
    # intronic for both, contained in both, 3' end in both: the smaller span (A) wins
    assert only_gene(idx.assign([(2150, 2850)], [])) == ("A", "U")
    # mostly beyond A's span: A is not a candidate, B is
    assert only_gene(idx.assign([(2950, 3600)], [])) == ("B", "U")


def test_three_prime_end_decides_when_neither_contains():
    # exon overlaps below gene_body_min_feature_bp (11 bp) count as no exonic evidence
    A = Tx("A", [(1000, 1099), (2990, 3000)])
    B = Tx("B", [(2500, 2510), (5000, 5099)])
    idx_p = GBA.GeneBodyIndex([A, B], "+")
    # read 2300..3300: neither span contains it; '+' 3' end (3300) lies in B only
    assert only_gene(idx_p.assign([(2300, 3300)], []))[0] == "B"
    Am = Tx("Am", [(1000, 1099), (2990, 3000)], strand="-")
    Bm = Tx("Bm", [(2500, 2510), (5000, 5099)], strand="-")
    idx_m = GBA.GeneBodyIndex([Am, Bm], "-")
    # same geometry on '-': the 3' end is the left end (2300), in Am only
    assert only_gene(idx_m.assign([(2300, 3300)], [], "-"))[0] == "Am"
    # containment outranks the 3' end: B2's span holds the whole read
    B2 = Tx("B2", [(2000, 2010), (5000, 5099)])
    idx2 = GBA.GeneBodyIndex([A, B2], "+")
    assert only_gene(idx2.assign([(2400, 3200)], []))[0] == "B2"


def test_unresolvable_tie_drop_and_split(monkeypatch):
    X = Tx("X", [(1000, 1999), (5000, 5999)])
    Y = Tx("Y", [(1000, 1999), (5000, 5999)])
    idx = GBA.GeneBodyIndex([Y, X], "+")
    tied = idx.assign([(3000, 3500)], [])
    resolved, rule = GBA.GeneBodyIndex.resolve(tied)
    assert rule == GBA.RULE_AMBIGUOUS
    assert ids_fracs(resolved) == [("X", 0.0), ("Y", 0.0)]
    # weights with nothing on either gene change nothing
    assert GBA.GeneBodyIndex.resolve(tied, {"Z": 5})[1] == GBA.RULE_AMBIGUOUS
    monkeypatch.setitem(LRAA_Globals.config, "gene_body_ambiguous_policy", "split")
    resolved, rule = GBA.GeneBodyIndex.resolve(tied)
    assert rule == GBA.RULE_AMBIGUOUS and ids_fracs(resolved) == [("X", 0.5), ("Y", 0.5)]
    # deterministic: same order on repeat, sorted by gene id
    assert [t.gene_id for t in idx.assign([(3000, 3500)], [])] == ["X", "Y"]


def test_tie_weighting_dominance_and_split(monkeypatch):
    # The ALS chr8 geometry: an annotated readthrough gene (MRPS28-TPD52) with a SMALLER
    # span than the gene whose 3' region it overlaps. Innermost alone hands it TPD52's
    # intronic reads; weighting by isoform-assigned reads gives them back.
    TPD52 = Tx("TPD52", [(34000, 34999), (60000, 60999), (90000, 90999)])
    RT = Tx("RT", [(20000, 20999), (60000, 60999)])
    idx = GBA.GeneBodyIndex([TPD52, RT], "+")
    tied = idx.assign([(40000, 41000)], [])
    assert [t.gene_id for t in tied] == ["RT", "TPD52"]
    assert only_gene(tied) == ("RT", "U")  # innermost, no weights
    assert GBA.GeneBodyIndex.resolve(tied)[1] == GBA.RULE_INNERMOST
    # dominance: >= 0.9 of the weight takes the whole read
    resolved, rule = GBA.GeneBodyIndex.resolve(tied, {"TPD52": 1536.0, "RT": 0.0})
    assert rule == GBA.RULE_DOMINANCE and ids_fracs(resolved) == [("TPD52", 1.0)]
    resolved, rule = GBA.GeneBodyIndex.resolve(tied, {"TPD52": 95.0, "RT": 5.0})
    assert rule == GBA.RULE_DOMINANCE and ids_fracs(resolved) == [("TPD52", 1.0)]
    resolved, rule = GBA.GeneBodyIndex.resolve(tied, {"TPD52": 90.0, "RT": 10.0})
    assert rule == GBA.RULE_DOMINANCE  # exactly at the cut
    # below the cut: fractional by weight
    resolved, rule = GBA.GeneBodyIndex.resolve(tied, {"TPD52": 3.0, "RT": 1.0})
    assert rule == GBA.RULE_WEIGHT_SPLIT
    assert ids_fracs(resolved) == [("RT", 0.25), ("TPD52", 0.75)]
    # the cut is a config key
    monkeypatch.setitem(LRAA_Globals.config, "gene_body_dominance_frac", 0.7)
    resolved, rule = GBA.GeneBodyIndex.resolve(tied, {"TPD52": 3.0, "RT": 1.0})
    assert rule == GBA.RULE_DOMINANCE and ids_fracs(resolved) == [("TPD52", 1.0)]
    monkeypatch.setitem(LRAA_Globals.config, "gene_body_dominance_frac", 1.0)
    resolved, rule = GBA.GeneBodyIndex.resolve(tied, {"TPD52": 95.0, "RT": 5.0})
    assert rule == GBA.RULE_WEIGHT_SPLIT


def test_isoform_weight_outranks_three_prime_end():
    # neither span contains the read; its '+' 3' end lies in B only, but A has the
    # isoform-assigned reads: weighting comes first
    A = Tx("A", [(1000, 1099), (2990, 3000)])
    B = Tx("B", [(2500, 2510), (5000, 5099)])
    idx = GBA.GeneBodyIndex([A, B], "+")
    tied = idx.assign([(2300, 3300)], [])
    assert [t.gene_id for t in tied] == ["A", "B"]
    assert only_gene(tied) == ("B", "U")  # 3' end, no weights
    assert GBA.GeneBodyIndex.resolve(tied)[1] == GBA.RULE_THREE_PRIME
    assert only_gene(tied, {"A": 50.0}) == ("A", "U")


class TaggedTx(Tx):
    def __init__(self, gene, exons, tags, strand="+"):
        super().__init__(gene, exons, strand)
        self._tags = frozenset(tags)

    def get_gtf_tags(self):
        return self._tags


def test_readthrough_genes_excluded_as_candidates_by_default(monkeypatch):
    TPD52 = Tx("TPD52", [(34000, 34999), (60000, 60999), (90000, 90999)])
    RT = TaggedTx("RT", [(20000, 20999), (60000, 60999)], ["basic", "readthrough_transcript"])
    idx = GBA.GeneBodyIndex([TPD52, RT], "+")
    assert idx.genes["RT"].readthrough and not idx.genes["TPD52"].readthrough
    # TPD52 intron: no tie at all, RT is not a candidate
    assert [t.gene_id for t in idx.assign([(40000, 41000)], [])] == ["TPD52"]
    # RT-only region (its own intron upstream of TPD52): no gene
    assert idx.assign([(25000, 26000)], []) is None
    # isoform rows still classify against RT
    cls, ev = idx.classify_for_gene("RT", [(20100, 20900)], [])
    assert cls == "A" and ev is not None
    monkeypatch.setitem(LRAA_Globals.config, "gene_body_exclude_readthrough_genes", False)
    assert [t.gene_id for t in idx.assign([(40000, 41000)], [])] == ["RT", "TPD52"]


def test_gene_is_readthrough_only_if_every_transcript_is():
    mixed = [
        TaggedTx("M", [(1000, 1999), (5000, 5999)], ["readthrough_transcript"]),
        TaggedTx("M", [(1000, 1999), (7000, 7999)], ["basic"]),
    ]
    assert not GBA.GeneBodyIndex(mixed, "+").genes["M"].readthrough
    gene_tag = TaggedTx("Q", [(1000, 1999)], ["readthrough_gene"])
    assert GBA.GeneBodyIndex([gene_tag], "+").genes["Q"].readthrough


def test_gtf_parser_keeps_repeated_tags(tmp_path):
    from Transcript import GTF_contig_to_transcripts

    gtf = tmp_path / "t.gtf"
    attrs = (
        'gene_id "G1"; transcript_id "T1"; gene_name "RTG"; tag "basic"; '
        'tag "readthrough_transcript"; tag "CCDS";'
    )
    lines = [
        "\t".join(["chr1", "HAVANA", "exon", "100", "200", ".", "+", ".", attrs]),
        "\t".join(["chr1", "HAVANA", "exon", "300", "400", ".", "+", ".", attrs]),
        "\t".join(["chr1", "HAVANA", "exon", "100", "200", ".", "+", ".",
                   'gene_id "G2"; transcript_id "T2";']),
    ]
    gtf.write_text("\n".join(lines) + "\n")
    parsed = GTF_contig_to_transcripts.parse_GTF_to_Transcripts(str(gtf))
    txs = {t.get_transcript_id(): t for t in parsed["chr1"]}
    assert txs["T1"].get_gtf_tags() == frozenset({"basic", "readthrough_transcript", "CCDS"})
    assert txs["T2"].get_gtf_tags() == frozenset()
    # the parsed meta holds what it held before (the last tag wins there)
    assert txs["T1"].get_meta()["tag"] == "CCDS"


def test_segments_from_cigar():
    # 10S 100M 5D 50M 1000N 20M 3I 30M 5S starting at 0-based 999
    cig = [(4, 10), (0, 100), (2, 5), (0, 50), (3, 1000), (0, 20), (1, 3), (0, 30), (4, 5)]
    segs, introns = GBA.segments_from_cigar(999, cig)
    assert segs == [(1000, 1154), (2155, 2204)]
    assert introns == [(1155, 2154)]


def test_recorder_isoform_rows_summed_per_gene_and_unassigned():
    idx = GBA.GeneBodyIndex([G, N], "+")
    fh = io.StringIO()
    rec = GBA.GeneBodyRecorder(idx, fh)
    rows = [
        ("G", "t1", "h", "3", "mp", 0.25, 1.0, "0"),
        ("G", "t2", "h", "3", "mp", 0.5, 1.0, "0"),
        ("N", "t3", "h", "2", "mp", 0.25, 1.0, "0"),
    ]
    rec.record_isoform("bc^umi^r1", [(5100, 5900)], [], rows)
    rec.record_unassigned("bc^umi^r2", [(3000, 3800)], [], "no_path")
    rec.record_unassigned("bc^umi^r3", [(20000, 20500)], [], "no_path")
    lines = [l.split("\t") for l in fh.getvalue().splitlines()]
    assert [(l[0], l[1], l[2], l[3], l[4], l[5]) for l in lines] == [
        ("G", "bc^umi^r1", "0.750000", "A", "isoform", "2"),
        ("N", "bc^umi^r1", "0.250000", "A", "isoform", "2"),
        ("G", "bc^umi^r2", "1.000000", "U", "gene_body:no_path", "1"),
    ]
    assert rec.counts["reads_no_gene"] == 1
    assert len(lines[0]) == len(GBA.GENE_BODY_TRACKING_HEADER)


def test_recorder_spills_ties_and_resolves_with_full_stream_counts(tmp_path):
    X = Tx("X", [(1000, 1999), (5000, 5999)])
    Y = Tx("Y", [(1000, 1999), (5000, 5999)])
    idx = GBA.GeneBodyIndex([X, Y], "+")
    fh = io.StringIO()
    rec = GBA.GeneBodyRecorder(idx, fh, spill_dir=str(tmp_path))
    # the tied read arrives BEFORE the isoform reads that decide it
    rec.record_unassigned("bc^u^tied", [(3000, 3500)], [], "unassignable")
    assert fh.getvalue() == ""
    rec.record_isoform("bc^u^r1", [(1100, 1900)], [], [("Y", "ty", "h", "2", "mp", 1.0, 1.0, "1")])
    rec.finish()
    rows = [l.split("\t") for l in fh.getvalue().splitlines()]
    assert [(r[0], r[1], r[2], r[3], r[4], r[5], r[7], r[8]) for r in rows] == [
        ("Y", "bc^u^r1", "1.000000", "A", "isoform", "1", "801", "0"),
        ("Y", "bc^u^tied", "1.000000", "U", "gene_body:unassignable", "1", "0", "501"),
    ]
    assert rec.counts["reads_gene_tied"] == 1
    assert rec.counts["tie_" + GBA.RULE_DOMINANCE] == 1


def test_chunked_merge_of_gene_body_sidecars(tmp_path):
    import gzip
    import ChunkedRun

    header = "\t".join(GBA.GENE_BODY_TRACKING_HEADER) + "\n"
    units = []
    for k in range(2):
        prefix = str(tmp_path / f"u{k}")
        with gzip.open(prefix + ChunkedRun.GENE_BODY_TRACKING_SUFFIX, "wt") as fh:
            fh.write("# LRAA version test\n" + header)
            fh.write(f"G{k}\tbc^u^r{k}\t1.000000\tU\tgene_body:no_path\t1\t10\t0\t10\n")
        units.append({"quant_prefix": prefix})
    assert ChunkedRun.GENE_BODY_TRACKING_SUFFIX == GBA.gene_body_tracking_suffix()
    out = ChunkedRun.merge_gene_body_tracking(str(tmp_path), units)
    lines = gzip.open(out, "rt").read().splitlines()
    assert lines[0].startswith("#") and lines[1] + "\n" == header
    assert [l.split("\t")[0] for l in lines[2:]] == ["G0", "G1"]
    absent = {"quant_prefix": str(tmp_path / "absent")}
    # no unit wrote one: nothing to publish; only some did: refused
    assert ChunkedRun.merge_gene_body_tracking(str(tmp_path), [absent]) is None
    with pytest.raises(ChunkedRun.PipelineError):
        ChunkedRun.merge_gene_body_tracking(str(tmp_path), units + [absent])
