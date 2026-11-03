#!/usr/bin/env python3
"""Reference identity in ref-guided runs (v0.50.0): the reference-splice-pattern reprieve
(last standing carrier of a reference chain, >= 1 full-bam read carrying it), single-exon
ref_supported, guide GTF class, and the splice_hash / ref_splice_pattern / ref_FSM_reads /
ref_supported GTF attributes."""

import re

import pysam
import pytest

import LRAA_Globals
import ReferenceSupport
import SplicePatternCollapse
import TranscriptFiltering
import Util_funcs
from Transcript import Transcript, GTF_contig_to_transcripts


def _vals(gtf_text, key):
    return re.findall(r'\b{} "([^"]*)"'.format(key), gtf_text.splitlines()[0])


def _parse(tmp_path, lines, name="in.gtf"):
    path = tmp_path / name
    path.write_text("\n".join(lines) + "\n")
    return GTF_contig_to_transcripts.parse_GTF_to_Transcripts(str(path))["chr1"]


def _guide_rows(source, extra="", tid="t1"):
    """One two-exon transcript (introns 201-299) as GTF rows from `source`."""
    return [
        f'chr1\t{source}\ttranscript\t100\t500\t.\t+\t.\tgene_id "g1"; transcript_id "{tid}";{extra}',
        f'chr1\t{source}\texon\t100\t200\t.\t+\t.\tgene_id "g1"; transcript_id "{tid}";',
        f'chr1\t{source}\texon\t300\t500\t.\t+\t.\tgene_id "g1"; transcript_id "{tid}";',
    ]


def _se_rows(source, extra="", tid="s1", lend=600, rend=700):
    return [
        f'chr1\t{source}\ttranscript\t{lend}\t{rend}\t.\t+\t.\tgene_id "g2"; transcript_id "{tid}";{extra}',
        f'chr1\t{source}\texon\t{lend}\t{rend}\t.\t+\t.\tgene_id "g2"; transcript_id "{tid}";',
    ]


def _model(tid="m", novel=False, ref_FSM_reads=1):
    t = Transcript("chr1", [[100, 200], [300, 500]], "+")
    t.set_transcript_id(tid)
    t.set_simple_path(["E:1", "I:1", "E:2"])
    t.set_has_novel_splice_pattern(novel)
    t.set_ref_FSM_reads(ref_FSM_reads)
    return t


def _se(tid, lend, rend, strand="+", guide_class=None, ref_supported=None):
    t = Transcript("chr1", [[lend, rend]], strand)
    t.set_transcript_id(tid)
    t._guide_class = guide_class
    t.set_ref_supported(ref_supported)
    return t


# ---- defaults and GTF attributes ------------------------------------------------------


def test_rules_on_by_default():
    assert LRAA_Globals.config["ref_splice_pattern_reprieve"] is True
    assert LRAA_Globals.config["ref_splice_pattern_reprieve_min_FSM_reads"] == 1
    assert LRAA_Globals.config["ref_supported_monoexonic"] is True
    assert LRAA_Globals.config["ref_supported_monoexonic_min_reciprocal_overlap"] == 0.8
    assert LRAA_Globals.config["allow_unflagged_LRAA_guide"] is False


def test_attributes_written_when_recorded():
    t = _model(ref_FSM_reads=7)
    assert _vals(t.to_GTF_format(), "ref_splice_pattern") == []  # not recorded: none
    t.set_ref_splice_pattern(True)
    text = t.to_GTF_format()
    assert _vals(text, "splice_hash") == [Util_funcs.get_hash_code(t.get_introns_string())]
    assert _vals(text, "ref_splice_pattern") == ["True"]
    assert _vals(text, "ref_FSM_reads") == ["7"]
    t.set_ref_splice_pattern(False)  # a novel chain carries no count
    text = t.to_GTF_format()
    assert _vals(text, "ref_splice_pattern") == ["False"] and _vals(text, "ref_FSM_reads") == []


def test_splice_hash_matches_read_hash():
    """The GTF splice_hash is the read-side hash the full-bam count matches."""
    t = _model()
    assert t.get_splice_hash() == Util_funcs.get_hash_code(t.get_introns_string())
    assert _se("s", 1, 100).get_splice_hash() is None


def test_single_exon_attribute():
    t = _se("s", 600, 700, ref_supported=True)
    t.set_simple_path(["E:9"])
    text = t.to_GTF_format()
    assert _vals(text, "ref_supported") == ["True"]
    assert _vals(text, "splice_hash") == [] and _vals(text, "ref_splice_pattern") == []


def test_attributes_round_trip_once(tmp_path):
    lines = _guide_rows(
        "LRAA", ' splice_hash "x"; ref_splice_pattern "True"; ref_FSM_reads "4";'
    ) + _se_rows("LRAA", ' ref_supported "False";')
    t, se = sorted(_parse(tmp_path, lines), key=lambda x: x.get_transcript_id())[::-1]
    assert t.get_ref_splice_pattern() is True and t.get_ref_FSM_reads() == 4
    assert se.get_ref_supported() is False
    t.set_simple_path(["E:1", "I:1", "E:2"])
    text = t.to_GTF_format()
    # re-emitted from the model, once; splice_hash recomputed rather than copied
    assert _vals(text, "ref_splice_pattern") == ["True"]
    assert _vals(text, "ref_FSM_reads") == ["4"]
    assert _vals(text, "splice_hash") == [t.get_splice_hash()]


def test_splice_pattern_collapse_keeps_identity():
    a, b = _model("a", ref_FSM_reads=3), _model("b", ref_FSM_reads=5)
    for t in (a, b):
        t.set_ref_splice_pattern(True)
        t.set_gene_id("g1")
    merged = SplicePatternCollapse.merge_isoforms([a, b])
    assert merged.get_ref_splice_pattern() is True and merged.get_ref_FSM_reads() == 5


# ---- guide GTF class ------------------------------------------------------------------


def test_reference_annotation_guide(tmp_path):
    (t,) = _parse(tmp_path, _guide_rows("HAVANA"))
    assert t.is_reference_annotation_model() and not t.is_LRAA_guide_model()
    assert t.is_reference_splice_pattern_source()
    assert ReferenceSupport.validate_guide_transcripts([t], "gencode") == (1, 0, 0)


@pytest.mark.parametrize(("flag", "expect"), [("True", True), ("False", False)])
def test_LRAA_guide_with_flags(tmp_path, flag, expect):
    (t,) = _parse(tmp_path, _guide_rows("LRAA", f' ref_splice_pattern "{flag}";'))
    assert t.is_LRAA_guide_model()
    assert t.is_reference_splice_pattern_source() is expect
    assert ReferenceSupport.validate_guide_transcripts([t], "lraa") == (0, 1, 0)


def test_unflagged_LRAA_guide_refused(tmp_path, monkeypatch):
    (t,) = _parse(tmp_path, _guide_rows("LRAA"))
    (se,) = _parse(tmp_path, _se_rows("LRAA"), "se.gtf")
    assert not t.is_reference_splice_pattern_source()
    for model in (t, se):
        with pytest.raises(RuntimeError, match="allow_unflagged_LRAA_guide"):
            ReferenceSupport.validate_guide_transcripts([model], "old_lraa.gtf")
    monkeypatch.setitem(LRAA_Globals.config, "allow_unflagged_LRAA_guide", True)
    assert ReferenceSupport.validate_guide_transcripts([t, se], "old_lraa.gtf") == (0, 2, 2)
    assert t.is_reference_splice_pattern_source()  # the pre-v0.50.0 behaviour
    assert se.is_reference_single_exon_source()


def test_validation_is_noop_with_rules_off(tmp_path, monkeypatch):
    monkeypatch.setitem(LRAA_Globals.config, "ref_splice_pattern_reprieve", False)
    monkeypatch.setitem(LRAA_Globals.config, "ref_supported_monoexonic", False)
    (t,) = _parse(tmp_path, _guide_rows("LRAA"))
    assert ReferenceSupport.validate_guide_transcripts([t], "x") is None


# ---- the reprieve rule ----------------------------------------------------------------


def test_reprieve_rule():
    f = TranscriptFiltering.ref_splice_pattern_reprieved
    assert f(_model(ref_FSM_reads=1), 1)
    assert not f(_model(ref_FSM_reads=5), 2)  # not the last standing carrier
    assert not f(_model(ref_FSM_reads=0), 1)  # no read carries the chain
    assert not f(_model(ref_FSM_reads=None), 1)  # not counted
    assert not f(_model(novel=True, ref_FSM_reads=5), 1)  # novel chain


def test_reprieve_threshold_and_filter_mode(monkeypatch):
    f = TranscriptFiltering.ref_splice_pattern_reprieved
    monkeypatch.setitem(LRAA_Globals.config, "ref_splice_pattern_reprieve_min_FSM_reads", 3)
    assert not f(_model(ref_FSM_reads=2), 1) and f(_model(ref_FSM_reads=3), 1)
    monkeypatch.setitem(LRAA_Globals.config, "ref_trans_filter_mode", "retain_filtered")
    assert not f(_model(ref_FSM_reads=3), 1)


def test_reprieve_uses_ref_supported_for_single_exon(monkeypatch):
    t = _se("x", 1000, 1999)
    t.init_quant_info()
    t.set_read_counts_assigned(5)
    t.set_source_reference_transcript_ids({"init_model"})  # contains a guide model
    t.set_ref_supported(False)
    assert not TranscriptFiltering.reference_model_reprieved(t)
    t.set_ref_supported(True)
    assert TranscriptFiltering.reference_model_reprieved(t)
    monkeypatch.setitem(LRAA_Globals.config, "ref_supported_monoexonic", False)
    t.set_ref_supported(False)
    assert TranscriptFiltering.reference_model_reprieved(t)  # containment, as v0.49.0


# ---- single-exon ref_supported --------------------------------------------------------


def _me_ref(tid, exons):
    t = Transcript("chr1", exons, "+")
    t.set_transcript_id(tid)
    t._guide_class = "reference_annotation"
    return t


def test_ref_supported_rule():
    ref = _se("refSE", 1000, 1999, guide_class="reference_annotation")
    me_ref = _me_ref("refME", [[3000, 3500], [4000, 4500]])
    models = {
        "match": _se("m1", 1050, 1999),  # reciprocal 0.95
        "short": _se("m2", 1000, 1500),  # half the reference
        "long": _se("m3", 900, 3000),  # model much longer
        "minus": _se("m4", 1000, 1999, "-"),  # opposite strand
        "ME_exon": _se("m5", 4000, 4400),  # inside an exon of a multi-exon reference
        "intron": _se("m6", 3400, 4100),  # spans a reference intron
    }
    assert ReferenceSupport.assign_ref_supported_monoexonic(list(models.values()), [ref, me_ref]) == 1
    assert {k for k, t in models.items() if t.is_ref_supported()} == {"match"}


def test_ref_supported_from_LRAA_guide_flags_and_no_guide():
    g_true = _se("g1", 1000, 1999, guide_class="LRAA_guide_model", ref_supported=True)
    g_false = _se("g2", 5000, 5999, guide_class="LRAA_guide_model", ref_supported=False)
    a, b = _se("a", 1000, 1990), _se("b", 5000, 5990)
    ReferenceSupport.assign_ref_supported_monoexonic([a, b], [g_true, g_false])
    assert a.get_ref_supported() is True and b.get_ref_supported() is False
    c = _se("c", 1000, 1990)
    ReferenceSupport.assign_ref_supported_monoexonic([c], None)  # de novo
    assert c.get_ref_supported() is False


# ---- full-bam reference chain counts --------------------------------------------------


def _bam(tmp_path, reads, contig_len=10000):
    """reads: (0-based start, cigar, flag, transcribed strand)."""
    header = {"HD": {"VN": "1.6", "SO": "coordinate"}, "SQ": [{"SN": "chr1", "LN": contig_len}]}
    unsorted, path = str(tmp_path / "u.bam"), str(tmp_path / "r.bam")
    with pysam.AlignmentFile(unsorted, "wb", header=header) as ofh:
        for i, (start, cigar, flag, strand) in enumerate(reads):
            a = pysam.AlignedSegment()
            a.query_name = f"r{i}"
            a.reference_id = 0
            a.reference_start = start
            a.cigarstring = cigar
            a.flag = flag
            a.mapping_quality = 60
            qlen = sum(int(n) for n, op in re.findall(r"(\d+)([MIS=X])", cigar))
            a.query_sequence = "A" * qlen
            a.query_qualities = pysam.qualitystring_to_array("I" * qlen)
            a.set_tag("NM", 0)
            a.set_tag("ts", strand)
            ofh.write(a)
    pysam.sort("-o", path, unsorted)
    pysam.index(path)
    return path


def _ref_chain():
    """Reference chain: exons 101-200, 301-500 (intron 201-300)."""
    ref = Transcript("chr1", [[101, 200], [301, 500]], "+")
    ref.set_transcript_id("R")
    ref._guide_class = "reference_annotation"
    return ref


def test_full_bam_counts_exact_chains_only(tmp_path):
    h = _ref_chain().get_splice_hash()
    bam = _bam(tmp_path, [
        (120, "80M100N150M", 0, "+"),  # exact chain
        (150, "50M100N100M", 0, "+"),  # exact chain, shorter exons
        (150, "50M150N100M", 0, "+"),  # another intron
        (120, "80M100N150M", 256, "+"),  # secondary
        (120, "80M100N150M", 0, "-"),  # other strand
        (300, "200M", 0, "+"),  # unspliced
        (120, "80M100D150M", 0, "+"),  # long deletion: an intron to quant, so it counts
    ])
    counts = ReferenceSupport.count_reads_by_reference_chain(bam, "chr1", "+", None, None, {h})
    assert counts == {h: 3}


def test_full_bam_counts_reads_crossing_a_work_unit_boundary(tmp_path):
    """A read carrying the chain is counted by a unit whose region it only overlaps.

    Chunk and --region work units fetch reads by overlap (Bam_alignment_extractor), so
    a read starting before the unit's left bound, or running past its right bound, still
    counts for the reference chain the unit holds; a unit the read does not reach does
    not count it. (Reads a chunked run severs at a cut are absent from both chunk bams
    and are counted by neither.)"""
    h = _ref_chain().get_splice_hash()
    bam = _bam(tmp_path, [
        (40, "160M100N150M", 0, "+"),  # starts at 41, before the unit at 150
        (120, "80M100N300M", 0, "+"),  # ends at 600, past the unit ending at 450
        (120, "80M100N150M", 0, "+"),  # inside
    ])
    count = ReferenceSupport.count_reads_by_reference_chain
    assert count(bam, "chr1", "+", 150, 450, {h}) == {h: 3}
    assert count(bam, "chr1", "+", 550, 900, {h}) == {h: 1}  # only the read reaching 600
    assert count(bam, "chr1", "+", 700, 900, {h}) == {h: 0}


def test_assign_ref_FSM_reads_scans_or_takes_pooled_counts(tmp_path):
    ref = _ref_chain()
    h = ref.get_splice_hash()
    bam = _bam(tmp_path, [(120, "80M100N150M", 0, "+"), (150, "50M100N100M", 0, "+")])
    known = Transcript("chr1", [[110, 200], [301, 480]], "+")  # same chain, other ends
    known.set_has_novel_splice_pattern(False)
    novel = Transcript("chr1", [[110, 180], [301, 480]], "+")
    novel.set_has_novel_splice_pattern(True)

    ReferenceSupport.assign_ref_FSM_reads([known, novel], [ref], {h}, bam, "chr1", "+", None, None)
    assert known.get_ref_FSM_reads() == 2 and novel.get_ref_FSM_reads() is None

    # an LRAA guide carrying pooled counts: used as they are, no scan
    pooled = Transcript("chr1", [[101, 200], [301, 500]], "+")
    pooled._guide_class = "LRAA_guide_model"
    pooled.set_ref_splice_pattern(True)
    pooled.set_ref_FSM_reads(9)
    ReferenceSupport.assign_ref_FSM_reads([known], [pooled], {h}, bam, "chr1", "+", None, None)
    assert known.get_ref_FSM_reads() == 9

    # any reference-chain guide row without a count: scan instead
    pooled.set_ref_FSM_reads(None)
    assert ReferenceSupport.pooled_ref_FSM_reads_from_guides([pooled]) is None
    se_guide = _se("x", 1, 2, guide_class="LRAA_guide_model", ref_supported=True)
    pooled.set_ref_FSM_reads(9)
    assert ReferenceSupport.pooled_ref_FSM_reads_from_guides([se_guide, pooled]) == {h: 9}


def test_oversimplified_contig_models_are_flagged():
    """Rows an oversimplified contig carries unchanged get their identity, so the GTF
    can guide a later run."""
    ref_me = _me_ref("refME", [[100, 200], [300, 500]])
    ref_se = _se("refSE", 600, 700, guide_class="reference_annotation")
    lraa_me = Transcript("chr1", [[100, 200], [300, 500]], "+")
    lraa_me._guide_class = "LRAA_guide_model"
    lraa_me.set_ref_splice_pattern(False)
    aggregate = _se("agg", 1, 900)  # built in-run (de novo aggregate)
    ReferenceSupport.record_carried_model_identity([ref_me, ref_se, lraa_me, aggregate])
    assert ref_me.get_ref_splice_pattern() is True and ref_se.get_ref_supported() is True
    assert lraa_me.get_ref_splice_pattern() is False
    assert aggregate.get_ref_supported() is False


def test_guide_gtf_file_checked_at_startup(tmp_path, monkeypatch):
    """The streamed file check LRAA runs once before any work unit: same verdicts as the
    parsed check, oversimplified contigs exempt, gzip read."""
    import gzip

    def write(name, lines, gz=False):
        path = tmp_path / name
        text = "\n".join(lines) + "\n"
        if gz:
            with gzip.open(path, "wt") as fh:
                fh.write(text)
        else:
            path.write_text(text)
        return str(path)

    check = ReferenceSupport.validate_guide_gtf_file
    assert check(write("gencode.gtf", _guide_rows("HAVANA") + _se_rows("HAVANA"))) == (2, 0, 0)
    flagged = _guide_rows("LRAA", ' ref_splice_pattern "False";') + _se_rows("LRAA", ' ref_supported "True";')
    assert check(write("flagged.gtf.gz", flagged, gz=True)) == (0, 2, 0)
    unflagged = write("old.gtf", _guide_rows("LRAA") + _se_rows("LRAA"))
    with pytest.raises(RuntimeError, match="allow_unflagged_LRAA_guide"):
        check(unflagged)
    assert check(unflagged, exempt_contigs={"chr1"}) == (0, 0, 0)
    monkeypatch.setitem(LRAA_Globals.config, "allow_unflagged_LRAA_guide", True)
    assert check(unflagged) == (0, 2, 2)
