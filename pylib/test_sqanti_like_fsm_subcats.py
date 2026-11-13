#!/usr/bin/env python3

"""FSM subcategories of the SQANTI-like classifier.

FSM stays the top-level category; extra fields record
  FSM_full_length_isoforms: matching ref isoforms whose both ends are within
      fsm_full_length_max_end_dist bases of the query's ends;
  FSM_ISM_of_isoforms: other ref isoforms whose intron chain properly contains the
      query's chain (the query would be their ISM), with their genes
      (FSM_ISM_of_genes) and whether those are the FSM's gene (FSM_ISM_of_gene_scope).
"""

import os
import pytest
from SQANTI_like_annotator import SQANTI_like_annotator
from MockTranscript import MockTranscript
import LRAA_Globals


def _exon_lines(tid, gene, exons, strand="+"):
    return "".join(
        'chr1\tsrc\texon\t{}\t{}\t.\t{}\t.\tgene_id "{}"; transcript_id "{}";\n'.format(
            lend, rend, strand, gene, tid
        )
        for lend, rend in exons
    )


# exons of the reference isoforms; introns 201-299, 401-499, 601-699
REF = {
    # A: two introns (201-299, 401-499)
    "A": [(100, 200), (300, 400), (500, 1000)],
    # A2: same intron chain as A, different ends
    "A2": [(150, 200), (300, 400), (500, 900)],
    # B: A's chain plus a downstream intron; A's chain is a proper subset of B's
    "B": [(100, 200), (300, 400), (500, 600), (700, 800)],
    # C: shares only A's first intron, skips exon 2: not a container of A's chain
    "C": [(100, 200), (500, 600)],
    # S: single-exon reference isoform
    "S": [(2000, 3000)],
}


@pytest.fixture
def annotator(tmp_path):
    gtf = os.path.join(tmp_path, "ref.gtf")
    with open(gtf, "w") as fh:
        for tid, exons in REF.items():
            fh.write(_exon_lines(tid, "g1", exons))
    return SQANTI_like_annotator(gtf, fsm_full_length_max_end_dist=50)


def classify(annotator, exons, name="q"):
    q = MockTranscript(name, "chr1", "+", exons)
    return annotator.classify_alignment_or_isoform("chr1", "+", name, q)


def test_fsm_of_A_is_also_ism_of_B(annotator):
    r = classify(annotator, [(100, 200), (300, 400), (500, 1000)])
    assert r["sqanti_cat"] == "FSM"
    assert r["matching_isoforms"] == "A,A2"
    # B contains the chain; C (shares one intron only) and the FSM matches do not count
    assert r["FSM_ISM_of_isoforms"] == "B"


def test_full_length_within_tolerance_of_the_same_ref(annotator):
    r = classify(annotator, [(150, 200), (300, 400), (500, 950)])
    assert r["sqanti_cat"] == "FSM"
    # A: |150-100|=50, |950-1000|=50 ; A2: |150-150|=0, |950-900|=50
    assert r["FSM_full_length_isoforms"] == "A,A2"


def test_ends_just_outside_tolerance(annotator):
    r = classify(annotator, [(151, 200), (300, 400), (500, 1051)])
    assert r["sqanti_cat"] == "FSM"
    # A: both ends off by 51 ; A2: lend off by 1 but rend off by 151
    assert r["FSM_full_length_isoforms"] == ""
    # one end within tolerance is not enough
    r = classify(annotator, [(100, 200), (300, 400), (500, 1051)])
    assert r["FSM_full_length_isoforms"] == ""


def test_multiple_ref_matches_only_those_within_tolerance(annotator):
    # exactly A's ends: full-length vs A (0/0); vs A2 off by 50 / 100 -> not
    r = classify(annotator, [(100, 200), (300, 400), (500, 1000)])
    assert r["matching_isoforms"] == "A,A2"
    assert r["FSM_full_length_isoforms"] == "A"


def test_full_length_and_ism_are_independent(annotator):
    # B itself: FSM full-length of B, and not contained in any other ref
    r = classify(annotator, [(100, 200), (300, 400), (500, 600), (700, 800)])
    assert r["sqanti_cat"] == "FSM"
    assert r["matching_isoforms"] == "B"
    assert r["FSM_full_length_isoforms"] == "B"
    assert r["FSM_ISM_of_isoforms"] == ""
    # a short FSM of A/A2 far from all ends: ISM of B but not full-length
    r = classify(annotator, [(190, 200), (300, 400), (500, 510)])
    assert r["FSM_full_length_isoforms"] == ""
    assert r["FSM_ISM_of_isoforms"] == "B"


def test_tolerance_is_a_parameter(tmp_path):
    gtf = os.path.join(tmp_path, "ref.gtf")
    with open(gtf, "w") as fh:
        fh.write(_exon_lines("A", "g1", REF["A"]))
    q = [(110, 200), (300, 400), (500, 990)]  # 10 b off at each end
    strict = SQANTI_like_annotator(gtf, fsm_full_length_max_end_dist=9)
    loose = SQANTI_like_annotator(gtf, fsm_full_length_max_end_dist=10)
    assert classify(strict, q)["FSM_full_length_isoforms"] == ""
    assert classify(loose, q)["FSM_full_length_isoforms"] == "A"


def test_default_tolerance_comes_from_lraa_globals(tmp_path):
    gtf = os.path.join(tmp_path, "ref.gtf")
    with open(gtf, "w") as fh:
        fh.write(_exon_lines("A", "g1", REF["A"]))
    assert LRAA_Globals.config["sqlike_fsm_full_length_max_end_dist"] == 50
    assert SQANTI_like_annotator(gtf).fsm_full_length_max_end_dist == 50


def test_subcats_empty_for_non_fsm(annotator):
    # ISM of A/A2/B: top-level ISM, FSM-only fields stay empty
    r = classify(annotator, [(150, 200), (300, 350)])
    assert r["sqanti_cat"] == "ISM"
    assert r["FSM_full_length_isoforms"] == ""
    assert r["FSM_ISM_of_isoforms"] == ""
    assert r["FSM_ISM_of_genes"] == r["FSM_ISM_of_gene_scope"] == ""


def test_monoexonic_se_fm_has_no_subcats(annotator):
    # mono-exonic match to S: se_FM; the subcategories are defined for spliced FSMs
    # only (no intron chain to contain), so both fields stay empty even at S's ends
    r = classify(annotator, [(2000, 3000)])
    assert r["sqanti_cat"] == "se_FM"
    assert r["matching_isoforms"] == "S"
    assert r["FSM_full_length_isoforms"] == ""
    assert r["FSM_ISM_of_isoforms"] == ""
    assert r["FSM_ISM_of_genes"] == r["FSM_ISM_of_gene_scope"] == ""


def _read(cigar, reference_start, is_reverse=False):
    import pysam

    header = pysam.AlignmentHeader.from_dict(
        {"HD": {"VN": "1.6"}, "SQ": [{"SN": "chr1", "LN": 100000}]}
    )
    a = pysam.AlignedSegment(header)
    a.query_name = "read"
    a.reference_id = 0
    a.reference_start = reference_start  # 0-based
    a.mapping_quality = 60
    a.flag = 16 if is_reverse else 0
    a.cigarstring = cigar
    qlen = sum(n for op, n in a.cigartuples if op in (0, 1, 4, 7, 8))
    a.query_sequence = "C" * qlen
    a.query_qualities = pysam.qualitystring_to_array("I" * qlen)
    return a


def test_read_mode_uses_the_alignment_span(annotator):
    # read mode: a read aligned to A's exons 130-200, 300-400, 500-960, with 40 b of
    # 5' soft clip that does not count towards its span. 30 b / 40 b from A's ends.
    from Pretty_alignment import Pretty_alignment

    read = _read("40S71M99N101M99N461M", reference_start=129)
    pa = Pretty_alignment.get_pretty_alignment(read)
    assert tuple(pa.get_alignment_span()) == (130, 960)
    r = annotator.classify_alignment_or_isoform("chr1", "+", "read", pa)
    assert r["sqanti_cat"] == "FSM"
    assert r["matching_isoforms"] == "A,A2"
    # A2 (150-900): 3' end off by 60
    assert r["FSM_full_length_isoforms"] == "A"
    assert r["FSM_ISM_of_isoforms"] == "B"


def _gene_gtf(tmp_path, isoforms):
    # isoforms: {tid: (gene_id, gene_name or None, exons)}
    gtf = os.path.join(tmp_path, "genes.gtf")
    with open(gtf, "w") as fh:
        for tid, (gene, name, exons) in isoforms.items():
            name_attr = ' gene_name "{}";'.format(name) if name else ""
            for lend, rend in exons:
                fh.write(
                    'chr1\tsrc\texon\t{}\t{}\t.\t+\t.\tgene_id "{}";{} transcript_id "{}";\n'.format(
                        lend, rend, gene, name_attr, tid
                    )
                )
    return SQANTI_like_annotator(gtf, fsm_full_length_max_end_dist=50)


A_EXONS = REF["A"]
B_EXONS = REF["B"]
# a readthrough-like isoform of another gene that also carries A's chain
R_EXONS = [(100, 200), (300, 400), (500, 600), (700, 800), (5000, 5100)]


def test_containing_gene_reported_same_gene(tmp_path):
    ann = _gene_gtf(
        tmp_path,
        {"A": ("g1", "G1", A_EXONS), "B": ("g1", "G1", B_EXONS)},
    )
    r = classify(ann, A_EXONS)
    assert r["FSM_ISM_of_isoforms"] == "B"
    assert r["FSM_ISM_of_genes"] == "G1^g1"
    assert r["FSM_ISM_of_gene_scope"] == "same_gene"


def test_containing_gene_reported_other_gene_only(tmp_path):
    # structural containment counts across genes; the scope column says so
    ann = _gene_gtf(
        tmp_path,
        {"A": ("g1", "G1", A_EXONS), "R": ("g2", None, R_EXONS)},
    )
    r = classify(ann, A_EXONS)
    assert r["FSM_ISM_of_isoforms"] == "R"
    assert r["FSM_ISM_of_genes"] == "g2"  # no gene_name: gene_id alone
    assert r["FSM_ISM_of_gene_scope"] == "other_gene"


def test_containing_genes_parallel_to_isoforms_mixed_scope(tmp_path):
    ann = _gene_gtf(
        tmp_path,
        {
            "A": ("g1", "G1", A_EXONS),
            "B": ("g1", "G1", B_EXONS),
            "R": ("g2", "G2", R_EXONS),
        },
    )
    r = classify(ann, A_EXONS)
    assert r["FSM_ISM_of_isoforms"] == "B,R"
    assert r["FSM_ISM_of_genes"] == "G1^g1,G2^g2"
    assert r["FSM_ISM_of_gene_scope"] == "same_and_other_gene"


def test_gene_fields_empty_without_containment(annotator):
    r = classify(annotator, [(100, 200), (300, 400), (500, 600), (700, 800)])
    assert r["sqanti_cat"] == "FSM"
    assert r["FSM_ISM_of_genes"] == ""
    assert r["FSM_ISM_of_gene_scope"] == ""


def test_end_offsets_signed_from_closest_ref(annotator):
    # vs A (100-1000): 5' +10 past, 3' -20 short ; vs A2 (150-900): -40 / +80.
    # A is closer (larger offset 20 vs 80)
    r = classify(annotator, [(90, 200), (300, 400), (500, 980)])
    assert r["matching_isoforms"] == "A,A2"
    assert (r["FSM_ref_5prime_offset"], r["FSM_ref_3prime_offset"]) == (10, -20)
    # full-length exactly when both offsets are within the tolerance
    assert r["FSM_full_length_isoforms"] == "A"


def test_end_offsets_minus_strand(tmp_path):
    # on the minus strand the 5' end is the right end
    gtf = os.path.join(tmp_path, "ref.gtf")
    with open(gtf, "w") as fh:
        fh.write(_exon_lines("M", "g1", REF["A"], strand="-"))
    ann = SQANTI_like_annotator(gtf, fsm_full_length_max_end_dist=50)
    q = MockTranscript("q", "chr1", "-", [(130, 200), (300, 400), (500, 1100)])
    r = ann.classify_alignment_or_isoform("chr1", "-", "q", q)
    assert r["sqanti_cat"] == "FSM"
    # 5' (right) end 100 past, 3' (left) end 30 short
    assert (r["FSM_ref_5prime_offset"], r["FSM_ref_3prime_offset"]) == (100, -30)


def test_end_offsets_empty_for_non_fsm(annotator):
    r = classify(annotator, [(150, 200), (300, 350)])
    assert r["sqanti_cat"] == "ISM"
    assert r["FSM_ref_5prime_offset"] == r["FSM_ref_3prime_offset"] == ""


def test_fsm_subcat_tallies():
    import importlib.util
    from collections import defaultdict

    path = os.path.join(
        os.path.dirname(__file__),
        "..",
        "util",
        "SQANTI-like_cats_for_reads_or_isoforms.py",
    )
    spec = importlib.util.spec_from_file_location("sqlike_script", path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    counter = defaultdict(int)
    base = dict(FSM_full_length_isoforms="", FSM_ISM_of_gene_scope="")
    mod.count_fsm_subcats(dict(base, sqanti_cat="ISM"), counter)
    assert not counter
    mod.count_fsm_subcats(
        dict(
            base,
            sqanti_cat="FSM",
            FSM_full_length_isoforms="A",
            FSM_ref_5prime_offset=-5,
            FSM_ref_3prime_offset=4999,
        ),
        counter,
    )
    mod.count_fsm_subcats(
        dict(
            base,
            sqanti_cat="FSM",
            FSM_ISM_of_gene_scope="other_gene",
            FSM_ref_5prime_offset=-1,
            FSM_ref_3prime_offset=-1001,
        ),
        counter,
    )
    assert counter[("subcat", "yes", "none")] == 1
    assert counter[("subcat", "no", "other_gene")] == 1
    # 10-nt bins by lower bound, clipped at +/-1000
    assert counter[("offset", "5prime", -10)] == 2
    assert counter[("offset", "3prime", 1000)] == 1
    assert counter[("offset", "3prime", -1000)] == 1
