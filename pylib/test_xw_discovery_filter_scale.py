#!/usr/bin/env python3

"""Discovery's isoform filters must measure support on ONE scale, for any bam.

Some of those filters are EM-derived and follow whatever quantity EM apportioned
(min_reads_novel_isoform, the TPM and isoform-fraction gates); the rest are integer
tallies over retained reads and cannot follow it (min_unique_reads_novel_isoform, the FSM
counts, the supporting-cell gate, the trellis score). If the pass those filters read were
weighted, the two families would disagree by the thinning factor and isoforms would be
selected under one notion of support and filtered under another.

XW weighting is no longer a mode. LRAA_Globals no longer carries a
use_XW_read_weights_for_quant key, both CLI flags are gone, and
LRAA._populate_read_multi_paths computes `use_XW_weights = True if weight_reads is None
else bool(weight_reads)` -- weights apply to every pass unless a caller says otherwise.
Discovery's pre-filter quantification is the one caller that says otherwise: it pins
weight_reads=False (LRAA, "Do an initial quant"). That pin is the whole guarantee, and
these tests are what holds it in place.

There is no longer a CLI route to this contract. A thinned, XW-tagged bam handed to
--bam is now refused outright by LRAA._require_no_thinning_weights, and there is no flag
left to toggle weighting on the surviving inputs, so a two-run CLI comparison cannot be
written any more. The pin is instead exercised where it lives: build a splice graph over
a thinned XW-tagged bam, then run LRAA.build_multipath_graph across it twice, once with
weight_reads=False and once with the default, and read the multipaths' support directly.

The fixture is one locus with two isoforms and deliberately UNEVEN thinning, because an
evenly thinned locus would hide the defect: EM's relative gates are scale-free, so a
uniform factor cancels out of them. The deep isoform sits far above the depth target and
is thinned to roughly a quarter, so its reads come back carrying XW ~4; the scarce
isoform sits below the target and is kept whole at XW 1. Weighted, the deep isoform's
surviving reads each stand for four and inflate the gene total the scarce isoform's share
is measured against, so the relative gates cut the scarce one -- while its own
unique-read and FSM tallies, being integer counts of retained reads, never moved.

That is not hypothetical. Deleting `weight_reads=False` from the call site and running
real discovery over this fixture dropped exactly the scarce model,
((1101,1200),(8001,8100)) -- the two-scales failure in one line of output.
"""

import os
import sys

import pysam
import pytest

REPO_ROOT = os.path.dirname(os.path.dirname(os.path.realpath(__file__)))

for _sub in ("util", "pylib"):
    _p = os.path.join(REPO_ROOT, _sub)
    if _p not in sys.path:
        sys.path.insert(0, _p)

import normalize_bam_by_strand as norm  # noqa: E402
import LRAA_Globals  # noqa: E402
# pylib/LRAA.py, holding the assembler class -- not the top-level CLI script of the same
# name, which is not importable and is no longer what this file drives.
import LRAA as LRAA_module  # noqa: E402
from Splice_graph import Splice_graph  # noqa: E402

CONTIG = "chr1"
CONTIG_LEN = 20000
# Deep enough that normalization thins it hard: the target is 1000, so a locus at several
# thousand reads gets acceptance probabilities well below 1 and XW values well above it.
N_DEEP = 4000
# Above the alt-splice frequency floor so the scarce junction survives into the graph
# (intron support is itself XW-weighted, so this is measured against ~4000, not against
# the deep isoform's surviving read count), yet still under the depth target so its own
# reads are kept whole. That gap between the two isoforms is the point of the fixture.
N_SCARCE = 200

# Multipath spans, 1-based inclusive, as get_coords() reports them.
DEEP = (1101, 3100)
SCARCE = (1101, 8100)


def _header():
    return pysam.AlignmentHeader.from_dict(
        {"HD": {"VN": "1.6", "SO": "coordinate"},
         "SQ": [{"SN": CONTIG, "LN": CONTIG_LEN}]}
    )


def _mk_read(hdr, name, ref_start, cigar):
    a = pysam.AlignedSegment(hdr)
    a.query_name = name
    a.reference_id = 0
    a.reference_start = ref_start
    a.mapping_quality = 60
    a.cigarstring = cigar
    span = sum(l for op, l in a.cigartuples if op in (0, 1))
    a.query_sequence = "A" * span
    a.query_qualities = pysam.qualitystring_to_array("I" * span)
    return a


def _build_bam(path):
    hdr = _header()
    reads = []
    # Deep isoform: exon1 -> exon2. Its junction is far above the target and gets thinned.
    for i in range(N_DEEP):
        reads.append(_mk_read(hdr, f"deep_{i}", 1100, "100M1800N100M"))
    # Scarce isoform: exon1 -> exon3. Below the target, so its junction is kept whole --
    # this is the mixture of acceptance rates within one locus that makes a weighted pass
    # and an unweighted one disagree at all.
    for i in range(N_SCARCE):
        reads.append(_mk_read(hdr, f"scarce_{i}", 1100, "100M6800N100M"))
    with pysam.AlignmentFile(str(path), "wb", header=hdr) as fh:
        for r in sorted(reads, key=lambda x: x.reference_start):
            fh.write(r)
    pysam.index(str(path))


def _genome_seq():
    """All A except canonical GT..AG at the two introns the reads span.

    The splice graph derives introns from the alignments and requires the canonical
    dinucleotides, so a uniform genome yields no junctions, no multipaths, and every
    comparison below would be a statement about an empty collection. Both reads share a
    donor at 0-based 1200; the acceptors are at 2998 (deep) and 7998 (scarce).
    """
    seq = bytearray(b"A" * CONTIG_LEN)
    seq[1200:1202] = b"GT"
    seq[2998:3000] = b"AG"
    seq[7998:8000] = b"AG"
    return seq.decode()


def _thin(true_bam, thinned_bam):
    """Thinned with the same tool and defaults LRAA's own --normalize_max_cov_level uses."""
    norm.sift_bam(
        str(true_bam),
        str(thinned_bam),
        normalize_max_cov_level=LRAA_Globals.config["normalize_max_cov_level"],
        depth_window=LRAA_Globals.config["chunk_depth_window"],
        random_seed=LRAA_Globals.config["chunk_random_seed"],
        window_origin=LRAA_Globals.config["chunk_grid_origin"],
        min_per_id=LRAA_Globals.config["min_per_id"],
        min_mapping_quality=LRAA_Globals.config["min_mapping_quality_for_final_quant"],
    )
    pysam.index(str(thinned_bam))


def _multipath_support(bam, contig_seq, weight_reads):
    """span -> (read count, read weight) for one build_multipath_graph pass.

    Keyed on the multipath's genomic span rather than on node ids, because the splice
    graph is rebuilt per pass and Exon ids are handed out by a running counter: the same
    structure is 'E:1,I:1,E:2' in the first pass and 'E:4,I:3,E:5' in the second.
    """
    splice_graph = Splice_graph()
    splice_graph.build_splice_graph_for_contig(
        CONTIG, "+", contig_seq, str(bam), None, None,
        input_transcripts=None, quant_mode=False,
    )
    assembler = LRAA_module.LRAA(splice_graph)
    mp_counter = assembler.build_multipath_graph(
        CONTIG, "+", contig_seq, str(bam), weight_reads=weight_reads
    )
    support = {}
    for pair in mp_counter.get_all_MultiPathCountPairs():
        multipath, _count = pair.get_multipath_and_count()
        support[multipath.get_coords()] = (
            multipath.get_read_count(),
            multipath.get_read_weight(),
        )
    return support


@pytest.fixture(scope="module")
def multipath_support(tmp_path_factory):
    tmp = tmp_path_factory.mktemp("xw_discovery_scale")
    true_bam = tmp / "true.bam"
    _build_bam(true_bam)
    thinned = tmp / "thinned.bam"
    _thin(true_bam, thinned)
    contig_seq = _genome_seq()

    # build_multipath_graph writes read-name/read-id stores and a pretty-alignment cache
    # relative to LRAA_TMP_DIR or the working directory, so both are pointed at the
    # tmpdir; otherwise a pytest run scatters __read_tracking under the repo root. The
    # stores are process-global and are dropped so they get created here rather than
    # inherited from whatever test module ran first.
    prior_cwd = os.getcwd()
    prior_tmp_dir = os.environ.get("LRAA_TMP_DIR")
    prior_stores = (LRAA_Globals.READ_NAME_STORE, LRAA_Globals.MP_READ_ID_STORE)
    os.chdir(tmp)
    os.environ["LRAA_TMP_DIR"] = str(tmp)
    LRAA_Globals.READ_NAME_STORE = None
    LRAA_Globals.MP_READ_ID_STORE = None
    try:
        results = {"thinned_bam": thinned}
        # weight_reads=False is discovery's pre-filter pass; None is what every other
        # caller gets, and is the default the pin exists to override.
        results["pinned"] = _multipath_support(thinned, contig_seq, False)
        results["default"] = _multipath_support(thinned, contig_seq, None)
    finally:
        LRAA_Globals.READ_NAME_STORE, LRAA_Globals.MP_READ_ID_STORE = prior_stores
        LRAA_Globals.reset_read_weight_registry()
        os.chdir(prior_cwd)
        if prior_tmp_dir is None:
            os.environ.pop("LRAA_TMP_DIR", None)
        else:
            os.environ["LRAA_TMP_DIR"] = prior_tmp_dir
    return results


def test_the_fixture_bam_really_carries_weights_above_one(multipath_support):
    """Vacuity guard. Without XW values above 1 in this bam, every other test here would
    pass by measuring nothing: an untagged read weighs 1, so weight and count would agree
    under both passes for reasons that have nothing to do with the pin.
    """
    weights = []
    with pysam.AlignmentFile(str(multipath_support["thinned_bam"])) as fh:
        for aln in fh:
            if aln.has_tag("XW"):
                weights.append(float(aln.get_tag("XW")))
    assert weights, "the thinned bam must carry XW tags at all"
    assert max(weights) > 1.5, (
        f"thinning must actually have thinned something; max XW was {max(weights)}"
    )


def test_both_isoforms_reach_the_multipath_counter(multipath_support):
    """Guards every comparison below against being a statement about an empty dict.

    Also pins the fixture's shape: two multipaths over one locus, present under both
    passes, so the pinned and default runs are compared over the same structures. The
    scarce junction is the one at risk -- alt-splice frequency filtering is measured on
    XW-weighted intron support, so dropping N_SCARCE much lower silently deletes it from
    the graph and leaves this file testing a single isoform.
    """
    for label in ("pinned", "default"):
        assert set(multipath_support[label]) == {DEEP, SCARCE}, (
            f"{label}: expected the deep and scarce multipaths, got "
            f"{sorted(multipath_support[label])}"
        )
        for span, (count, _weight) in multipath_support[label].items():
            assert count > 0, f"{label}: multipath at {span} carries no reads"


def test_the_pinned_pre_filter_pass_ignores_the_weights_it_was_handed(multipath_support):
    """The pin itself: weight_reads=False must make support a plain read tally.

    Every multipath must report get_read_weight() == get_read_count(), on a bam whose
    reads demonstrably carry XW above 1. That equality is what puts the EM-derived gates
    and the integer tallies on one scale. Delete `weight_reads=False` from the call site
    in LRAA and this is the assertion that fails first.
    """
    for span, (count, weight) in sorted(multipath_support["pinned"].items()):
        assert weight == pytest.approx(count), (
            f"the pinned pass followed the XW tags at {span}: "
            f"weight {weight} against count {count}"
        )


def test_the_default_pass_does_follow_the_weights(multipath_support):
    """The converse, so the test above proves the parameter does something.

    Without this, a pinned pass reporting weight == count would be equally consistent
    with weighting having been removed from LRAA altogether, or with the fixture bam
    having lost its tags somewhere between thinning and intake.

    Measured on this fixture: the deep multipath comes back with count 995 and weight
    4079.5, recovering the ~4000 reads it stands for.
    """
    count, weight = multipath_support["default"][DEEP]
    assert weight > count * 1.5, (
        f"the default pass did not weight the thinned isoform: "
        f"weight {weight} against count {count}"
    )


def test_the_default_pass_reweights_the_locus_unevenly(multipath_support):
    """Why a uniform thinning factor would not have caught this.

    EM's relative gates are scale-free, so if every read in the locus carried the same
    weight the weighted and unweighted passes would agree and the pin would be
    unnecessary. The damage comes from the MIXTURE: here the deep isoform is inflated
    ~4x while the scarce one, kept whole at XW 1, is not touched -- so weighting moves
    the scarce isoform's share of the gene total without moving any of the integer
    tallies applied beside it. Measured: scarce count 200 and weight 200.0 against the
    deep isoform's 995 and 4079.5.
    """
    scarce_count, scarce_weight = multipath_support["default"][SCARCE]
    deep_count, deep_weight = multipath_support["default"][DEEP]
    assert scarce_weight == pytest.approx(scarce_count), (
        f"the scarce isoform was supposed to be kept whole at XW 1; "
        f"weight {scarce_weight} against count {scarce_count}"
    )
    assert (deep_weight / deep_count) > 2 * (scarce_weight / scarce_count), (
        "the two isoforms must be reweighted by visibly different factors, or this "
        "fixture cannot distinguish a weighted pass from an unweighted one"
    )


@pytest.fixture(scope="module")
def gate_inputs(tmp_path_factory):
    """A deep locus, its thinned counterpart, and a genome, built once.

    Returns the full and retained record counts alongside the paths, because the claim
    below is precisely that the gate sees the former and not the latter.
    """
    tmp = tmp_path_factory.mktemp("gate_scale")
    true_bam = tmp / "true.bam"
    _build_bam(true_bam)

    genome = tmp / "g.fa"
    seq = _genome_seq()
    with open(genome, "w") as fh:
        fh.write(f">{CONTIG}\n")
        for i in range(0, len(seq), 60):
            fh.write(seq[i : i + 60] + "\n")

    thinned = tmp / "thinned.bam"
    _thin(true_bam, thinned)

    def _count(path):
        with pysam.AlignmentFile(str(path)) as fh:
            return sum(1 for _ in fh)

    full, retained = _count(true_bam), _count(thinned)
    assert retained < 0.5 * full, (
        f"the fixture must actually be thinned for these tests to discriminate: "
        f"{retained} of {full}"
    )
    return {
        "dir": tmp,
        "bam": true_bam,
        "thinned": thinned,
        "genome": genome,
        "full": full,
        "retained": retained,
    }


# The configurations that change WHICH bam the draft quant reads. The draft mirrors the
# final quant's estimating pass by the same expression run_quant_only uses --
# bam_file_for_sg under streaming, bam_file_for_quant without it -- so streaming and
# normalization are the settings that matter, plus the pre-normalized --bam_for_sg
# composition every chunk worker runs (ChunkedRun hands workers a full mini-bam as --bam
# and the normalized one as --bam_for_sg with --no_norm). --no_norm alone is not a point
# in the space: streaming refuses to run when nothing thins its first pass.
#
# `thinned_draft` records whether the draft reads a THINNED bam in that configuration,
# which is what decides whether weighting has to do the work of recovering library scale.
GATE_MODALITIES = [
    pytest.param([], True, id="stock_normalized_and_streaming"),
    pytest.param(["--no_norm", "--no_stream_reads"], False, id="unnormalized_single_pass"),
    pytest.param(
        ["--bam_for_sg", "THINNED", "--no_norm"], True, id="chunk_worker_composition"
    ),
]


@pytest.mark.parametrize("extra_flags,thinned_draft", GATE_MODALITIES)
def test_the_filters_see_library_scale_abundance_and_literal_unique_reads(
    gate_inputs, extra_flags, thinned_draft
):
    """Both halves of the draft quant's contract, in each modality, from one real run.

    The draft quant now reads the same bam the final quant's estimating pass will read and
    weights it, so the two families of gate get deliberately different quantities:

      min_reads_novel_isoform   an ABUNDANCE estimate, on the scale of the original
                                library. Under streaming the draft reads the THINNED bam,
                                so this number is only library-scale because XW weights
                                divide the acceptance rates back out -- the retained read
                                count would be several times smaller.
      min_unique_reads_novel_isoform
                                a count of OBSERVATIONS: literal reads of that same bam.
                                Not weighted, and so not library-scale.

    Asserting both together is the point. Either alone is satisfiable by the wrong design:
    library-scale abundance alone was equally true when the draft read the full bam
    unweighted (which is what this test used to pin), and a literal unique count alone
    would be equally true if nothing were weighted anywhere.

    Measured on this fixture at 3002 full / 1033 retained: abundance 3095.0 (+3.1%, which
    is the 1/sqrt(1033) sampling error inverse-probability weighting carries), unique count
    1033 exactly.

    --no_parallelize_contigs is required rather than incidental: per-contig workers
    configure their own logging, and these lines never reach this process otherwise.
    """
    import re
    import subprocess

    flags = [
        str(gate_inputs["thinned"]) if f == "THINNED" else f for f in extra_flags
    ]
    # Its OWN directory per modality, and not incidental: --debug writes fixed-name
    # artifacts (__MPGN_components_described.<contig>.<strand>.<n>.bed) and
    # MultiPathGraph refuses to overwrite one, so two runs sharing a cwd make the second
    # fail on the first's leftovers.
    rundir = gate_inputs["dir"] / ("run_" + str(abs(hash(tuple(flags))))[:8])
    rundir.mkdir(exist_ok=True)
    out_prefix = rundir / "out"
    # Split the two thresholds, which --min_reads_novel would otherwise tie together:
    # min_reads must PASS so the run reaches the unique gate, and the unique gate must
    # FAIL so it logs the count it compared.
    cfg = rundir / "cfg.json"
    cfg.write_text(
        '{"min_reads_novel_isoform": 2, "min_unique_reads_novel_isoform": 99999999, '
        '"min_frac_gene_unique_reads": 0.0}'
    )
    cmd = [sys.executable, os.path.join(REPO_ROOT, "LRAA"),
           "--bam", str(gate_inputs["bam"]),
           "--genome", str(gate_inputs["genome"]),
           "--output_prefix", str(out_prefix),
           "--no_chunk", "--no_parallelize_contigs", "--debug",
           "--config_update", str(cfg),
           "--no_rescue_unassigned_reads_via_transcriptome_alignment"] + flags
    r = subprocess.run(cmd, capture_output=True, text=True, cwd=str(rundir))
    assert r.returncode == 0, r.stdout + r.stderr
    out = r.stdout + r.stderr

    abundances = [
        float(m)
        for m in re.findall(r"is a novel isoform with ([0-9.]+) read support", out)
    ]
    assert abundances, "the min_reads gate must have examined at least one novel isoform"
    total_abundance = sum(abundances)

    unique_counts = [
        int(m) for m in re.findall(r"too few unique reads: ([0-9]+)", out)
    ]
    assert unique_counts, "the unique gate must have logged the count it compared"
    total_unique = sum(unique_counts)

    # Half one: abundance is on the library's scale, whichever bam it came from.
    assert total_abundance == pytest.approx(gate_inputs["full"], rel=0.05), (
        f"min_reads was fed {total_abundance}, not the library's "
        f"{gate_inputs['full']}"
    )

    # Half two: the unique tally counts reads of the bam the draft actually read, and is
    # NOT rescaled to the library.
    expected_unique = (
        gate_inputs["retained"] if thinned_draft else gate_inputs["full"]
    )
    assert total_unique == expected_unique, (
        f"the unique gate compared {total_unique}; the reads available to the draft "
        f"in this configuration number {expected_unique}"
    )

    if thinned_draft:
        # And the two really are different quantities here, so neither assertion above is
        # accidentally restating the other.
        assert total_abundance > 2 * total_unique, (
            f"weighting did no work: abundance {total_abundance} is not meaningfully "
            f"above the {total_unique} reads it was estimated from"
        )