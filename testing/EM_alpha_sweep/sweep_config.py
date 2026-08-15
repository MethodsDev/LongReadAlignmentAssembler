"""Corpus and arm definitions for the EM_alpha x 3'-weighting factorial sweep.

Measurement scaffolding only -- nothing here is imported by LRAA itself.

Corpus is the 10 truth-bearing quant-only samples named in the sweep charter:
SIRV E1/E2 across four cell lines (SIRV E0 is deliberately excluded -- its truth
is equimolar, so MARD there rewards pushing mass toward uniform and recommends
the opposite of every truth-bearing sample), plus the two MORF minigenome
samples.  Arabidopsis and mouse are out of scope for this phase.
"""

import os

# ---------------------------------------------------------------- paths

LRAA_HOME = "/home/unix/bhaas/projects/SingleCellOverhaul/LRAA-quant-alpha"
LRAA = os.path.join(LRAA_HOME, "LRAA")

INPUTS = "/home/unix/bhaas/projects/LRAA_PAPER_Analyses/__LRAA_local_runs/v0.18.3_batch/inputs"
SIRV_REF = "/home/unix/bhaas/projects/LRAA_PAPER_Analyses/Houlins_compiled_data/LRAA_manuscript_benchmarking/references/2.sirvs_pb"
MORF_REF = "/home/unix/bhaas/projects/LRAA_PAPER_Analyses/MISC/morfs"
QUANT_ONLY = "/home/unix/bhaas/projects/LRAA_PAPER_Analyses/QUANT_ONLY_outputs"

BASE = "/home/unix/bhaas/projects/LRAA_PAPER_Analyses/__LRAA_local_runs/alpha_grid_2026-08-15"
WORK = os.path.join(BASE, "work")
SCORE = os.path.join(BASE, "score")

# ---------------------------------------------------------------- corpus

SIRV_GENOME = os.path.join(SIRV_REF, "SIRV_isoforms_multi-fasta_170612a.fasta")
SIRV_GTF = os.path.join(SIRV_REF, "SIRV_isoforms_multi-fasta-annotation.expressed.gtf")
SIRV_TRUTH_GTF = os.path.join(
    QUANT_ONLY, "SIRVs/reference_data/SIRV_isoforms_multi-fasta-annotation_C_170612a.gtf"
)

MORF_GENOME = os.path.join(MORF_REF, "minigenome.fa")
MORF_GTF = os.path.join(MORF_REF, "minigenome.UTRs_trimmed.wGeneName.restructured.gtf")


def _sirv(cell, elevel):
    return dict(
        sample=f"CL_{cell}_{elevel}_sirv",
        corpus="SIRV",
        platform="HiFi",
        hifi=True,
        genome=SIRV_GENOME,
        gtf=SIRV_GTF,
        bam=os.path.join(INPUTS, f"{cell}_{elevel}_merged_sirv_sorted.bam"),
        truth_gtf=SIRV_TRUTH_GTF,
        truth_quant=os.path.join(
            QUANT_ONLY,
            f"SIRVs/reference_data/{cell}_{elevel}_merged_sirv_sorted_groundtruth_{elevel}.tsv",
        ),
    )


def _morf(plat, hifi):
    d = os.path.join(QUANT_ONLY, f"MORFs/morf2_{plat}_merged_annot_compat_1isoform/reference_data")
    return dict(
        sample=f"morf2_{plat}",
        corpus="MORF",
        platform="HiFi" if hifi else "ONT",
        hifi=hifi,
        genome=MORF_GENOME,
        gtf=MORF_GTF,
        bam=os.path.join(INPUTS, f"morf2_{plat}_merged_annot_compat_sorted.bam"),
        # The two MORF benchmark dirs name their truth gtf differently: the
        # pacbio one carries a _pacbio suffix.  Picking the wrong one fails
        # loudly (FileNotFoundError inside the notebook), not silently.
        truth_gtf=os.path.join(
            d,
            "minigenome.UTRs_trimmed_1isoformref_expressed_pacbio.gtf" if hifi
            else "minigenome.UTRs_trimmed_1isoformref_expressed.gtf",
        ),
        truth_quant=os.path.join(d, f"morf2_{plat}_merged_annot_compat_sorted_tn_counts.tsv"),
    )


SAMPLES = [_sirv(c, e) for e in ("E1", "E2") for c in ("BT474", "HG002", "K562", "UHRR")]
SAMPLES += [_morf("pacbio", True), _morf("ont", False)]

# ------------------------------------------------- depth probe (scale test)

# Both terms in the M-step -- the fractional read counts and alpha *
# ambiguous_read_counts -- scale linearly with depth, so alpha is predicted to
# be dimensionless and its optimum should not move with library size.  That is
# a claim to test, not to assume: these are deterministic samtools subsamples
# (seed 42) of two SIRV BAMs at 50% and 25%.  Truth is unchanged, because the
# scorer renormalizes every quant file to 1e6 before comparing.

DEPTH_INPUTS = os.path.join(BASE, "depth_inputs")


def _sirv_sub(cell, elevel, pct):
    s = _sirv(cell, elevel)
    s["sample"] = f"CL_{cell}_{elevel}_sirv_sub{pct}"
    s["bam"] = os.path.join(DEPTH_INPUTS, f"{cell}_{elevel}_sub{pct}.bam")
    s["corpus"] = "SIRV_depth"
    return s


DEPTH_SAMPLES = [
    _sirv_sub(c, "E2", p) for c in ("HG002", "K562") for p in ("50", "25")
]

# The depth probe is a side experiment; it is addressable by name but is not
# part of the 10-sample corpus the recommendation is computed on.
SAMPLES_BY_NAME = {s["sample"]: s for s in SAMPLES + DEPTH_SAMPLES}

# ---------------------------------------------------------------- arms

# alpha grid spans zero plus four decades at ~half-decade spacing.  0.01 is the
# shipped default.
#
# alpha=0 is the natural lower boundary, not an arbitrary edge: the M-step adds
# alpha * ambiguous_read_count as a pseudocount, so a negative alpha subtracts
# mass hardest from the most ambiguous transcripts and drives the unnormalized
# abundance vector negative (pinned in pylib/test_EM_alpha_semantics.py).  There
# is nothing below zero to measure.
#
# The top end is NOT a boundary, so it was extended once measurement demanded
# it: the SIRV E1 mixture puts its MARD optimum at alpha=1.0, which is a
# truncation rather than an optimum, so 3 and 10 were added.  They are a
# separate list because the MORF samples were cut to five alpha points for
# wall-clock and so cannot supply them; the analysis balances each stratum
# independently rather than discarding the extension corpus-wide.
CORE_ALPHAS = [0.0, 0.0003, 0.001, 0.003, 0.01, 0.03, 0.1, 0.3, 1.0]
EXTENDED_ALPHAS = [3.0, 10.0]
ALPHAS = CORE_ALPHAS + EXTENDED_ALPHAS

# The five points MORF was originally cut to for wall-clock.  The cut was
# reversed once the SIRV strata turned out to be confounded and MORF became the
# only corpus that can answer the question; kept because the first MORF arms
# were scheduled against it.
MORF_ALPHAS = [0.0, 0.003, 0.01, 0.1, 1.0]

# alpha=1e6 is not a grid point, it is the LIMIT.  The M-step is
# theta ~ counts + alpha*ambiguous, so at large enough alpha the counts term is
# numerically irrelevant and theta IS the normalized ambiguity profile.  One run
# at this value therefore reports what the alpha->infinity prior actually
# predicts, which is NOT the profile itself: 88-92% of read mass is compatible
# with exactly one transcript and lands there regardless of theta, so the
# reported output is (unique mass, fixed) + (ambiguous mass split by profile).
LIMIT_ALPHA = 1e6

DEFAULT_ALPHA = 0.01

# 3'-end agreement weighting: on is the shipped default.  The grid is factorial
# because alpha and the 3' weights redistribute the same ambiguous reads inside
# the same EM, so an alpha optimum measured only at 3p=on is not known to hold
# at 3p=off.
THREEP = [True, False]


def alpha_token(a):
    """Filesystem- and regex-safe token for an alpha value ('0.01' -> 'a0p01').

    The '+' that %g emits in an exponent is a regex quantifier and these tokens
    go straight into the scorer's quant_pattern, so it is stripped: 1e6 becomes
    'a1e06', not 'a1e+06', which would match 'ae06' and never the real filename.
    """
    s = ("%g" % a).replace(".", "p").replace("-", "m").replace("+", "")
    return "a" + s


def arm_name(a, threep):
    return f"{alpha_token(a)}_w{1 if threep else 0}"


ARMS = [
    dict(arm=arm_name(a, w), alpha=a, threep=w, tool=f"alpha_{arm_name(a, w)}")
    for w in THREEP
    for a in ALPHAS + [LIMIT_ALPHA]
]

ARMS_BY_NAME = {x["arm"]: x for x in ARMS}
DEFAULT_ARM = arm_name(DEFAULT_ALPHA, True)
