#!/usr/bin/env python3

"""`uniq_reads`/`uniq_FSM_reads` in quant.expr, and `is_FSM` in quant.tracking.

Three properties, all of which a plausible implementation gets wrong:

  1. A reported unique read is EXCLUSIVELY ASSIGNABLE -- compatible with exactly one
     isoform. Not "EM gave it almost all the mass", which is what
     unique_read_report_min_frac tested: a read compatible with two isoforms reaches
     1.0 whenever the competitor holds zero mass under the current theta.
  2. `uniq_FSM_reads` is the subset of those whose intron chain is EXACTLY the
     isoform's, so it can never exceed `uniq_reads`.
  3. A monoexonic model has no intron chain to match, so its count is 0 -- a
     measurement, not a placeholder.

The expr columns are aggregates of the per-read `is_FSM` the tracking file carries, so
the two files must agree transcript by transcript. That cross-file identity is asserted
here because it is what makes either column trustworthy: checking it by hand on a whole
run is possible but not repeatable, and an assertion is what keeps it true.

The fractional assignment map is built directly rather than by running EM, so the test
controls exactly which isoforms each multipath is compatible with -- which is the
variable under test and is not otherwise observable.
"""

import io
import sys
from pathlib import Path

import pytest

REPO_ROOT = Path(__file__).resolve().parents[1]
if str(REPO_ROOT / "pylib") not in sys.path:
    sys.path.insert(0, str(REPO_ROOT / "pylib"))

import LRAA_Globals
from GenomeFeature import Exon, Intron
from MultiPath import MultiPath
from Quantify import Quantify
from Splice_graph import Splice_graph

CONTIG = "chr1"
#          0            1            2            3
EXONS = [(100, 150), (200, 250), (400, 450), (600, 650)]


@pytest.fixture(autouse=True)
def _isolate():
    keys = ("weight_reads_by_3prime_agreement", "num_total_reads")
    saved = {k: LRAA_Globals.config[k] for k in keys}
    LRAA_Globals.config["weight_reads_by_3prime_agreement"] = False
    # get_TPM divides by this; unset in a bare process, and RPM_total_reads is not
    # what this test is about.
    LRAA_Globals.config["num_total_reads"] = 1_000_000
    LRAA_Globals.reset_read_weight_registry()
    try:
        yield
    finally:
        LRAA_Globals.config.update(saved)
        LRAA_Globals.reset_read_weight_registry()


def _graph():
    """Exons AND introns as graph nodes.

    The chain test filters a path for its "I:" nodes, so a graph of exons alone makes
    every path look monoexonic, no read could ever be an FSM, and the FSM assertions
    below would pass vacuously.
    """
    Exon.reset_counter()
    sg = Splice_graph()
    exons = [Exon(CONTIG, lend, rend, "+", 1) for lend, rend in EXONS]
    for e in exons:
        sg._node_id_to_node[e.get_id()] = e
    # 0->1, 1->2, and the long 1->3 skip that makes t_other a different chain
    introns = [Intron(CONTIG, lend, rend, "+", 10)
               for lend, rend in [(151, 199), (251, 399), (251, 599)]]
    for i in introns:
        sg._node_id_to_node[i.get_id()] = i
        sg._intron_objs["{}:{}".format(*i.get_coords())] = i
    sg._contig_acc = CONTIG
    sg._contig_strand = "+"
    return sg, [e.get_id() for e in exons], [i.get_id() for i in introns]


def _tx(sg, path, tid, gene, counts):
    t = MultiPath(sg, [list(path)]).toTranscript()
    t.set_gene_id(gene)
    t.set_transcript_id(tid)
    t.set_read_counts_assigned(counts)
    t.set_isoform_fraction(1.0)
    return t


def _mp(sg, path, names):
    return MultiPath(
        sg, [list(path)],
        read_types={"PacBio"}, read_names=set(names), read_count=len(names),
    )


def _report(transcripts, assignment_map):
    q = Quantify(run_EM=False, max_EM_iterations=1)
    expr, track = io.StringIO(), io.StringIO()
    q.report_quant_results(transcripts, assignment_map, expr, track)
    exprs = {}
    for line in expr.getvalue().strip().split("\n"):
        f = line.split("\t")
        # (uniq_reads, uniq_FSM_reads): col 2 and the trailing column
        exprs[f[1]] = (int(f[2]), int(f[-1]))
    tracking = []
    for line in track.getvalue().strip().split("\n"):
        if not line:
            continue
        f = line.split("\t")
        tracking.append({"tid": f[1], "read": f[5], "is_unique": f[8], "is_FSM": f[9]})
    return exprs, tracking


def _fixture():
    """Three isoforms, and three kinds of read.

    The shared read is given fraction EXACTLY 1.0 for t_multi and 0.0 for t_other.
    That is deliberate and is what makes this test able to catch a regression: the old
    rule was `frac >= unique_read_report_min_frac` with that key defaulting to 1.0, so
    a shared read at 0.9999 would have been excluded by BOTH rules and the test would
    prove nothing. At exactly 1.0 the old rule counts it as unique for t_multi while
    exclusivity refuses it -- the `frac1_but_shared` case, which a competitor holding
    zero mass produces.

    exclusive_fsm  -- only t_multi can take it, chain identical            -> unique + FSM
    exclusive_prefix -- only t_multi can take it, chain is a PREFIX        -> unique, NOT FSM
    shared         -- t_multi AND t_other can take it                     -> neither
    mono_read      -- only t_mono can take it, no chain at all            -> unique, NOT FSM
    """
    sg, E, I = _graph()
    chain_multi = [E[0], I[0], E[1], I[1], E[2]]   # introns 0,1
    chain_other = [E[0], I[0], E[1], I[2], E[3]]   # intron 0 then the long skip
    prefix_path = [E[0], I[0], E[1]]               # intron 0 only

    t_multi = _tx(sg, chain_multi, "t_multi", "g1", 5.0)
    t_other = _tx(sg, chain_other, "t_other", "g1", 2.0)
    t_mono = _tx(sg, [E[0]], "t_mono", "g2", 3.0)

    exclusive_fsm = _mp(sg, chain_multi, ["f1", "f2", "f3"])
    exclusive_prefix = _mp(sg, prefix_path, ["p1"])
    shared = _mp(sg, prefix_path, ["s1", "s2"])
    mono_read = _mp(sg, [E[0]], ["m1", "m2"])

    t_multi.add_multipaths_evidence_assigned([exclusive_fsm, exclusive_prefix, shared])
    t_other.add_multipaths_evidence_assigned([shared])
    t_mono.add_multipaths_evidence_assigned([mono_read])

    # Built from each transcript's ACTUAL assigned evidence, because toTranscript()
    # also registers the transcript's own defining multipath as evidence and
    # report_quant_results looks up every one of them. That path carries no read names,
    # so it contributes no rows and no counts -- but it must be present in the map.
    intended = {
        "t_multi": {exclusive_fsm: 1.0, exclusive_prefix: 1.0, shared: 1.0},
        "t_other": {shared: 0.0},
        "t_mono": {mono_read: 1.0},
    }
    amap = {}
    for t in (t_multi, t_other, t_mono):
        tid = t.get_transcript_id()
        amap[tid] = {
            mp: intended[tid].get(mp, 0.0)
            for mp in t.get_multipaths_evidence_assigned()
        }
    return [t_multi, t_other, t_mono], amap


def test_unique_means_exclusively_assignable_not_well_allocated():
    """`shared` holds exactly 1.0 for t_multi yet t_other can also take it: not unique.

    Fails if the column reverts to `frac >= unique_read_report_min_frac`, which at its
    default of 1.0 counts this read for t_multi.
    """
    transcripts, amap = _fixture()
    exprs, _ = _report(transcripts, amap)
    assert exprs["t_multi"][0] == 4, (
        "expected the 3 FSM reads plus the 1 prefix read, all exclusive to t_multi; "
        "a count of 6 means the 2 shared reads were counted as unique"
    )
    assert exprs["t_other"][0] == 0, (
        "t_other's only read is shared with t_multi, so it has no unique reads"
    )


def test_uniq_FSM_is_the_exact_chain_subset():
    """Of t_multi's 4 unique reads only the 3 with its exact chain are FSM."""
    transcripts, amap = _fixture()
    exprs, _ = _report(transcripts, amap)
    assert exprs["t_multi"][1] == 3, (
        "the prefix read shares no full chain with t_multi and must not be FSM; "
        "a count of 4 means compatibility is being accepted for identity"
    )


def test_uniq_FSM_never_exceeds_uniq_reads():
    """Invariant, not an observation: FSM-unique is a subset of unique."""
    transcripts, amap = _fixture()
    exprs, _ = _report(transcripts, amap)
    for tid, (uniq, uniq_fsm) in exprs.items():
        assert uniq_fsm <= uniq, f"{tid}: uniq_FSM_reads {uniq_fsm} > uniq_reads {uniq}"


def test_monoexonic_model_reports_zero_FSM():
    """No intron chain to match. Its unique reads are real; its FSM count is 0."""
    transcripts, amap = _fixture()
    exprs, _ = _report(transcripts, amap)
    assert exprs["t_mono"][0] == 2, "t_mono's 2 reads are exclusively its own"
    assert exprs["t_mono"][1] == 0, (
        "a monoexonic model has no chain, so no read can be a full splice match to it"
    )


def test_tracking_is_FSM_aggregates_to_the_expr_column():
    """The cross-file identity: expr counts are aggregates of the per-read flags.

    Both expr columns are recomputable from the tracking file alone by reading the
    explicit is_unique/is_FSM flags. If they ever disagree,
    one of the two files is wrong and neither can be trusted.
    """
    transcripts, amap = _fixture()
    exprs, tracking = _report(transcripts, amap)

    # Keyed on the explicit is_unique flag, NOT on how many rows a read produced:
    # that inference is invalid on the --oversimplify paths, which emit one row for a
    # read whose best overlap tied and which is rightly not counted as unique.
    recomputed = {tid: [0, 0] for tid in exprs}
    for row in tracking:
        if row["is_unique"] != "1":
            continue
        recomputed[row["tid"]][0] += 1
        if row["is_FSM"] == "1":
            recomputed[row["tid"]][1] += 1

    for tid in exprs:
        assert tuple(recomputed[tid]) == exprs[tid], (
            f"{tid}: tracking implies {tuple(recomputed[tid])} but expr reports "
            f"{exprs[tid]} -- the aggregate and the per-read flag disagree"
        )
