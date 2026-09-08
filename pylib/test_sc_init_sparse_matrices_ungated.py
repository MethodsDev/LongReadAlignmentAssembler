#!/usr/bin/env python3

"""The initial pass's single-cell matrices must not be gated on clustering.

`LRAA-singlecell.wdl` built its `init_sc_*` sparse matrices INSIDE the
`run_clustering_phase` conditional, which made them an intermediate of Seurat
rather than a deliverable of the initial LRAA run. The consequence was silent
and easy to mistake for a code failure: supply
`precomputed_cluster_assignments_tsv` -- which says only "do not re-cluster" --
and the run emitted `init_gtf`, `init_quant_expr` and `init_quant_tracking` while
every `init_sc_*` output came back null, with the tracking file those matrices
are derived from sitting in the same output set.

Nothing else catches this. Both builds are optional-output calls in a
conditional, so `miniwdl check` and womtool are satisfied either way -- they
prove syntax and types, not which tasks a given input set reaches -- and a run
that skips them is a clean success.

Two layers of assertion, because either alone has a hole:

  STRUCTURE (`test_*_gated*`) -- neither build sits under `run_clustering_phase`,
  while the empty-droplet filter and Seurat still do. Catches a future edit that
  re-nests the builds under the clustering gate. Blind, though, to the gate being
  folded into a boolean the builds already depend on.

  SEMANTICS (`test_gate_truth_table`, `test_*_branch`) -- the workflow's own gate
  expressions, evaluated by miniwdl over concrete input scenarios, must select
  exactly one matrix source per scenario and skip clustering exactly when
  assignments were supplied. Catches the folded-gate case and any change to the
  boolean chain, since it re-derives `run_initial_phase`,
  `run_clustering_phase`, `use_sc_sparse_from_shards` and
  `init_quant_tracking_file` from the file rather than restating them.

Static parse and expression evaluation only: no docker, no fixtures, no
execution.
"""

import pytest

WDL = pytest.importorskip("WDL", reason="miniwdl not installed")

from pathlib import Path

from WDL import Env, StdLib, Value

REPO = Path(__file__).resolve().parents[1]
WDL_DIR = REPO / "WDL"

SC_WDL = "LRAA-singlecell.wdl"

# The two mutually exclusive sources of the initial matrices: the streaming merge
# of per-contig shard artifacts, and the single library-wide pass over the merged
# tracking file.
SHARD_MERGE = "merge_sc_from_shards"
LIBRARY_BUILD = "build_sc_from_init_tracking"
INIT_SPARSE_CALLS = (LIBRARY_BUILD, SHARD_MERGE)

# Steps that legitimately stay behind the clustering gate.
CLUSTERING_CALLS = ("filter_good_cells", "cluster_cells")

CLUSTERING_GATE = "run_clustering_phase"

_STDLIB = StdLib.Base("1.0")


def _workflow(rel=SC_WDL):
    return WDL.load(str(WDL_DIR / rel), path=[str(WDL_DIR)]).workflow


def _idents(expr):
    """Every identifier named anywhere in ``expr``."""
    names = []
    if isinstance(expr, WDL.Expr.Ident):
        names.append(str(expr.name))
    for child in getattr(expr, "children", []) or []:
        if isinstance(child, WDL.Expr.Base):
            names.extend(_idents(child))
    return names


def _elements(node, conditions=()):
    """Yield (element, enclosing if/scatter expressions) for the whole body."""
    for child in getattr(node, "body", []) or []:
        if isinstance(child, (WDL.Tree.Conditional, WDL.Tree.Scatter)):
            yield from _elements(child, conditions + (child.expr,))
        else:
            yield child, conditions


def _call_conditions(workflow, call_name):
    for element, conditions in _elements(workflow):
        if isinstance(element, WDL.Tree.Call) and element.name == call_name:
            return conditions
    raise AssertionError(f"{SC_WDL}: no call named {call_name}")


def _call_condition_idents(workflow, call_name):
    names = []
    for condition in _call_conditions(workflow, call_name):
        names.extend(_idents(condition))
    return names


# --------------------------------------------------------------------------
# Structure
# --------------------------------------------------------------------------


@pytest.mark.parametrize("call_name", INIT_SPARSE_CALLS)
def test_init_sparse_build_not_gated_on_clustering(call_name):
    conditions = _call_condition_idents(_workflow(), call_name)
    assert CLUSTERING_GATE not in conditions, (
        f"{SC_WDL}: {call_name} is nested under {CLUSTERING_GATE}, so supplying "
        "precomputed_cluster_assignments_tsv drops the initial matrices"
    )


@pytest.mark.parametrize("call_name", CLUSTERING_CALLS)
def test_clustering_steps_remain_gated(call_name):
    conditions = _call_condition_idents(_workflow(), call_name)
    assert CLUSTERING_GATE in conditions, (
        f"{SC_WDL}: {call_name} runs regardless of {CLUSTERING_GATE}, so supplied "
        "cluster assignments would be recomputed"
    )


def test_shard_merge_requires_the_initial_pass_to_have_run():
    """The shard merge's input is `LRAA_init.scShardSparse`.

    With the initial pass skipped there are no shard artifacts, so a gate that
    ignores `run_initial_phase` merges an empty shard list and publishes matrices
    built from nothing. That is what the old
    `select_first([LRAA_init.scShardSparse, []])` default made silent.
    """
    workflow = _workflow()

    decls = [
        element
        for element, _ in _elements(workflow)
        if isinstance(element, WDL.Tree.Decl)
        and element.name == "use_sc_sparse_from_shards"
    ]
    assert len(decls) == 1, f"{SC_WDL}: expected one use_sc_sparse_from_shards decl"

    gate = _idents(decls[0].expr)
    for required in ("sc_sparse_from_shards", "scattering_init", "run_initial_phase"):
        assert required in gate, (
            f"{SC_WDL}: use_sc_sparse_from_shards ignores {required}; the merge would "
            "run on an absent or empty shard list"
        )

    assert "use_sc_sparse_from_shards" in _call_condition_idents(workflow, SHARD_MERGE)


# Workflow output -> the declaration it must publish. The outputs existed before
# this change and were null half the time, so a name check alone proves nothing:
# what matters is that each one still binds the top-level declaration the builds
# now feed, rather than some other expression reintroducing a gate of its own.
INIT_SPARSE_OUTPUTS = {
    "init_sc_gene_sparse_tar_gz": "init_sc_gene_tgz",
    "init_sc_isoform_sparse_tar_gz": "init_sc_isoform_tgz",
    "init_sc_splice_pattern_sparse_tar_gz": "init_sc_splice_tgz",
    "init_sc_gene_transcript_splicehash_mapping": "init_sc_mapping",
}


def test_init_matrices_are_workflow_outputs():
    workflow = _workflow()
    outputs = {output.name: output for output in workflow.outputs}

    top_level = {
        element.name
        for element, conditions in _elements(workflow)
        if isinstance(element, WDL.Tree.Decl) and not conditions
    }

    for name, decl_name in INIT_SPARSE_OUTPUTS.items():
        assert name in outputs, f"{SC_WDL}: {name} is not a workflow output"
        assert str(outputs[name].expr) == decl_name, (
            f"{SC_WDL}: output {name} publishes {outputs[name].expr!s}, not {decl_name}"
        )
        assert decl_name in top_level, (
            f"{SC_WDL}: {decl_name} is not declared at workflow scope, so it is "
            "conditional on whatever block encloses it"
        )


# --------------------------------------------------------------------------
# Semantics: the workflow's own gate expressions, evaluated
# --------------------------------------------------------------------------

_FILE = Value.File("/dev/null")
_NULL = Value.Null()

# Booleans the scenarios below assert on, all derived from the file.
_DERIVED = (
    "run_initial_phase",
    CLUSTERING_GATE,
    "use_sc_sparse_from_shards",
    "init_quant_tracking_file",
)


def _scenario_inputs(
    *,
    precomputed_init_gtf=False,
    precomputed_init_quant_tracking=False,
    precomputed_cluster_assignments_tsv=False,
    scattering_init="by_chromosome",
    sc_sparse_from_shards=True,
):
    return {
        "scattering_init": Value.String(scattering_init),
        "sc_sparse_from_shards": Value.Boolean(sc_sparse_from_shards),
        "enable_filter_good_cells": Value.Boolean(True),
        "precomputed_init_gtf": _FILE if precomputed_init_gtf else _NULL,
        "precomputed_init_quant_tracking": (
            _FILE if precomputed_init_quant_tracking else _NULL
        ),
        "precomputed_cluster_assignments_tsv": (
            _FILE if precomputed_cluster_assignments_tsv else _NULL
        ),
    }


def _fixpoint(inputs, required=()):
    """Bind ``inputs``, then evaluate the TOP-LEVEL declarations that resolve.

    Top-level only, and that is a correctness bound rather than a shortcut: a
    declaration inside `if (...)` is evaluated by the engine ONLY when its
    condition holds, so evaluating one here would be running code the engine
    would have skipped -- and `gene_sparse_for_clustering`'s
    `select_first([...])` raises on a false branch's nulls, which would fail a
    scenario for a reason that has nothing to do with it. Every gate asserted on
    is declared at workflow scope, so nothing is lost by refusing to leave it.
    """
    env = Env.Bindings()
    for name, value in inputs.items():
        env = env.bind(name, value)

    workflow = _workflow()
    decls = [
        element
        for element, conditions in _elements(workflow)
        if isinstance(element, WDL.Tree.Decl) and not conditions
    ]

    # Declarations depending on calls this pass does not bind stay unbound; the
    # caller asserts that everything it needs resolved.
    bound = set(inputs)
    progress = True
    while progress:
        progress = False
        for decl in decls:
            if decl.name in bound:
                continue
            if all(name in bound for name in _idents(decl.expr)):
                env = env.bind(decl.name, decl.expr.eval(env, _STDLIB))
                bound.add(decl.name)
                progress = True

    missing = [name for name in required if name not in bound]
    assert not missing, (
        f"{SC_WDL}: could not derive {missing} from the workflow's declarations; "
        "the gate chain has changed shape"
    )
    return env


def _resolve(**scenario):
    """Evaluate the workflow's top-level declarations for one input scenario.

    TWO PASSES, so that nothing here restates the workflow's own logic. Pass one
    binds the scenario's inputs alone and lets the file derive
    `run_initial_phase`. Pass two then binds `LRAA_init`'s outputs the way
    miniwdl does -- present iff that call ran, since it sits inside
    `if (run_initial_phase)` -- using the value pass one produced. A stub that
    were unconditional would make a skipped initial pass look as though it left a
    tracking file behind, which is the distinction these scenarios turn on; a
    stub computed from a local copy of `has_precomputed_init` would assert
    against this file's model of the gate instead of the gate.
    """
    inputs = _scenario_inputs(**scenario)

    first = _fixpoint(inputs, required=("run_initial_phase",))
    init_output = _FILE if first["run_initial_phase"].value else _NULL

    return _fixpoint(
        {
            **inputs,
            "LRAA_init.mergedQuantTracking": init_output,
            "LRAA_init.mergedQuantExpr": init_output,
            "LRAA_init.mergedGTF": init_output,
        },
        required=_DERIVED,
    )


def _reached(env, workflow, call_name):
    """True when every conditional enclosing ``call_name`` holds in ``env``."""
    for condition in _call_conditions(workflow, call_name):
        if not condition.eval(env, _STDLIB).value:
            return False
    return True


# (label, scenario, init pass runs, clustering runs, matrix source)
_SCENARIOS = [
    ("stock_basic", {}, True, True, SHARD_MERGE),
    # THE branch this test exists for: assignments supplied, initial pass still
    # runs, so the matrices are still built and only clustering is skipped.
    (
        "precomputed_clusters_only",
        {"precomputed_cluster_assignments_tsv": True},
        True,
        False,
        SHARD_MERGE,
    ),
    # Initial pass skipped, tracking supplied: no shards exist, so the
    # library-wide build over the precomputed tracking file is the only source.
    (
        "precomputed_init_and_clusters",
        {
            "precomputed_init_gtf": True,
            "precomputed_init_quant_tracking": True,
            "precomputed_cluster_assignments_tsv": True,
        },
        False,
        False,
        LIBRARY_BUILD,
    ),
    # The documented null case: nothing to build from, so neither source runs.
    (
        "precomputed_gtf_and_clusters_no_tracking",
        {"precomputed_init_gtf": True, "precomputed_cluster_assignments_tsv": True},
        False,
        False,
        None,
    ),
    # One whole-genome invocation emits no per-contig shards, whatever
    # sc_sparse_from_shards asks for.
    ("scattering_off", {"scattering_init": "off"}, True, True, LIBRARY_BUILD),
]


@pytest.mark.parametrize(
    "label,scenario,init_runs,clustering_runs,source",
    _SCENARIOS,
    ids=[entry[0] for entry in _SCENARIOS],
)
def test_gate_truth_table(label, scenario, init_runs, clustering_runs, source):
    workflow = _workflow()
    env = _resolve(**scenario)

    assert env["run_initial_phase"].value is init_runs, label
    assert env[CLUSTERING_GATE].value is clustering_runs, label

    reached = {name: _reached(env, workflow, name) for name in INIT_SPARSE_CALLS}

    if source is None:
        assert not any(reached.values()), (
            f"{label}: a matrix build runs with no init tracking file: {reached}"
        )
    else:
        assert reached[source], f"{label}: expected {source} to run, got {reached}"
        assert sum(reached.values()) == 1, (
            f"{label}: the two matrix sources are not exclusive: {reached}"
        )

    for call_name in CLUSTERING_CALLS:
        assert _reached(env, workflow, call_name) is clustering_runs, (
            f"{label}: {call_name} reachability disagrees with {CLUSTERING_GATE}"
        )


def test_precomputed_clusters_branch_keeps_matrices_and_skips_clustering():
    """The reported defect, stated as one assertion rather than a table row."""
    workflow = _workflow()
    env = _resolve(precomputed_cluster_assignments_tsv=True)

    assert env["run_initial_phase"].value is True
    assert env[CLUSTERING_GATE].value is False
    assert not isinstance(env["init_quant_tracking_file"], Value.Null)

    assert _reached(env, workflow, SHARD_MERGE), (
        "precomputed cluster assignments must not suppress the initial matrices"
    )
    for call_name in CLUSTERING_CALLS:
        assert not _reached(env, workflow, call_name), (
            f"{call_name} must not run when cluster assignments were supplied"
        )


def test_precomputed_tracking_forces_the_library_wide_build():
    """A skipped initial pass has no shards, so the shard merge must not run."""
    workflow = _workflow()
    env = _resolve(
        precomputed_init_gtf=True,
        precomputed_init_quant_tracking=True,
        precomputed_cluster_assignments_tsv=True,
    )

    assert env["run_initial_phase"].value is False
    assert env["use_sc_sparse_from_shards"].value is False
    assert not isinstance(env["init_quant_tracking_file"], Value.Null)

    assert _reached(env, workflow, LIBRARY_BUILD)
    assert not _reached(env, workflow, SHARD_MERGE), (
        "the shard merge would run on an empty shard list"
    )


def test_gene_symbol_sources_follow_the_filter_having_run():
    """`enable_filter_good_cells` is an input; the filter RUNNING is a gate.

    filter_good_cells lives behind the clustering gate, so with cluster
    assignments supplied it does not run while that input keeps its default of
    true. Selecting the filtered matrices on the input alone therefore yielded a
    null at all three levels, and the `defined(...)` conjunction guarding
    add_gene_symbols turned that into ref_annot_gtf_source_gene_symbols doing
    nothing at all -- in precisely the configuration that now has initial
    matrices to annotate. So the three selectors must test both conditions,
    which `have_filtered_matrices` is.
    """
    workflow = _workflow()

    decls = [
        element
        for element, conditions in _elements(workflow)
        if isinstance(element, WDL.Tree.Decl)
        and element.name == "have_filtered_matrices"
        and not conditions
    ]
    assert len(decls) == 1, f"{SC_WDL}: expected one have_filtered_matrices decl"

    gate = _idents(decls[0].expr)
    assert CLUSTERING_GATE in gate and "enable_filter_good_cells" in gate, (
        f"{SC_WDL}: have_filtered_matrices is {gate}; the filter's OUTPUTS exist "
        f"only when both {CLUSTERING_GATE} and enable_filter_good_cells hold"
    )

    for level in ("gene", "isoform", "splice_pattern"):
        name = f"{level}_sparse_for_symbols"
        decl = next(
            element
            for element, _ in _elements(workflow)
            if isinstance(element, WDL.Tree.Decl) and element.name == name
        )
        referenced = _idents(decl.expr)
        assert "have_filtered_matrices" in referenced, (
            f"{SC_WDL}: {name} chooses the filtered matrices on {referenced}, not "
            "on the filter having run"
        )
        assert f"init_sc_{'splice' if level == 'splice_pattern' else level}_tgz" in (
            referenced
        ), f"{SC_WDL}: {name} does not fall back to the initial matrices"
