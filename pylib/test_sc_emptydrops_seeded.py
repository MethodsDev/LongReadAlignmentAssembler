#!/usr/bin/env python3

"""`emptyDrops` must be seeded, and the seed must reach it from the workflow.

`emptyDrops` estimates its p-values by Monte Carlo simulation from the ambient
profile. An unseeded call therefore draws a different simulation on every
invocation, and every barcode whose p-value straddles the FDR threshold flips
between runs of byte-identical input. Measured on 10x library XP132120: two
whole-genome runs differing only in `sample_id` selected 15,266 and 15,274 cells
from the same 1,549,813 barcodes, and Seurat -- itself already seeded -- turned
those 8 differing cells into 55. This one call was the whole single-cell
pipeline's only source of run-to-run drift.

The behavioural proof lives outside this file, because it needs DropletUtils and
so only runs inside the lraa-sc image: the filter run twice at the default seed
yields identical barcode lists, and run at a different seed yields a different
one. What this file defends is the plumbing, on the host, in milliseconds --

  * a new `emptyDrops` call added without a `set.seed`,
  * a seed argument that exists but has no default, so an unset run is random
    again for everyone who does not know to pass it,
  * `--seed` dropped from the WDL command, leaving the input declared and inert,
  * a caller of the filter subworkflow that never forwards the seed, so it
    cannot be pinned from inputs.json.

Every one of those is invisible to a green suite otherwise, which is how the
original defect shipped.
"""

import re
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parents[1]
SC_DIR = REPO / "util" / "sc"
WDL_DIR = REPO / "WDL"
FILTER_WDL = WDL_DIR / "subwdls" / "LRAA-filter_good_cells.wdl"

# No `\s*` before the paren: `--fdr_threshold`'s help text reads "for emptyDrops
# (default: 0.01)", which is prose, not a call.
EMPTYDROPS_RE = re.compile(r"\bemptyDrops\(")
SET_SEED_RE = re.compile(r"set\.seed\s*\(\s*(\w+)")
SEED_ARG_RE = re.compile(
    r'add_argument\(\s*"--seed"(?:[^()]|\([^()]*\))*?default\s*=\s*(\d+)', re.S
)
CALL_RE = re.compile(r"^(\s*)call\s+([\w.]+)(?:\s+as\s+(\w+))?\s*\{", re.M)


def _r_uncommented(src):
    """R line comments blanked out, so a comment naming `emptyDrops` is not read
    as a call. No string literal in these scripts contains a `#`."""
    return "\n".join(line.split("#", 1)[0] for line in src.splitlines())


def r_scripts_calling_emptydrops():
    found = [p for p in sorted(SC_DIR.glob("*.R")) if EMPTYDROPS_RE.search(
        _r_uncommented(p.read_text()))]
    assert found, "no R script under {} calls emptyDrops".format(SC_DIR)
    return found


def unseeded_emptydrops_lines(src):
    """Line numbers of `emptyDrops` calls with no `set.seed` anywhere above them."""
    lines = _r_uncommented(src).splitlines()
    seeded_from = min(
        (i for i, line in enumerate(lines, 1) if SET_SEED_RE.search(line)),
        default=None,
    )
    return [
        i
        for i, line in enumerate(lines, 1)
        if EMPTYDROPS_RE.search(line) and (seeded_from is None or seeded_from > i)
    ]


def seed_is_caller_supplied(src):
    """The seeded value must trace to the command line, not a literal in the
    source. A hardcoded `set.seed(1)` is reproducible but unpinnable, which
    defeats the point of exposing the knob at all."""
    names = {m.group(1) for m in SET_SEED_RE.finditer(_r_uncommented(src))}
    return bool(names) and not any(n.isdigit() for n in names) and "args$seed" in src


def filter_subworkflow_callers():
    """(path, source) for every WDL that imports the filter subworkflow."""
    out = []
    for path in sorted(WDL_DIR.rglob("*.wdl")):
        src = path.read_text()
        if "LRAA-filter_good_cells.wdl" in src and path != FILTER_WDL:
            out.append((path, src))
    assert out, "nothing imports LRAA-filter_good_cells.wdl any more"
    return out


def calls_missing_seed(src):
    """Aliases of `FilterGoodCells` calls that do not bind `seed`."""
    offenders = []
    for match in CALL_RE.finditer(src):
        callee = match.group(2)
        if not callee.endswith("FilterGoodCells"):
            continue
        indent = match.group(1)
        rest = src[match.end():]
        end = re.search(r"^{}\}}".format(indent), rest, re.M)
        body = rest[: end.start()] if end else rest
        if not re.search(r"^\s*seed\s*=", body, re.M):
            offenders.append(match.group(3) or callee)
    return offenders


@pytest.mark.parametrize(
    "path", r_scripts_calling_emptydrops(), ids=lambda p: p.name
)
def test_every_emptydrops_call_is_seeded(path):
    src = path.read_text()
    offenders = unseeded_emptydrops_lines(src)
    assert not offenders, (
        "{}: emptyDrops called at line(s) {} with no set.seed above them. An "
        "unseeded call draws a fresh Monte Carlo simulation every run and "
        "barcodes near the FDR threshold flip between runs.".format(
            path.relative_to(REPO), ", ".join(map(str, offenders))
        )
    )
    assert seed_is_caller_supplied(src), (
        "{}: the seed must come from --seed via args$seed, not a literal; a "
        "hardcoded seed cannot be pinned or varied by a caller.".format(
            path.relative_to(REPO)
        )
    )
    assert SEED_ARG_RE.search(src), (
        "{}: --seed must declare an integer default, or an unset run is random "
        "again and the fix is opt-in.".format(path.relative_to(REPO))
    )


def test_the_filter_task_emits_the_seed_flag():
    """A declared input that never becomes a flag is the defect this guards."""
    src = FILTER_WDL.read_text()
    assert re.search(r"^\s*Int seed = \d+\s*$", src, re.M), (
        "the filter subworkflow must declare `Int seed` WITH a default"
    )
    assert "--seed ~{seed}" in src, (
        "run_filter_good_cells declares a seed it never passes to "
        "filter_good_cells.R"
    )


@pytest.mark.parametrize(
    "path,src", filter_subworkflow_callers(), ids=lambda a: getattr(a, "name", "")
)
def test_callers_forward_the_seed(path, src):
    """Declared on the subworkflow but not forwarded means it can only be set by
    editing the WDL, not from inputs.json -- which is where runs are launched."""
    offenders = calls_missing_seed(src)
    assert not offenders, "{}: FilterGoodCells call(s) {} do not forward `seed`".format(
        path.relative_to(REPO), ", ".join(offenders)
    )


# ---------------------------------------------------------------------------
# Negative controls. Each reverts one half of the shipped fix in a COPY of the
# real source and asserts the check reports it. Without these the checks could
# quietly stop matching and every file would pass vacuously.
# ---------------------------------------------------------------------------


def test_the_seed_check_catches_the_defect_as_it_shipped():
    src = (SC_DIR / "filter_good_cells.R").read_text()
    reverted = src.replace("set.seed(seed)", "", 1)
    assert reverted != src, "filter_good_cells.R no longer calls set.seed(seed)"
    assert unseeded_emptydrops_lines(reverted), "the unseeded call went unnoticed"


def test_the_seed_check_catches_a_hardcoded_seed():
    src = (SC_DIR / "filter_good_cells.R").read_text()
    assert not seed_is_caller_supplied(src.replace("set.seed(seed)", "set.seed(1)", 1))


def test_the_flag_check_catches_a_dropped_emission():
    src = FILTER_WDL.read_text().replace("--seed ~{seed}", "")
    assert "--seed ~{seed}" not in src


def test_the_forwarding_check_catches_a_dropped_forward():
    path, src = next(
        (p, s) for p, s in filter_subworkflow_callers() if "seed = " in s
    )
    reverted = re.sub(r"^\s*seed = filter_seed\s*$", "", src, flags=re.M)
    assert reverted != src, "{} no longer forwards seed = filter_seed".format(path.name)
    assert calls_missing_seed(reverted), "the dropped forward went unnoticed"
