#!/usr/bin/env python3

"""A declared input that never becomes a flag is the defect class this repo keeps
shipping. This is the static check for it.

Nine instances so far, all the same shape -- an input is declared, documented,
validated, and then silently dropped one layer before the command line:

  * ``--bam_for_sg`` was absent from ``LRAA``'s ``ChunkedRun.default_args(...)``
    enumeration, so every chunk manufactured its own splice-graph evidence.
  * six chunk resource knobs in ``WDL/LRAA.wdl`` stopped one call short of the
    subworkflow that consumes them.
  * ``String? oversimplify`` was declared by ``process_chunk`` in
    ``WDL/subwdls/LRAA_chunk_scatter.wdl``, passed to it by the workflow, and
    never emitted onto the leaf's ``ChunkedRun.py`` command line.
  * ``--oversimplify`` was missing from that same ``default_args`` enumeration,
    so it was inert on EVERY WDL path -- ``subwdls/LRAA_runner.wdl`` always
    passes ``--chunk``.
  * ``main_chromosomes`` validated with ``scattering=off`` and was then dropped,
    because the ``LRAA_direct`` call never set the ``contig`` input its callee
    had always declared -- so a run asked for chr21 and got the whole genome.
  * five negated flags -- ``--no_weight_reads_by_3prime_agreement``,
    ``--no_filter_internal_priming``, ``--no_infer_TSS``, ``--no_infer_PolyA``
    and ``--no_use_community_clustering`` -- were missing from
    ``_NEGATED_CONFIG_FLAGS``, so no chunk worker ever received them.

Every one of them produced a run that completed, wrote plausible output, and
reported nothing. None was caught by a test, by ``miniwdl check`` (which reports
``UnusedDeclaration`` only as an advisory and still exits 0), or by ``womtool
validate`` (an unreferenced declaration is legal WDL). So the check has to be
here, and it has to be an assertion rather than a lint.

Four surfaces, one rule -- a declared input must reach a flag:

  A. every WDL task that actually RUNS ``LRAA`` or ``ChunkedRun.py`` must
     reference each of its declared inputs from its command block, directly or
     through its own task-level declarations.
  B. every CLI flag that BOTH ``LRAA`` and ``pylib/ChunkedRun.py`` declare must
     appear in ``LRAA``'s ``ChunkedRun.default_args(...)`` call. That call is an
     allowlist: an unlisted flag is parsed off the driver's command line and
     dropped, and the chunked arm silently uses ChunkedRun's own default.
  C. two calls to the SAME subworkflow in one file must set the same inputs.
     Divergence is how ``LRAA_direct`` came to ignore ``main_chromosomes`` while
     ``LRAA_scatter``, three hundred lines up, passed ``contig`` correctly.
  D. every ``config[KEY] = not args.DEST`` in ``LRAA`` must have its pair in
     ``_NEGATED_CONFIG_FLAGS``. The chunk-worker override collector routes on
     ``dest in config``, and a negated flag's dest is never a config key, so an
     unmapped pair is dropped -- which is where
     ``--no_weight_reads_by_3prime_agreement`` and four others were found.

Purely static: parses text and ASTs, imports nothing from the pipeline, needs no
cluster, no docker and no fixture data.
"""

import ast
import re
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parents[1]
WDL_DIR = REPO / "WDL"
LRAA_DRIVER = REPO / "LRAA"
CHUNKED_RUN = REPO / "pylib" / "ChunkedRun.py"

TASK_RE = re.compile(r"^task\s+(\w+)\s*\{", re.M)
BOUNDARY_RE = re.compile(r"^(task|workflow)\s+\w+\s*\{", re.M)
INPUT_BLOCK_RE = re.compile(r"^(\s*)input\s*\{(.*?)^\1\}", re.M | re.S)
COMMAND_RE = re.compile(r"^\s*command\s*<<<(.*?)^\s*>>>", re.M | re.S)
RUNTIME_RE = re.compile(r"^(\s*)runtime\s*\{(.*?)^\1\}", re.M | re.S)
OUTPUT_RE = re.compile(r"^(\s*)output\s*\{(.*?)^\1\}", re.M | re.S)
# A WDL declaration: a type (possibly optional/compound) then an identifier.
DECL_RE = re.compile(
    r"^[ \t]*(?:Array|Map|Pair)\[[^\]]*\][?+]*[ \t]+(\w+)"
    r"|^[ \t]*(?:String|File|Int|Float|Boolean)\??[ \t]+(\w+)",
    re.M,
)
PLACEHOLDER_RE = re.compile(r"~\{(.*?)\}", re.S)
IDENT_RE = re.compile(r"\b([A-Za-z_]\w*)\b")
# The driver / module invoked as a command word at the start of a shell line,
# with or without a leading path. `$(LRAA --version)` mid-line is deliberately
# NOT a match: a version probe emits no pipeline flags and has no inputs to drop.
INVOCATION_RE = re.compile(r"^[ \t]*(?:\S*/)?(?:LRAA|ChunkedRun\.py)[ \t]", re.M)
CALL_RE = re.compile(r"^(\s*)call\s+([\w.]+)(?:\s+as\s+(\w+))?\s*\{", re.M)


# ---------------------------------------------------------------------------
# Allowlists. Each entry is an input consumed somewhere a command line cannot
# see it. Keep this SHORT: every addition is a hole, and a long allowlist turns
# the whole file decorative.
# ---------------------------------------------------------------------------

# Consumed by the `runtime` block, which is not a command line at all.
RUNTIME_ONLY_INPUTS = {
    "docker",  # the container image
    "memoryGB",  # task sizing
    "diskSizeGB",  # task sizing
    "chunkPreemptible",  # per-leaf preemptibility
}

# Calls to one callee that legitimately differ, because the WDL layer above them
# forbids the combination. Keyed (file, callee) -> {input: why}.
DIVERGENT_CALL_INPUTS = {
    ("LRAA.wdl", "LRAA_runner"): {
        # An unscattered run has no shard, so there is no index to name one.
        "shardno": "identifies a chromosome shard; scattering=off has none",
        # validate_scattering REFUSES a plan with scattering=off (LRAA.wdl):
        # a plan is geometry shared between sibling runs and off has no siblings.
        "chunk_plan": "refused with scattering=off by validate_scattering",
        # validate_scattering REQUIRES scattering=off for a region, so the
        # scattered call cannot carry one.
        "region": "requires scattering=off by validate_scattering",
    },
}


# ---------------------------------------------------------------------------
# Surface A -- WDL task inputs vs. the command block that runs the pipeline
# ---------------------------------------------------------------------------


def _strip_comments(text):
    """WDL line comments. Nothing in this repo's WDL puts a `#` inside a string
    in a declaration or call block, and the command blocks are read whole rather
    than stripped, so a line-wise strip is sufficient and does not lose a flag.
    """

    return re.sub(r"#[^\n]*", "", text)


def wdl_tasks(src):
    """(name, body) for every top-level task in ``src``."""

    out = []
    for match in TASK_RE.finditer(src):
        after = src[match.end() :]
        nxt = BOUNDARY_RE.search(after)
        end = match.end() + (nxt.start() if nxt else len(after))
        out.append((match.group(1), src[match.start() : end]))
    return out


def _declared_names(block):
    return [m.group(1) or m.group(2) for m in DECL_RE.finditer(_strip_comments(block))]


def _placeholder_idents(text):
    """Identifiers appearing inside `~{...}` placeholders."""

    names = set()
    for expr in PLACEHOLDER_RE.findall(text):
        names.update(IDENT_RE.findall(expr))
    return names


def unreached_task_inputs(src):
    """Declared inputs of pipeline-invoking tasks that no command line can see.

    Reachability is TRANSITIVE through the task's own declarations, because that
    is how the real tasks consume half their inputs: ``LRAA_runner_task`` turns
    ``no_EM`` into ``no_EM_flag`` and ``sample_id`` into ``output_prefix_use``
    before the command interpolates either. A direct-reference-only rule would
    report those and the allowlist would have to swallow them, which is exactly
    how a guard stops guarding.

    The ``runtime`` and ``output`` blocks deliberately do NOT count: an input
    read only there never becomes a flag, which is the whole question.
    """

    offenders = {}
    for name, body in wdl_tasks(src):
        cmd_match = COMMAND_RE.search(body)
        if not cmd_match or not INVOCATION_RE.search(cmd_match.group(1)):
            continue
        input_match = INPUT_BLOCK_RE.search(body)
        if not input_match:
            continue
        declared = _declared_names(input_match.group(2))

        # Task-level declarations: everything in the body that is not an input,
        # runtime or output block. `Int x = ...` here is an intermediate the
        # command may interpolate instead of the raw input.
        rest = body[: input_match.start()] + body[input_match.end() :]
        for pat in (RUNTIME_RE, OUTPUT_RE):
            rest = pat.sub("", rest)
        rest = _strip_comments(rest)
        intermediates = {}
        for line in rest.splitlines():
            decl = DECL_RE.match(line)
            if decl and "=" in line:
                lhs = decl.group(1) or decl.group(2)
                intermediates[lhs] = set(IDENT_RE.findall(line.split("=", 1)[1]))

        reached = _placeholder_idents(cmd_match.group(1))
        # Transitive closure: an intermediate the command uses pulls in whatever
        # its own right-hand side referenced.
        changed = True
        while changed:
            changed = False
            for lhs, refs in intermediates.items():
                if lhs in reached and not refs <= reached:
                    reached |= refs
                    changed = True

        missing = [
            d for d in declared if d not in reached and d not in RUNTIME_ONLY_INPUTS
        ]
        if missing:
            offenders[name] = missing
    return offenders


# ---------------------------------------------------------------------------
# Surface B -- LRAA's flags vs. what it forwards into ChunkedRun's namespace
# ---------------------------------------------------------------------------


def argparse_dests(src):
    """dest -> first long option, for every ``add_argument`` in ``src``.

    Read from the AST rather than by regex so a multi-line call, a group
    (``primary.add_argument``) and an explicit ``dest=`` are all handled the same
    way -- and so a flag added inside an ``if`` or a loop is still seen.
    """

    dests = {}
    for node in ast.walk(ast.parse(src)):
        if not isinstance(node, ast.Call):
            continue
        func = node.func
        if not (isinstance(func, ast.Attribute) and func.attr == "add_argument"):
            continue
        longs = [
            a.value
            for a in node.args
            if isinstance(a, ast.Constant)
            and isinstance(a.value, str)
            and a.value.startswith("--")
        ]
        dest = None
        for kw in node.keywords:
            if kw.arg == "dest" and isinstance(kw.value, ast.Constant):
                dest = kw.value.value
        if dest is None:
            if not longs:
                continue
            dest = longs[0][2:].replace("-", "_")
        dests.setdefault(dest, longs[0] if longs else "--" + dest)
    return dests


def forwarded_to_chunked_run(src):
    """Names ``LRAA`` passes into ``ChunkedRun.default_args(...)``.

    ``**kwargs`` is resolved rather than ignored: the single-cell tags travel as
    ``**sc_tag_overrides``, a comprehension over a literal tuple of key names, so
    the string constants in that assignment are the keys. A dict built some way
    this cannot read yields no names, and its flags are then reported as
    unforwarded -- the safe direction to be wrong in.
    """

    tree = ast.parse(src)
    call = None
    for node in ast.walk(tree):
        if (
            isinstance(node, ast.Call)
            and isinstance(node.func, ast.Attribute)
            and node.func.attr == "default_args"
        ):
            call = node
            break
    assert call is not None, "LRAA no longer calls ChunkedRun.default_args()"

    names = {kw.arg for kw in call.keywords if kw.arg}
    star_names = [
        kw.value.id
        for kw in call.keywords
        if kw.arg is None and isinstance(kw.value, ast.Name)
    ]
    for target in star_names:
        for node in ast.walk(tree):
            if isinstance(node, ast.Assign) and any(
                isinstance(t, ast.Name) and t.id == target for t in node.targets
            ):
                names.update(
                    n.value
                    for n in ast.walk(node.value)
                    if isinstance(n, ast.Constant) and isinstance(n.value, str)
                )
    return names, call.lineno


def unforwarded_shared_flags(lraa_src, chunked_src):
    lraa = argparse_dests(lraa_src)
    chunked = argparse_dests(chunked_src)
    forwarded, _ = forwarded_to_chunked_run(lraa_src)
    shared = sorted(set(lraa) & set(chunked))
    assert len(shared) > 15, (
        "only {} shared flag(s) found -- the parsers or this extraction "
        "changed shape and the check has gone vacuous".format(len(shared))
    )
    return {d: lraa[d] for d in shared if d not in forwarded}


# ---------------------------------------------------------------------------
# Surface C -- sibling calls to one callee
# ---------------------------------------------------------------------------


def _call_bodies(src):
    """(callee, alias, set of input names) for every call in ``src``."""

    out = []
    for match in CALL_RE.finditer(src):
        i = match.end() - 1
        depth = 0
        while i < len(src):
            if src[i] == "{":
                depth += 1
            elif src[i] == "}":
                depth -= 1
                if depth == 0:
                    break
            i += 1
        body = _strip_comments(src[match.end() : i])
        names = set(re.findall(r"(?:^|,)\s*(\w+)\s*=", body, re.M))
        out.append((match.group(2).split(".")[-1], match.group(3) or match.group(2), names))
    return out


def divergent_sibling_calls(src, filename):
    """Inputs one call to a callee sets and another call to it does not."""

    by_callee = {}
    for callee, alias, names in _call_bodies(src):
        by_callee.setdefault(callee, []).append((alias, names))

    offenders = {}
    for callee, instances in by_callee.items():
        if len(instances) < 2:
            continue
        union = set().union(*[n for _, n in instances])
        allowed = set(DIVERGENT_CALL_INPUTS.get((filename, callee), {}))
        for alias, names in instances:
            missing = sorted((union - names) - allowed)
            if missing:
                offenders["{}.{}".format(callee, alias)] = missing
    return offenders


# ---------------------------------------------------------------------------
# Surface D -- negated flags vs. LRAA's _NEGATED_CONFIG_FLAGS map
# ---------------------------------------------------------------------------


def _config_key(node):
    """``KEY`` when ``node`` is ``config["KEY"]`` or ``LRAA_Globals.config["KEY"]``."""

    if not isinstance(node, ast.Subscript):
        return None
    value = node.value
    owner = value.id if isinstance(value, ast.Name) else getattr(value, "attr", None)
    if owner != "config":
        return None
    index = node.slice
    if isinstance(index, ast.Constant) and isinstance(index.value, str):
        return index.value
    return None


def _negated_arg_dest(value):
    """The args dest when ``value`` is EXACTLY the boolean negation of one of them.

    Narrow on purpose. Both spellings the driver actually uses are accepted --
    ``not args.x`` (with or without a ``bool()`` around it) and
    ``False if args.x else True`` -- and anything else is left alone. A looser
    matcher would start reporting derivations such as
    ``config["min_PolyA_iso_fraction"] = resolve(args.min_isoform_fraction, ...)``,
    which are NOT negations and are routed correctly by re-derivation in the
    worker. Reporting those would force allowlist entries, and an allowlist is
    how a guard stops guarding.
    """

    dests = {
        node.attr
        for node in ast.walk(value)
        if isinstance(node, ast.Attribute)
        and isinstance(node.value, ast.Name)
        and node.value.id == "args"
    }
    if len(dests) != 1:
        return None
    dest = next(iter(dests))
    if isinstance(value, ast.UnaryOp) and isinstance(value.op, ast.Not):
        return dest
    if (
        isinstance(value, ast.IfExp)
        and isinstance(value.body, ast.Constant)
        and isinstance(value.orelse, ast.Constant)
        and value.body.value is False
        and value.orelse.value is True
    ):
        return dest
    return None


def negated_config_assignments(src):
    """dest -> config key, for every ``config[KEY] = not args.DEST`` in ``src``."""

    found = {}
    for node in ast.walk(ast.parse(src)):
        if not isinstance(node, ast.Assign):
            continue
        for target in node.targets:
            key = _config_key(target)
            if key is None:
                continue
            dest = _negated_arg_dest(node.value)
            if dest is not None:
                found.setdefault(dest, key)
    return found


def declared_negated_map(src):
    """LRAA's ``_NEGATED_CONFIG_FLAGS`` literal, read without importing LRAA."""

    for node in ast.walk(ast.parse(src)):
        if isinstance(node, ast.Assign) and any(
            isinstance(t, ast.Name) and t.id == "_NEGATED_CONFIG_FLAGS"
            for t in node.targets
        ):
            return {
                k.value: v.value
                for k, v in zip(node.value.keys, node.value.values)
                if isinstance(k, ast.Constant) and isinstance(v, ast.Constant)
            }
    raise AssertionError("LRAA no longer defines _NEGATED_CONFIG_FLAGS")


def unmapped_negated_flags(src):
    found = negated_config_assignments(src)
    assert found, "no `config[...] = not args....` assignments found at all"
    declared = declared_negated_map(src)
    return {
        dest: key
        for dest, key in found.items()
        if declared.get(dest) != key
    }


# ---------------------------------------------------------------------------
# Tests
# ---------------------------------------------------------------------------


def _wdl_files():
    files = sorted(WDL_DIR.rglob("*.wdl"))
    assert files, "no WDL files found under {}".format(WDL_DIR)
    return files


@pytest.mark.parametrize("path", _wdl_files(), ids=lambda p: p.name)
def test_pipeline_task_inputs_reach_the_command_line(path):
    """Surface A. A task input no command line sees is an input that does nothing."""

    offenders = unreached_task_inputs(path.read_text())
    assert not offenders, (
        "{} declares task input(s) that never reach the LRAA/ChunkedRun command "
        "line. The workflow above passes them, the run completes, and the "
        "setting is simply ignored -- which is how `String? oversimplify` in "
        "process_chunk shipped. Emit the flag, or add the input to "
        "RUNTIME_ONLY_INPUTS with a reason if it is genuinely consumed "
        "elsewhere:\n  {}".format(
            path.name,
            "\n  ".join("{}: {}".format(t, ", ".join(v)) for t, v in offenders.items()),
        )
    )


def test_shared_cli_flags_are_forwarded_to_chunked_run():
    """Surface B. LRAA's default_args enumeration is an allowlist, so an omission
    is silent: the flag is parsed off the driver's command line, dropped, and the
    chunked arm runs on ChunkedRun's own default instead."""

    offenders = unforwarded_shared_flags(
        LRAA_DRIVER.read_text(), CHUNKED_RUN.read_text()
    )
    assert not offenders, (
        "LRAA declares these flag(s), pylib/ChunkedRun.py accepts the same "
        "destination, and LRAA's ChunkedRun.default_args(...) call does not "
        "forward them -- so they are parsed and silently dropped in chunked "
        "mode, which subwdls/LRAA_runner.wdl makes EVERY WDL path:\n  {}".format(
            "\n  ".join(
                "{} ({})".format(dest, flag) for dest, flag in sorted(offenders.items())
            )
        )
    )


@pytest.mark.parametrize("path", _wdl_files(), ids=lambda p: p.name)
def test_sibling_calls_to_one_callee_do_not_diverge(path):
    """Surface C. Two calls to one subworkflow that set different inputs are
    two different pipelines, and the difference is invisible in the outputs."""

    offenders = divergent_sibling_calls(path.read_text(), path.name)
    assert not offenders, (
        "{} calls one callee more than once and the calls disagree about which "
        "inputs to set. Either pass it, or record it in "
        "DIVERGENT_CALL_INPUTS with the validation rule that forbids it:\n"
        "  {}".format(
            path.name,
            "\n  ".join(
                "{} does not set: {}".format(k, ", ".join(v))
                for k, v in sorted(offenders.items())
            ),
        )
    )


def test_negated_config_flags_are_mapped_for_chunk_workers():
    """Surface D. `_collect_explicit_config_overrides` routes a changed flag by
    `dest in config`, and a negated flag's dest is NEVER a config key -- it is the
    negation of one. Unmapped, the flag falls through every branch and no chunk
    worker hears about it, so a chunked run accepts the flag and ignores it."""

    offenders = unmapped_negated_flags(LRAA_DRIVER.read_text())
    assert not offenders, (
        "LRAA sets these config key(s) from the NEGATION of a flag and does not "
        "map the pair in _NEGATED_CONFIG_FLAGS, so the setting is dropped on "
        "every chunked run -- which is EVERY WDL path:\n  {}".format(
            "\n  ".join(
                'args.{} -> config["{}"]'.format(dest, key)
                for dest, key in sorted(offenders.items())
            )
        )
    )


def test_the_pipeline_task_set_is_what_we_think_it_is():
    """Guard the guard for surface A: pin WHICH tasks it inspects.

    Surface A only looks at tasks that run the pipeline, so a regression in that
    detection would make it inspect nothing and pass everywhere. This names the
    four tasks that must be in scope; a new one is a deliberate edit here.
    """

    found = set()
    for path in _wdl_files():
        for name, body in wdl_tasks(path.read_text()):
            cmd = COMMAND_RE.search(body)
            if cmd and INVOCATION_RE.search(cmd.group(1)):
                found.add("{}::{}".format(path.name, name))

    assert found == {
        "LRAA_chunk_scatter.wdl::make_chunks",
        "LRAA_chunk_scatter.wdl::process_chunk",
        "LRAA_runner.wdl::LRAA_runner_task",
        "LRAA_quant_by_cluster.wdl::emit_shared_chunk_plan",
    }, sorted(found)


# ---------------------------------------------------------------------------
# Negative controls. Each reverts one shipped fix in a COPY of the real source
# and asserts the check reports it. Without these the checks could quietly stop
# matching anything and every file would pass vacuously -- which is how all the
# original defects survived a green suite.
# ---------------------------------------------------------------------------


def test_surface_a_catches_the_dropped_oversimplify_flag():
    """F1: delete process_chunk's `--oversimplify` emission and it must be found."""

    src = (WDL_DIR / "subwdls" / "LRAA_chunk_scatter.wdl").read_text()
    reverted, n = re.subn(
        r"^.*--oversimplify.*\n", "", src, count=1, flags=re.M
    )
    assert n == 1, "process_chunk no longer emits --oversimplify at all"

    assert not unreached_task_inputs(src), "the real file must be clean"
    offenders = unreached_task_inputs(reverted)
    assert offenders.get("process_chunk") == ["oversimplify"], offenders


def test_surface_b_catches_the_unforwarded_oversimplify_kwarg():
    """F2: delete `oversimplify=args.oversimplify` from the forwarding block."""

    src = LRAA_DRIVER.read_text()
    reverted, n = re.subn(
        r"^\s*oversimplify=args\.oversimplify,\n", "", src, count=1, flags=re.M
    )
    assert n == 1, "LRAA no longer forwards oversimplify=args.oversimplify"

    chunked = CHUNKED_RUN.read_text()
    assert not unforwarded_shared_flags(src, chunked), "the real driver must be clean"
    offenders = unforwarded_shared_flags(reverted, chunked)
    assert sorted(offenders) == ["oversimplify"], offenders


def test_surface_c_catches_the_dropped_contig_restriction():
    """F3: delete `contig = direct_contig` from the LRAA_direct call."""

    src = (WDL_DIR / "LRAA.wdl").read_text()
    reverted, n = re.subn(
        r"^\s*contig = direct_contig,\n", "", src, count=1, flags=re.M
    )
    assert n == 1, "LRAA_direct no longer sets contig"

    assert not divergent_sibling_calls(src, "LRAA.wdl"), "the real file must be clean"
    offenders = divergent_sibling_calls(reverted, "LRAA.wdl")
    assert offenders == {"LRAA_runner.LRAA_direct": ["contig"]}, offenders


def test_surface_d_catches_the_unmapped_3prime_agreement_flag():
    """F4: drop the `no_weight_reads_by_3prime_agreement` entry from the map."""

    src = LRAA_DRIVER.read_text()
    reverted, n = re.subn(
        r'^\s*"no_weight_reads_by_3prime_agreement": '
        r'"weight_reads_by_3prime_agreement",\n',
        "",
        src,
        count=1,
        flags=re.M,
    )
    assert n == 1, "_NEGATED_CONFIG_FLAGS no longer maps the 3' agreement flag"

    assert not unmapped_negated_flags(src), "the real driver must be clean"
    offenders = unmapped_negated_flags(reverted)
    assert offenders == {
        "no_weight_reads_by_3prime_agreement": "weight_reads_by_3prime_agreement"
    }, offenders


def test_surface_a_catches_the_unemitted_3prime_agreement_input():
    """F4, WDL side, on BOTH leaf shapes.

    ``process_chunk`` reaches the setting through a config-override file rather
    than a flag, because ChunkedRun routes config that way on purpose -- so the
    fix is two command lines, the conditional write and the
    ``--worker_config_json`` that names it. Reverting both is reverting the fix.
    """

    src = (WDL_DIR / "subwdls" / "LRAA_chunk_scatter.wdl").read_text()
    reverted, n = re.subn(
        r'^.*if \[\[ "~\{no_weight_reads_by_3prime_agreement\}" == "true" \]\]; then\n',
        "    if false; then\n",
        src,
        count=1,
        flags=re.M,
    )
    assert n == 1, "process_chunk no longer branches on the 3' agreement input"
    reverted, n = re.subn(
        r"^.*--worker_config_json.*\n", "", reverted, count=1, flags=re.M
    )
    assert n == 1, "process_chunk no longer passes --worker_config_json"

    offenders = unreached_task_inputs(reverted)
    assert offenders.get("process_chunk") == [
        "no_weight_reads_by_3prime_agreement"
    ], offenders

    runner = (WDL_DIR / "subwdls" / "LRAA_runner.wdl").read_text()
    reverted, n = re.subn(
        r"^.*--no_weight_reads_by_3prime_agreement.*\n", "", runner, count=1, flags=re.M
    )
    assert n == 1, "LRAA_runner_task no longer emits the flag"
    offenders = unreached_task_inputs(reverted)
    assert offenders.get("LRAA_runner_task") == [
        "no_weight_reads_by_3prime_agreement"
    ], offenders
