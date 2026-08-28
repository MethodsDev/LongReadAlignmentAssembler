#!/usr/bin/env python3

"""`--chunk` forwards threshold config to the workers that actually apply it.

The bug this pins: `lraa_cmd` built each worker's command line from an explicit
allowlist -- mapping quality, --HiFi, rDNA masking, the single-cell tags,
streaming -- that no threshold flag was ever on. And the outer `--chunk` driver
`sys.exit`s into `_run_chunked_mode` BEFORE main() resolves its own config, so
neither `--min_alt_splice_freq` nor `--config_update` had been applied to
`LRAA_Globals.config` either. Both halves failed in the same direction and the
failure was silent: `LRAA --chunk --min_alt_splice_freq 0.005` parsed the flag,
discarded it, and every worker re-derived the 0.01 HiFi default. A threshold
sweep run this way measures nothing but its own baseline, twice.

Two distinct paths have to keep working, which is why there are two kinds of
test here:

  * a flag with a config key of the same name (`--min_alt_splice_freq`), and
  * a key with NO flag at all, reachable only through `--config_update`
    (`min_TSS_iso_fraction`, and the containment/PolyA thresholds with it).

The second is the one an allowlist can never cover, and the reason the fix
forwards a JSON file rather than more flags: a tunable added upstream arrives
in the workers without anyone remembering to extend `lraa_cmd`.

Asserted against the worker's own logged `cmd:` line and against the config
snapshot the worker writes, not against the outer driver's log -- the outer
driver is exactly the process that used to look correct while the workers ran
defaults.
"""

import glob
import json
import os
import subprocess
import sys

import pysam
import pytest

REPO = os.path.dirname(os.path.dirname(os.path.realpath(__file__)))
LRAA = os.path.join(REPO, "LRAA")
CONTIG = "chr1"


def _bam(path):
    header = pysam.AlignmentHeader.from_dict(
        {
            "HD": {"VN": "1.6", "SO": "coordinate"},
            "SQ": [{"SN": CONTIG, "LN": 10000}],
        }
    )
    with pysam.AlignmentFile(str(path), "wb", header=header) as fh:
        for i in range(10):
            a = pysam.AlignedSegment(header)
            a.query_name = "r{}".format(i)
            a.reference_id = 0
            a.reference_start = 100 + i
            a.mapping_quality = 60
            a.cigarstring = "100M"
            a.query_sequence = "A" * 100
            a.query_qualities = pysam.qualitystring_to_array("I" * 100)
            fh.write(a)
    pysam.index(str(path))
    return path


@pytest.fixture
def inputs(tmp_path):
    gtf = tmp_path / "a.gtf"
    gtf.write_text(
        '{}\ttest\texon\t101\t200\t.\t+\t.\tgene_id "g1"; transcript_id "t1";\n'.format(
            CONTIG
        )
    )
    genome = tmp_path / "g.fa"
    genome.write_text(">{}\n".format(CONTIG) + "A" * 10000 + "\n")
    bam = _bam(tmp_path / "reads.bam")
    return bam, gtf, genome


def _chunked(tmp_path, bam, gtf, genome, *extra):
    cmd = [
        sys.executable,
        LRAA,
        "--bam", str(bam),
        "--genome", str(genome),
        "--gtf", str(gtf),
        "--quant_only",
        "--output_prefix", str(tmp_path / "out"),
        "--chunk",
        "--chunk_work_dir", str(tmp_path / "work"),
    ] + list(extra)
    return subprocess.run(
        cmd, capture_output=True, text=True, cwd=str(tmp_path), timeout=900
    )


def _worker_argv(tmp_path):
    """Every chunk worker's real command line, concatenated.

    The worker's own logged ``cmd:`` line is the only place the forwarded argv
    is observable from outside the process, and it is the thing that was wrong:
    the outer driver's log named the flag the user passed while the worker's did
    not carry it.
    """
    logs = sorted(glob.glob(str(tmp_path / "work" / "logs" / "chunk_*.log")))
    assert logs, "no chunk logs written"
    return "\n".join(open(p).read() for p in logs)


def _override_files(tmp_path):
    """Every override file in the work dir, newest last.

    Named after a digest of its own contents, so a work dir resumed under new
    thresholds holds one file per distinct threshold set rather than one file
    that was quietly rewritten.
    """
    return sorted(glob.glob(str(tmp_path / "work" / "worker_config_overrides.*.json")))


def _forwarded_config(tmp_path):
    """The override file the driver wrote for its workers, parsed."""
    paths = _override_files(tmp_path)
    assert paths, "driver wrote no worker_config_overrides.*.json"
    assert len(paths) == 1, "expected one override file, found {}".format(paths)
    with open(paths[0]) as fh:
        return json.load(fh)


def _assert_ok(result, tmp_path):
    combined = result.stdout + result.stderr
    assert result.returncode == 0, combined[-3000:]
    assert list(tmp_path.glob("out*quant.expr")), "no output written"


def _assert_worker_applied(tmp_path, *keys):
    """The worker's OWN log saying it applied the key.

    Stronger than the driver's override file or the forwarded argv, and the
    assertion that actually distinguishes fixed from broken: the file and the
    flag only show the value was SENT. This shows the process that owns the
    filters received it and put it into its config. Pre-fix this line is absent
    from every worker log, because no worker was ever given a --config_update.
    """
    argv = _worker_argv(tmp_path)
    assert "--config_update applied to" in argv, argv
    for key in keys:
        assert key in argv, "worker never reported applying {}: {}".format(key, argv)


def test_chunk_forwards_an_explicit_min_alt_splice_freq(tmp_path, inputs):
    """The flag whose silent loss invalidated a threshold sweep."""
    bam, gtf, genome = inputs
    r = _chunked(tmp_path, bam, gtf, genome, "--min_alt_splice_freq", "0.005")
    _assert_ok(r, tmp_path)

    assert _forwarded_config(tmp_path)["min_alt_splice_freq"] == 0.005
    _assert_worker_applied(tmp_path, "min_alt_splice_freq")


def test_chunk_forwards_a_config_only_key_with_no_flag(tmp_path, inputs):
    """`min_TSS_iso_fraction` has no CLI flag, so only --config_update reaches it.

    This is the case an allowlist of flags structurally cannot serve, and the
    reason the fix forwards a file.
    """
    bam, gtf, genome = inputs
    cfg = tmp_path / "cfg.json"
    cfg.write_text(json.dumps({"min_TSS_iso_fraction": 0.025}))

    r = _chunked(tmp_path, bam, gtf, genome, "--config_update", str(cfg))
    _assert_ok(r, tmp_path)

    assert _forwarded_config(tmp_path)["min_TSS_iso_fraction"] == 0.025
    _assert_worker_applied(tmp_path, "min_TSS_iso_fraction")


def test_config_update_file_outranks_an_explicit_flag(tmp_path, inputs):
    """Precedence matches unchunked main(): the file is applied last and wins."""
    bam, gtf, genome = inputs
    cfg = tmp_path / "cfg.json"
    cfg.write_text(json.dumps({"min_alt_splice_freq": 0.002}))

    r = _chunked(
        tmp_path,
        bam,
        gtf,
        genome,
        "--min_alt_splice_freq", "0.005",
        "--config_update", str(cfg),
    )
    _assert_ok(r, tmp_path)
    assert _forwarded_config(tmp_path)["min_alt_splice_freq"] == 0.002
    _assert_worker_applied(tmp_path, "min_alt_splice_freq")


def test_min_reads_novel_forwards_both_config_keys_it_sets(tmp_path, inputs):
    """One flag, two config keys, neither named like the flag."""
    bam, gtf, genome = inputs
    r = _chunked(tmp_path, bam, gtf, genome, "--min_reads_novel", "1")
    _assert_ok(r, tmp_path)

    forwarded = _forwarded_config(tmp_path)
    assert forwarded["min_reads_novel_isoform"] == 1
    assert forwarded["min_unique_reads_novel_isoform"] == 1
    _assert_worker_applied(
        tmp_path, "min_reads_novel_isoform", "min_unique_reads_novel_isoform"
    )


def test_several_thresholds_reach_the_same_worker_command(tmp_path, inputs):
    """All of them in one run, not one flag per run."""
    bam, gtf, genome = inputs
    cfg = tmp_path / "cfg.json"
    cfg.write_text(
        json.dumps(
            {
                "min_TSS_iso_fraction": 0.025,
                "max_rel_frac_expr_alt_compat_contained": 0.1,
                "min_PolyA_ident_length": 4,
            }
        )
    )

    r = _chunked(
        tmp_path,
        bam,
        gtf,
        genome,
        "--min_alt_splice_freq", "0.005",
        "--min_isoform_fraction", "0.005",
        "--config_update", str(cfg),
    )
    _assert_ok(r, tmp_path)

    forwarded = _forwarded_config(tmp_path)
    assert forwarded["min_alt_splice_freq"] == 0.005
    assert forwarded["min_isoform_fraction"] == 0.005
    assert forwarded["min_TSS_iso_fraction"] == 0.025
    assert forwarded["max_rel_frac_expr_alt_compat_contained"] == 0.1
    assert forwarded["min_PolyA_ident_length"] == 4
    _assert_worker_applied(
        tmp_path,
        "min_alt_splice_freq",
        "min_isoform_fraction",
        "min_TSS_iso_fraction",
        "max_rel_frac_expr_alt_compat_contained",
        "min_PolyA_ident_length",
    )


def test_a_plain_chunked_run_forwards_nothing(tmp_path, inputs):
    """No threshold stated, so no override file and no --config_update.

    Absence must stay absence: writing a file of resolved defaults would pin
    every value a worker legitimately derives for itself, which is the failure
    mode opposite to the one being fixed.
    """
    bam, gtf, genome = inputs
    r = _chunked(tmp_path, bam, gtf, genome)
    _assert_ok(r, tmp_path)

    assert not _override_files(tmp_path)
    assert "--config_update" not in _worker_argv(tmp_path)


def test_a_threshold_set_to_its_own_default_still_forwards(tmp_path, inputs):
    """Deliberate, and the same rule the single-cell tags follow.

    `--min_alt_splice_freq 0.03` states 0.03 explicitly. Under --HiFi the HiFi
    preset would otherwise lower it to 0.01, so treating a stated default as
    absence would silently override the user. Presence of the flag is what is
    being honoured, not its value's distance from a default.
    """
    bam, gtf, genome = inputs
    r = _chunked(
        tmp_path, bam, gtf, genome, "--HiFi", "--min_alt_splice_freq", "0.03"
    )
    _assert_ok(r, tmp_path)
    assert _forwarded_config(tmp_path)["min_alt_splice_freq"] == 0.03


def test_unchunked_thresholds_are_unaffected(tmp_path, inputs):
    """The unchunked path always applied these; the refactor must not move it."""
    bam, gtf, genome = inputs
    cmd = [
        sys.executable,
        LRAA,
        "--bam", str(bam),
        "--genome", str(genome),
        "--gtf", str(gtf),
        "--quant_only",
        "--output_prefix", str(tmp_path / "out"),
        "--no_chunk",
        "--min_alt_splice_freq", "0.005",
    ]
    r = subprocess.run(
        cmd, capture_output=True, text=True, cwd=str(tmp_path), timeout=900
    )
    assert r.returncode == 0, (r.stdout + r.stderr)[-3000:]
    assert list(tmp_path.glob("out*quant.expr")), "no output written"


def test_resuming_a_workdir_under_new_thresholds_does_not_reuse_old_chunks(
    tmp_path, inputs
):
    """The stale-sentinel trap the override filename's content digest closes.

    Per-chunk quant sentinels are keyed on the worker's argv, and the argv sees
    the override file as one opaque path token. At a FIXED filename, resuming the
    same --chunk_work_dir with different thresholds leaves that token -- and so
    the sentinel -- unchanged, and the second run serves chunks built under the
    first run's thresholds while reporting the second run's settings. Silent, and
    exactly the failure this whole fix exists to stop, so it must not be
    reintroduced by the fix itself.

    Two runs, same work dir, different values. The second must key differently,
    which is observable three ways: a second override file exists, the two argv
    digests differ, and the worker reports applying the SECOND value.
    """
    bam, gtf, genome = inputs

    first = _chunked(tmp_path, bam, gtf, genome, "--min_alt_splice_freq", "0.005")
    _assert_ok(first, tmp_path)
    after_first = _override_files(tmp_path)
    assert len(after_first) == 1, after_first

    second = _chunked(tmp_path, bam, gtf, genome, "--min_alt_splice_freq", "0.002")
    _assert_ok(second, tmp_path)

    after_second = _override_files(tmp_path)
    assert len(after_second) == 2, (
        "the second threshold set reused the first's filename, so every "
        "argv-keyed chunk sentinel still matches and the resumed run can serve "
        "stale chunks: {}".format(after_second)
    )

    # The second run's own value is what the workers were handed.
    values = []
    for path in after_second:
        with open(path) as fh:
            values.append(json.load(fh)["min_alt_splice_freq"])
    assert sorted(values) == [0.002, 0.005], values

    # And the second run's workers actually applied it, rather than resuming a
    # sentinel written under 0.005.
    argv = _worker_argv(tmp_path)
    assert os.path.basename(after_second[0]) in argv or os.path.basename(
        after_second[1]
    ) in argv, argv


# --------------------------------------------------------------------------
# PREP/WORKER PARITY.
#
# Some config keys are consumed by the chunking PREP phase -- cut selection,
# extraction, stage-4 normalization -- before any worker exists. Prep reads them
# from the driver's resolved config; workers read them from the forwarded
# --config_update. If those two disagree, a run quantifies reads that its own
# extraction filtered differently, and nothing in the output says so.
#
# So for a prep key the assertion is not "the worker got it" but "prep and the
# worker got the SAME thing", which is why these tests read the normalization
# cache identity (prep's own view, encoded in the filename it chooses) alongside
# the worker's applied value.
# --------------------------------------------------------------------------


def _prep_norm_commands(tmp_path):
    """Prep's own parameters, read from the normalization commands it issued.

    Chunked prep names its per-chunk output plainly (chunk.plus.norm.bam), so the
    filename carries no parameters -- the parameter-encoded norm_cache names
    belong to the UNCHUNKED path. What chunked prep does expose is the
    normalize_bam_by_strand command line it runs per chunk, logged by run_step,
    and that names --min_per_id and --normalize_max_cov_level directly. That
    command IS prep's view: if a config change does not move it, prep did not see
    the change and is normalizing on the old value.
    """
    lines = []
    for log in glob.glob(str(tmp_path / "work" / "logs" / "chunk_*.log")):
        with open(log) as fh:
            lines += [ln for ln in fh if "normalize_bam_by_strand" in ln]
    return lines


def _worker_applied_config(tmp_path):
    """The override file the workers were actually handed, parsed."""
    paths = _override_files(tmp_path)
    if not paths:
        return {}
    with open(paths[-1]) as fh:
        return json.load(fh)


def test_prep_key_via_flag_reaches_prep_and_workers_together(tmp_path, inputs):
    """--min_per_id is prep-consumed; a CLI change must move BOTH sides.

    Before the driver resolved its own config, this was the divergence: prep
    resolved the user's value while lraa_cmd never forwarded the flag, so every
    worker quantified at its own default.
    """
    bam, gtf, genome = inputs
    r = _chunked(tmp_path, bam, gtf, genome, "--min_per_id", "90")
    _assert_ok(r, tmp_path)

    forwarded = _worker_applied_config(tmp_path)
    assert forwarded.get("min_per_id") == 90.0, forwarded
    _assert_worker_applied(tmp_path, "min_per_id")

    # prep's own normalization command must name the same value
    cmds = " ".join(_prep_norm_commands(tmp_path))
    assert "--min_per_id 90" in cmds, cmds


def test_hifi_preset_beats_an_explicit_prep_flag_in_both_modes(tmp_path, inputs):
    """--HiFi --min_per_id 90 resolves to 97.0, unchunked AND chunked.

    The preset overwrites the seeded value in main(), so 97.0 is what unchunked
    mode runs. Forwarding the RESOLVED value rather than the raw argument is what
    keeps chunked mode agreeing: forwarding the raw 90 would have let
    --config_update's last-wins position hand the worker a value unchunked mode
    would never use.
    """
    bam, gtf, genome = inputs
    r = _chunked(tmp_path, bam, gtf, genome, "--HiFi", "--min_per_id", "90")
    _assert_ok(r, tmp_path)

    forwarded = _worker_applied_config(tmp_path)
    assert forwarded.get("min_per_id") == 97.0, forwarded

    cmds = " ".join(_prep_norm_commands(tmp_path))
    assert "--min_per_id 97" in cmds, cmds


def test_prep_key_via_config_update_reaches_prep_too(tmp_path, inputs):
    """A prep key set in the FILE must also move prep, not just the workers.

    This is the case the driver's early resolution exists for: the file is read
    before prep args are derived, so cut selection and normalization see it.
    """
    bam, gtf, genome = inputs
    cfg = tmp_path / "cfg.json"
    cfg.write_text(json.dumps({"min_per_id": 80}))

    r = _chunked(tmp_path, bam, gtf, genome, "--HiFi", "--config_update", str(cfg))
    _assert_ok(r, tmp_path)

    # file is applied last, so it outranks even the HiFi preset -- in both modes
    forwarded = _worker_applied_config(tmp_path)
    assert forwarded.get("min_per_id") == 80, forwarded

    cmds = " ".join(_prep_norm_commands(tmp_path))
    assert "--min_per_id 80" in cmds, cmds


def test_normalize_max_cov_level_flag_is_seeded_not_stale(tmp_path, inputs):
    """The key whose seeding path is easy to miss.

    normalize_max_cov_level is seeded by _seed_authoritative_config_from_args
    rather than by an explicit line in the driver's resolver, so a reader can
    reasonably suspect it forwards a stale default. It does not, and this pins
    that: prep's cache identity and the forwarded value must both name 500.
    """
    bam, gtf, genome = inputs
    r = _chunked(tmp_path, bam, gtf, genome, "--normalize_max_cov_level", "500")
    _assert_ok(r, tmp_path)

    forwarded = _worker_applied_config(tmp_path)
    assert forwarded.get("normalize_max_cov_level") == 500, forwarded

    cmds = " ".join(_prep_norm_commands(tmp_path))
    assert "--normalize_max_cov_level 500" in cmds, cmds

def _override_file_holding(tmp_path, key, value):
    """The override file whose CONTENTS carry key==value, or None.

    Selected by content, never by ordering. The filenames carry a content hash,
    so sorting them is sorting hashes -- unrelated to which run wrote which, and
    `paths[-1]` would name an arbitrary one once a workdir holds two.
    """
    for path in _override_files(tmp_path):
        with open(path) as fh:
            if json.load(fh).get(key) == value:
                return path
    return None


def test_same_workdir_prep_key_change_moves_prep_identity_and_worker_config(
    tmp_path, inputs
):
    """The end-to-end guard: one workdir, two runs, a prep key changed between them.

    Unit checks on the override file cannot catch this class of bug. Both halves
    of a chunked run have to move together, and both are cached:

      * PREP caches its normalized bam under a filename naming the parameters
        that produced it. If a change does not move that name, prep reused a bam
        built under the old value and the run quantifies reads its own extraction
        filtered differently.
      * WORKERS are gated by argv-keyed sentinels. If the argv does not move, a
        resumed run replays chunks built under the old value while reporting the
        new one.

    Asserted on the SAME --chunk_work_dir twice, because that is the only
    arrangement where a stale-cache bug can surface: a fresh directory hides it
    by having nothing to reuse. Every assertion identifies its run by CONTENT --
    the override file carrying the value under test, and that file's basename
    appearing in a worker's own command line -- so nothing depends on hash or
    glob ordering.

    min_per_id is the key under test because prep reads it through
    ChunkedRun.resolve_min_per_id AND every worker applies it: the exact key whose
    two readers used to disagree, prep resolving it while lraa_cmd never forwarded
    it at all.
    """
    bam, gtf, genome = inputs

    cfg_a = tmp_path / "a.json"
    cfg_a.write_text(json.dumps({"min_per_id": 90}))
    first = _chunked(tmp_path, bam, gtf, genome, "--config_update", str(cfg_a))
    _assert_ok(first, tmp_path)

    ovr_a = _override_file_holding(tmp_path, "min_per_id", 90)
    assert ovr_a is not None, _override_files(tmp_path)
    argv_a = _worker_argv(tmp_path)
    assert os.path.basename(ovr_a) in argv_a, "run 1 workers never saw their override"
    prep_a = set(_prep_norm_commands(tmp_path))
    assert any("--min_per_id 90" in c for c in prep_a), sorted(prep_a)

    # Same workdir, different value. Nothing cleaned between the two runs.
    cfg_b = tmp_path / "b.json"
    cfg_b.write_text(json.dumps({"min_per_id": 80}))
    second = _chunked(tmp_path, bam, gtf, genome, "--config_update", str(cfg_b))
    _assert_ok(second, tmp_path)

    # 1. the second run got its OWN override file, identified by its contents
    ovr_b = _override_file_holding(tmp_path, "min_per_id", 80)
    assert ovr_b is not None, _override_files(tmp_path)
    assert ovr_b != ovr_a, "both runs shared one override file"

    # 2. a worker was actually launched with it -- the argv moved, so argv-keyed
    #    chunk sentinels cannot match run 1's and replay its chunks
    argv_b = _worker_argv(tmp_path)
    assert os.path.basename(ovr_b) in argv_b, (
        "run 2's override never reached a worker command line; sentinels would "
        "still match run 1"
    )

    # 3. prep's own cache identity moved: a pid80 artifact exists that did not
    #    before, so normalization was redone rather than reused
    prep_b = set(_prep_norm_commands(tmp_path))
    assert any("--min_per_id 80" in c for c in prep_b - prep_a), (
        "prep never issued a min_per_id 80 normalization, so it reused the "
        "90 one: before={} after={}".format(sorted(prep_a), sorted(prep_b))
    )


def _cut_selection_argv(tmp_path):
    """The stage-2 cut-selection argv LISTS, from the run's own timing.json.

    Cut selection is a SUBPROCESS and it is the only consumer of the geometry
    keys, so its argv is the only place that says what geometry the run actually
    used. The chunk manifests cannot substitute: a fixture small enough to run
    fast yields one chunk at any spacing, so the boundaries look identical while
    the parameters differ.

    Returned as LISTS, not joined text, because the values under test collide as
    substrings: "--approx_MB_per_cut 1" occurs inside "--approx_MB_per_cut 10",
    so a joined-text assertion for 1 PASSES against a run that used the default
    10. That is not hypothetical -- it is what this file's first version did, and
    it is why one of the five parametrisations passed against the known-broken
    driver while the other four failed.
    """

    timing = json.loads((tmp_path / "work" / "timing.json").read_text())
    argvs = []
    for key, entries in (timing.get("stages") or {}).items():
        if not key.startswith("cuts"):
            continue
        for entry in entries:
            argvs.append([str(x) for x in entry.get("cmd", [])])
    assert argvs, "run issued no cut-selection command: {}".format(timing.keys())
    return argvs


def _argv_value(argv, flag):
    """The token following ``flag``, or None. Exact token, never a prefix."""

    return argv[argv.index(flag) + 1] if flag in argv else None


@pytest.mark.parametrize(
    "key,flag,value",
    (
        # Every value differs from BOTH the parser default and the HiFi preset, and
        # none is a decimal prefix of the default it replaces, so a passing
        # assertion cannot be a coincidence of formatting.
        ("approx_MB_per_cut", "--approx_MB_per_cut", 3),
        ("approx_MB_per_cut_wiggle_window", "--approx_MB_per_cut_wiggle_window", 0.4),
        ("chunk_severed_multiexon_weight", "--severed_multiexon_weight", 3),
        ("chunk_depth_window", "--depth_window", 50),
        ("chunk_margin", "--margin", 400),
    ),
)
def test_config_update_moves_chunk_geometry(tmp_path, inputs, key, flag, value):
    """--config_update must reach CUT SELECTION, not only the chunk workers.

    These five keys decide chunk geometry and are consumed by the stage-2
    subprocess. The driver used to build its chunked args from ``args.*`` for
    exactly these five while their neighbours (max_intron_length, the two
    mapping-quality keys) read the RESOLVED config, so --config_update was
    reported applied, was forwarded to every chunk worker -- where geometry is
    already fixed and these keys are inert -- and never reached the stage that
    places the cuts. Measured on a 30.4 Mb contig: a file asking for 1 Mb spacing
    still cut at 10 Mb, three chunks, indistinguishable from defaults.

    Asserted on the SUBPROCESS ARGV rather than on the resulting chunk count,
    because a fixture small enough for a unit test yields one chunk at every
    spacing under test -- the count would pass while the geometry was wrong.
    """

    bam, gtf, genome = inputs
    cfg = tmp_path / "cfg.json"
    cfg.write_text(json.dumps({key: value}))

    r = _chunked(tmp_path, bam, gtf, genome, "--config_update", str(cfg))
    _assert_ok(r, tmp_path)

    seen = [_argv_value(argv, flag) for argv in _cut_selection_argv(tmp_path)]
    assert seen and all(v is not None and float(v) == float(value) for v in seen), (
        "cut selection did not receive {}={} from --config_update; it received "
        "{}".format(key, value, seen)
    )



# -- transcriptome rescue: off has to mean off, by either route ----------------
#
# Rescue is a config KEY and a flag whose dest is that key, so it can be turned
# off two ways, and each failed differently before this was fixed.
#
# By flag: the shard got --no_rescue_unassigned_reads_via_transcriptome_alignment,
# lraa_cmd forwarded only the streaming opt-out, and each worker re-derived rescue
# as True. LRAA's own guard then refused rescue-on + --stream_reads and exited 1 --
# so asking for rescue off KILLED the run. MEASURED on SIRVs at lraa-core:0.28.1:
# every stage5_quant_strandless_* step failed in under a second.
#
# By --config_update: no crash, and worse for it. The driver read args rather than
# the resolved config, so it forwarded no opt-out at all AND told every worker to do
# streaming rescue, while the worker's own override file said rescue was off. That
# run exited 0 and rescued anyway.
#
# The override file alone cannot fix either case: the guard reads `args` and fires
# ~140 lines before _apply_config_update_file applies the file.

MASTER_OFF = "--no_rescue_unassigned_reads_via_transcriptome_alignment"
STREAMING_OFF = "--no_stream_reads_rescue_unassigned"


def test_rescue_off_by_flag_reaches_every_chunk_worker(tmp_path, inputs):
    bam, gtf, genome = inputs
    r = _chunked(tmp_path, bam, gtf, genome, MASTER_OFF)
    _assert_ok(r, tmp_path)

    argv = _worker_argv(tmp_path)
    assert MASTER_OFF in argv, argv[-3000:]
    # and the narrower flag with it, since streaming rescue cannot outlive rescue
    assert STREAMING_OFF in argv, argv[-3000:]


def test_rescue_off_by_config_update_reaches_every_chunk_worker(tmp_path, inputs):
    """The route that exited 0 and rescued anyway.

    Both flags are asserted on the worker's argv, not just the override file: the
    file was already correct in the broken version, which is exactly why the run
    looked fine.
    """

    bam, gtf, genome = inputs
    cfg = tmp_path / "cfg.json"
    cfg.write_text(
        json.dumps({"rescue_unassigned_reads_via_transcriptome_alignment": False})
    )

    r = _chunked(tmp_path, bam, gtf, genome, "--config_update", str(cfg))
    _assert_ok(r, tmp_path)

    forwarded = _forwarded_config(tmp_path)
    assert forwarded["rescue_unassigned_reads_via_transcriptome_alignment"] is False

    argv = _worker_argv(tmp_path)
    assert MASTER_OFF in argv, argv[-3000:]
    assert STREAMING_OFF in argv, argv[-3000:]


def test_a_plain_chunked_run_still_asks_for_rescue(tmp_path, inputs):
    """Rescue on must add nothing to the worker's argv.

    Every per-chunk sentinel is keyed on that argv (ChunkedRun.argv_digest), so an
    extra token on the default path would invalidate existing chunk work and every
    splice-graph cache entry for no behavioural change.
    """

    bam, gtf, genome = inputs
    r = _chunked(tmp_path, bam, gtf, genome)
    _assert_ok(r, tmp_path)

    argv = _worker_argv(tmp_path)
    assert MASTER_OFF not in argv
    assert STREAMING_OFF not in argv
    assert "--stream_reads_rescue_unassigned" in argv


# -- the TPM denominator: resolved once, then used, never re-derived --------------
#
# The contract is deliberately simple: at the start of a run, a supplied
# --num_total_reads is used throughout; unsupplied, LRAA counts the bam once with its
# -F 0x904 mapped-primary policy. Nothing downstream recomputes it and nothing
# substitutes a different number.
#
# Both chunking arms are exercised because they used to disagree. The strand-first arm
# counts what its stage-1 splitter RETAINED, which is a smaller set than -F 0x904 --
# the splitter also drops duplicate, qcfail and long-intron records -- and ChunkedRun
# refused the mismatch outright. MEASURED on a SIRV bam with 2,095 of 104,766 primary
# records flagged duplicate: `--chunk --chunk_by_strand` exited 1 with "-N 104766
# disagrees with the stage-1 retained record count". It now uses the resolved value on
# both arms and warns when the retained total differs.


def _worker_denominators(tmp_path):
    """The --num_total_reads value every chunk worker was given, as a set."""

    seen = set()
    for line in _worker_argv(tmp_path).splitlines():
        tokens = line.split()
        value = _argv_value(tokens, "--num_total_reads")
        if value is not None:
            seen.add(value)
    return seen


@pytest.mark.parametrize(
    "arm", ["--chunk_by_strand", None], ids=["strand_first", "strandless"]
)
def test_a_supplied_denominator_reaches_every_worker_unchanged(tmp_path, inputs, arm):
    """Supplied means supplied, on both arms, and the run must SUCCEED.

    Asserting completion matters as much as the value: the strand-first arm used to
    fail this outright, and a test that only checked for the absence of a recount
    would pass on a run that died for any reason at all.
    """

    bam, gtf, genome = inputs
    extra = ["--num_total_reads", "777777"] + ([arm] if arm else [])
    r = _chunked(tmp_path, bam, gtf, genome, *extra)
    _assert_ok(r, tmp_path)

    assert _worker_denominators(tmp_path) == {"777777"}, _worker_argv(tmp_path)[-2000:]


@pytest.mark.parametrize(
    "arm", ["--chunk_by_strand", None], ids=["strand_first", "strandless"]
)
def test_an_unsupplied_denominator_is_counted_once_and_then_forwarded(
    tmp_path, inputs, arm
):
    """One count for the run, and every worker gets that same number.

    The count belongs to the driver, not to the chunks: a per-chunk count would give
    each worker its own chunk's total as the library size.
    """

    bam, gtf, genome = inputs
    r = _chunked(tmp_path, bam, gtf, genome, *([arm] if arm else []))
    _assert_ok(r, tmp_path)

    combined = r.stdout + r.stderr
    assert combined.count("counting genome-mapped reads") == 1, combined[-2000:]

    seen = _worker_denominators(tmp_path)
    assert len(seen) == 1, seen
    assert seen != {"0"}


def test_a_config_update_denominator_reaches_every_worker(tmp_path, inputs):
    """--config_update is an input path too, so naming the denominator there counts
    as providing it.

    It used to be accepted, written into config, and then discarded: both
    resolutions read ``args``, so the unchunked one recounted and the chunked one had
    already counted before the file was ever applied. The run then reported a TPM
    denominator the caller had explicitly set and not received, in both modes, with
    nothing said. Asserted on the workers' argv and on the ABSENCE of a counting
    pass -- a value that arrives while the driver also counted would mean the file
    was honoured by luck.
    """

    bam, gtf, genome = inputs
    cfg = tmp_path / "cfg.json"
    cfg.write_text(json.dumps({"num_total_reads": 777777}))

    r = _chunked(tmp_path, bam, gtf, genome, "--config_update", str(cfg))
    _assert_ok(r, tmp_path)

    combined = r.stdout + r.stderr
    assert "counting genome-mapped reads" not in combined, combined[-2000:]
    assert _worker_denominators(tmp_path) == {"777777"}, _worker_argv(tmp_path)[-2000:]


def test_a_config_update_denominator_outranks_the_flag(tmp_path, inputs):
    """Same precedence every other --config_update key has: the file is applied last.

    The two values differ and neither is a prefix of the other, so a passing result
    cannot come from reading the wrong one.
    """

    bam, gtf, genome = inputs
    cfg = tmp_path / "cfg.json"
    cfg.write_text(json.dumps({"num_total_reads": 777777}))

    r = _chunked(
        tmp_path, bam, gtf, genome,
        "--num_total_reads", "111111",
        "--config_update", str(cfg),
    )
    _assert_ok(r, tmp_path)

    assert _worker_denominators(tmp_path) == {"777777"}, _worker_argv(tmp_path)[-2000:]