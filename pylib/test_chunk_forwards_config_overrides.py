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
