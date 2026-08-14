#!/usr/bin/env python3

"""Tracking is always written gzipped, and the encoding must never change the content.

The tracking file is the one artifact that scales with library size, so compressing it is
what makes a very large run affordable. The risk is not compression itself but the six
places the path is constructed -- per-shard temp path, shard writer, merge source, merge
destination, the completeness check that decides a shard can be reused, and the cleanup
list. If any one of them disagrees about the suffix the failure is quiet: a finished shard
looks missing and is recomputed, or a temp file is left behind, or the merge silently
appends nothing.

These tests assert the two properties that matter: the flag off leaves the filename and
bytes exactly as before, and the flag on produces a file that decompresses to the same
rows.
"""

import gzip
import os
import shutil
import subprocess
import sys

import pytest

sys.path.insert(0, os.path.dirname(os.path.realpath(__file__)))

LRAA = os.path.join(os.path.dirname(os.path.realpath(__file__)), "..", "LRAA")


def _write_inputs(tmp_path):
    """A two-exon transcript with reads over it: enough to produce tracking rows."""
    genome = tmp_path / "g.fa"
    seq = "ACGT" * 500
    genome.write_text(">chr1\n" + "\n".join(seq[i:i+60] for i in range(0, len(seq), 60)) + "\n")

    gtf = tmp_path / "a.gtf"
    rows = []
    for lend, rend in ((101, 200), (301, 400)):
        rows.append(
            f'chr1\ttest\texon\t{lend}\t{rend}\t.\t+\t.\t'
            f'gene_id "g1"; transcript_id "t1";'
        )
    gtf.write_text("\n".join(rows) + "\n")
    return gtf, genome


def _write_bam(path, gtf, genome):
    import pysam

    header = {"HD": {"VN": "1.6", "SO": "coordinate"},
              "SQ": [{"SN": "chr1", "LN": 2000}]}
    with pysam.AlignmentFile(str(path), "wb", header=header) as fh:
        for i in range(40):
            a = pysam.AlignedSegment()
            a.query_name = f"read{i}"
            a.query_sequence = "A" * 200
            a.flag = 0
            a.reference_id = 0
            a.reference_start = 100
            a.mapping_quality = 60
            a.cigarstring = "100M100N100M"
            a.query_qualities = pysam.qualitystring_to_array("I" * 200)
            fh.write(a)
    pysam.index(str(path))
    return path


def _run(tmp_path, prefix, extra):
    """Run LRAA on the fixture."""
    gtf, genome = _write_inputs(tmp_path)
    bam = _write_bam(tmp_path / "r.bam", gtf, genome)
    cmd = [sys.executable, LRAA, "--quant_only",
           "--bam", str(bam), "--gtf", str(gtf), "--genome", str(genome),
           "--output_prefix", str(tmp_path / prefix),
           "--no_cleanup", "--no_parallelize_contigs"]
    cmd += extra
    r = subprocess.run(cmd, capture_output=True, text=True, cwd=str(tmp_path))
    return r


def _content_rows(path):
    """Non-comment lines, so the embedded command-line comment does not count as content."""
    opener = gzip.open if str(path).endswith(".gz") else open
    with opener(path, "rt") as fh:
        return [l for l in fh if not l.startswith("#")]


def test_tracking_is_always_written_gzipped(tmp_path):
    """No flag: the compressed name is the only name LRAA writes.

    The uncompressed name must not appear at all. Anything that expected it -- the WDL's own
    compensating `gzip` step, a pipeline glob -- has to see either the compressed file or
    nothing, never both.
    """
    r = _run(tmp_path, "out", [])
    assert r.returncode == 0, f"run failed:\n{r.stdout[-2000:]}\n{r.stderr[-2000:]}"

    gz = tmp_path / "out.LRAA.quant-only.quant.tracking.gz"
    plain = tmp_path / "out.LRAA.quant-only.quant.tracking"
    assert gz.exists(), "tracking was not written compressed"
    assert not plain.exists(), (
        "an uncompressed tracking file was also written: a consumer could read either"
    )
    with gzip.open(gz, "rt") as fh:
        assert fh.readline(), "compressed tracking is empty"


def test_compressed_tracking_decompresses_to_stable_content(tmp_path):
    """The compression must be transparent: readable, and the same across runs.

    This used to compare a plain run against a gzipped one, which is no longer possible now
    that the uncompressed form is not written at all. What survives from that test is the part
    that mattered: the bytes on disk decompress to well-formed tracking rows, and two runs over
    the same input agree, so compression is not perturbing content or ordering.
    """
    (tmp_path / "a").mkdir()
    (tmp_path / "b").mkdir()
    r1 = _run(tmp_path / "a", "out", [])
    r2 = _run(tmp_path / "b", "out", [])
    assert r1.returncode == 0, r1.stderr[-2000:]
    assert r2.returncode == 0, r2.stderr[-2000:]

    rows_a = _content_rows(tmp_path / "a" / "out.LRAA.quant-only.quant.tracking.gz")
    rows_b = _content_rows(tmp_path / "b" / "out.LRAA.quant-only.quant.tracking.gz")

    assert rows_a, "fixture produced no tracking rows; the test proves nothing"
    assert rows_a[0].startswith("gene_id\t"), (
        f"first content row is not the tracking header: {rows_a[0]!r}"
    )
    assert rows_b == rows_a, "two runs over the same input produced different tracking content"



def test_no_shard_tracking_of_either_encoding_survives_a_cleanup_run(tmp_path):
    """After a cleaning run, no shard tracking survives under either encoding.

    Named for what it can actually observe. It does NOT isolate the per-file cleanup list's
    stale-suffix entry: dropping that entry leaves this test green, because two other
    mechanisms remove the same file first -- the shard writer's own stale-sibling unlink
    (planting a sibling makes the shard ambiguous, which vetoes reuse and forces a recompute,
    and the writer clears it on the way through) and the whole-tree removal. That entry is
    defence in depth with no reachable path that exercises it alone; saying so here is better
    than a test whose name implies otherwise.

    What this does pin: the end state a user sees. --clean_parallel_tmp enables the per-file
    loop (it defaults off, which is why earlier versions of this test observed nothing at all)
    and --no_cleanup disables the tree removal, so the run exercises the narrower path and the
    assertion still holds.
    """
    gtf, genome = _write_inputs(tmp_path)
    bam = _write_bam(tmp_path / "r.bam", gtf, genome)
    # --clean_parallel_tmp turns ON the per-file cleanup list (it defaults OFF, which is why
    # earlier versions of this test could not observe it at all), and --no_cleanup turns OFF
    # the whole-tree removal that would otherwise delete the evidence regardless. Together
    # they isolate the one code path this test is named for.
    cmd = [sys.executable, LRAA, "--quant_only",
           "--bam", str(bam), "--gtf", str(gtf), "--genome", str(genome),
           "--output_prefix", str(tmp_path / "out"),
           "--no_parallelize_contigs",
           "--clean_parallel_tmp", "--no_cleanup"]
    # First pass with --no_cleanup so a shard tree exists, then plant a plain sibling beside
    # each .gz so the cleanup list's stale-suffix entry has something to remove. Without this
    # only .gz files exist and the test passes even if cleanup handles one encoding.
    seed_cmd = [c for c in cmd if c != "--clean_parallel_tmp"]
    seed = subprocess.run(seed_cmd, capture_output=True, text=True, cwd=str(tmp_path))
    assert seed.returncode == 0, seed.stderr[-2000:]
    planted = []
    for root, _dirs, files in os.walk(tmp_path):
        if "contigtmp" not in root:
            continue
        for f in files:
            if f.endswith(".quant.tracking.gz"):
                sib = os.path.join(root, f[: -len(".gz")])
                with open(sib, "wt") as fh:
                    fh.write("forged\tplain\tsibling\n")
                planted.append(sib)
    assert planted, "no shard tracking to plant beside; fixture is inert"

    r = subprocess.run(cmd, capture_output=True, text=True, cwd=str(tmp_path))
    assert r.returncode == 0, r.stderr[-2000:]

    for sib in planted:
        assert not os.path.exists(sib), (
            f"cleanup left the other encoding's shard tracking behind: {sib}"
        )

    assert (tmp_path / "out.LRAA.quant-only.quant.tracking.gz").exists(), (
        "final output must survive cleanup"
    )

    shard_tracking = [
        os.path.join(root, f)
        for root, _dirs, files in os.walk(tmp_path)
        for f in files
        if ".quant.tracking" in f and "contigtmp" in root
    ]
    assert not shard_tracking, (
        f"per-shard tracking files survived cleanup: {shard_tracking}"
    )


def test_legacy_uncompressed_shard_is_cleared_when_the_shard_is_rewritten(tmp_path):
    """An older LRAA's uncompressed shard tracking must not outlive the run that rewrites it.

    Nothing writes the uncompressed name any more, so this state can only arrive by upgrading
    over a shard tree an older version left -- which means seeding it rather than producing it,
    since the flag that used to create it is gone.

    It matters because the reuse check looks for the name THIS version writes. A leftover under
    the old name is invisible to that check but visible to the ambiguity guard, so leaving it
    behind means every later resume refuses the shard and recomputes it forever.
    """
    gtf, genome = _write_inputs(tmp_path)
    bam = _write_bam(tmp_path / "r.bam", gtf, genome)
    base = [sys.executable, LRAA, "--quant_only",
            "--bam", str(bam), "--gtf", str(gtf), "--genome", str(genome),
            "--output_prefix", str(tmp_path / "out"),
            "--no_cleanup", "--no_parallelize_contigs"]

    r1 = subprocess.run(base, capture_output=True, text=True, cwd=str(tmp_path))
    assert r1.returncode == 0, r1.stderr[-2000:]

    gz_shards = [os.path.join(root, f)
                 for root, _dirs, files in os.walk(tmp_path)
                 for f in files if f.endswith(".quant.tracking.gz") and "contigtmp" in root]
    assert gz_shards, "first run left no shard tracking; fixture is inert"

    # seed the legacy name beside each shard, as an older version would have left it
    legacy = []
    for t in gz_shards:
        plain = t[: -len(".gz")]
        with open(plain, "wt") as fh:
            fh.write("legacy\tuncompressed\trow\n")
        legacy.append(plain)

    r2 = subprocess.run(base, capture_output=True, text=True, cwd=str(tmp_path))
    assert r2.returncode == 0, r2.stderr[-2000:]

    survivors = [q for q in legacy if os.path.exists(q)]
    assert not survivors, (
        "legacy uncompressed shard tracking survived the rewrite, so every later resume would "
        f"see an ambiguous shard and recompute it: {survivors}"
    )


def test_ambiguous_shard_with_both_encodings_is_not_reused(tmp_path):
    """Both suffixes present plus `.ok` -- reuse must be refused, not resolved by guesswork.

    The two-run test above covers the case where THIS run rewrites the shard, because the
    deletion happens at write time. It cannot cover the case that matters more: a shard the
    reuse check decides to SKIP never reaches the writer, so the deletion never runs. Both
    files can coexist from an interrupted flag flip, a killed run, or a version predating the
    deletion.

    The `.ok` checkpoint records no encoding, so with both files present nothing on disk says
    which run wrote which. The only safe answer is to recompute. Seeded directly rather than
    produced by two runs, because arriving at this state through the CLI requires the very
    deletion under test to be absent.
    """
    gtf, genome = _write_inputs(tmp_path)
    bam = _write_bam(tmp_path / "r.bam", gtf, genome)
    prefix = tmp_path / "out"

    # one real run to establish a genuine shard tree and its .ok checkpoints
    base = [sys.executable, LRAA, "--quant_only",
            "--bam", str(bam), "--gtf", str(gtf), "--genome", str(genome),
            "--output_prefix", str(prefix),
            "--no_cleanup", "--no_parallelize_contigs"]
    r1 = subprocess.run(base + [], capture_output=True, text=True,
                        cwd=str(tmp_path))
    assert r1.returncode == 0, r1.stderr[-2000:]

    gz_shards = [os.path.join(root, f)
                 for root, _dirs, files in os.walk(tmp_path)
                 for f in files
                 if f.endswith(".quant.tracking.gz") and "contigtmp" in root]
    assert gz_shards, "first run produced no gzipped shard tracking; fixture is inert"

    # forge the ambiguity: a plain sibling beside each .gz, with .ok left in place
    forged = []
    for gz in gz_shards:
        plain = gz[: -len(".gz")]
        with open(plain, "wt") as fh:
            fh.write("forged\tstale\trow\n")
        forged.append(plain)
        assert os.path.exists(os.path.join(os.path.dirname(gz),
                                           os.path.basename(gz).replace(
                                               ".quant.tracking.gz", ".ok"))), (
            "fixture expects an .ok checkpoint beside the shard"
        )

    # rerun WITH gzip: the current suffix matches a real file, and .ok is present, so the
    # only thing that can prevent reuse is noticing the plain sibling
    r2 = subprocess.run(base + [], capture_output=True, text=True,
                        cwd=str(tmp_path))
    assert r2.returncode == 0, r2.stderr[-2000:]

    out = r1.stdout + r1.stderr + r2.stdout + r2.stderr
    assert "invalidating its checkpoint" in out, (
        "ambiguous shard was reused silently; the merge could carry rows from either run"
    )
    for plain in forged:
        assert not os.path.exists(plain), (
            f"forged stale sibling survived the recompute: {plain}"
        )


def test_legacy_uncompressed_output_is_removed_from_the_prefix(tmp_path):
    """The deliverable must never exist under both names at once.

    Upgrading over a directory an older LRAA wrote leaves out.quant.tracking beside the
    out.quant.tracking.gz this version writes, one of them a run old, with nothing marking
    which is current. Anything globbing for either name -- a converter, a pipeline step, a WDL
    gather -- could read the stale one and never know.

    Seeded rather than produced, because no flag creates the uncompressed name any more.
    """
    gtf, genome = _write_inputs(tmp_path)
    bam = _write_bam(tmp_path / "r.bam", gtf, genome)
    prefix = tmp_path / "out"
    base = [sys.executable, LRAA, "--quant_only",
            "--bam", str(bam), "--gtf", str(gtf), "--genome", str(genome),
            "--output_prefix", str(prefix),
            "--no_parallelize_contigs"]

    plain = tmp_path / "out.LRAA.quant-only.quant.tracking"
    plain.write_text("legacy\tuncompressed\tfinal\trow\n")

    r = subprocess.run(base, capture_output=True, text=True, cwd=str(tmp_path))
    assert r.returncode == 0, r.stderr[-2000:]

    gz = tmp_path / "out.LRAA.quant-only.quant.tracking.gz"
    assert gz.exists(), "run produced no compressed tracking"
    assert not plain.exists(), (
        "the older version's uncompressed tracking survived beside the compressed one: a "
        "consumer globbing *.quant.tracking would read a stale file with no indication"
    )


def test_ambiguous_new_root_is_not_overridden_by_complete_legacy_root(tmp_path):
    """Ambiguity in either temp tree must veto reuse, not be outvoted by the other tree.

    Resume looks in a new temp root and falls back to a legacy one. Returning a single
    boolean per root let this happen: the new root is ambiguous so it declines, the legacy
    root is complete so it accepts, `complete_new or complete_old` is true, the shard is
    reused, the writer never runs, and the ambiguity survives to warn on every later run
    while resolving nothing. Completeness ORs across roots; ambiguity has to VETO across
    them, which a single boolean cannot express.

    Asserted on CONTENT, not on log text: a sentinel row is planted in the legacy shard, and
    the merged output must not contain it. If the shard were reused from the legacy tree the
    sentinel would be merged verbatim.
    """
    gtf, genome = _write_inputs(tmp_path)
    bam = _write_bam(tmp_path / "r.bam", gtf, genome)
    prefix = tmp_path / "out"
    base = [sys.executable, LRAA, "--quant_only",
            "--bam", str(bam), "--gtf", str(gtf), "--genome", str(genome),
            "--output_prefix", str(prefix),
            "--no_cleanup", "--no_parallelize_contigs"]

    r1 = subprocess.run(base, capture_output=True, text=True, cwd=str(tmp_path))
    assert r1.returncode == 0, r1.stderr[-2000:]

    # LRAA builds the new root as "__" + the output prefix it derived, which is an ABSOLUTE
    # path here, so the directory lands at <cwd>/__<abs prefix>.contigtmp. The legacy root is
    # the same name without the "__", i.e. a plain sibling of the outputs. Discover the new
    # one, and name the legacy one the way the resume check does -- deriving it by string
    # surgery on the discovered path produces a directory nothing ever looks at, and the test
    # then passes whatever the code does.
    roots = [os.path.join(r, d)
             for r, dirs, _f in os.walk(tmp_path)
             for d in dirs if d.endswith(".contigtmp")]
    assert roots, f"no .contigtmp temp root found under {tmp_path}"
    new_root = roots[0]
    old_root = str(tmp_path / "out.LRAA.quant-only.contigtmp")
    assert os.path.basename(new_root) == os.path.basename(old_root), (
        f"legacy root name does not match the discovered new root: {new_root} vs {old_root}"
    )

    # legacy tree: a complete copy of the finished shards
    shutil.copytree(new_root, old_root)

    # The sentinel must go where the MERGE reads from, or nothing can observe the
    # difference. The merge prefers the new root whenever it holds .ok/.quant.expr, so a
    # sentinel planted in the legacy tree is invisible whether the shard is reused or not.
    # Plant it in the NEW root and make the LEGACY root the ambiguous one: reuse then merges
    # the sentinel, while a veto recomputes the shard and overwrites it.
    SENTINEL = "SENTINEL_REUSED_ROW"
    new_tracks = [os.path.join(root, f)
                  for root, _dirs, files in os.walk(new_root)
                  for f in files if f.endswith(".quant.tracking.gz")]
    assert new_tracks, "new root holds no shard tracking; fixture is inert"
    for t in new_tracks:
        with gzip.open(t, "rt") as fh:
            lines = fh.readlines()
        with gzip.open(t, "wt") as fh:
            fh.writelines(lines[:1] + [f"{SENTINEL}\t0\t0\t0\t0\t0\t0\n"] + lines[1:])

    # legacy tree: ambiguous, via a plain sibling beside each .gz
    legacy_tracks = [os.path.join(root, f)
                     for root, _dirs, files in os.walk(old_root)
                     for f in files if f.endswith(".quant.tracking.gz")]
    assert legacy_tracks, "legacy copy holds no shard tracking; fixture is inert"
    for t in legacy_tracks:
        with open(t[: -len(".gz")], "wt") as fh:
            fh.write("forged\tstale\trow\n")

    r2 = subprocess.run(base, capture_output=True, text=True, cwd=str(tmp_path))
    assert r2.returncode == 0, r2.stderr[-2000:]

    merged = tmp_path / "out.LRAA.quant-only.quant.tracking.gz"
    assert merged.exists(), "second run produced no merged tracking"
    body = "".join(_content_rows(merged))
    assert SENTINEL not in body, (
        "the merged tracking carries the sentinel planted in the reused shard, so an "
        "ambiguous legacy root did not veto reuse of the new root"
    )


def test_shard_without_ok_is_refused_even_when_tracking_matches(tmp_path):
    """The state an interrupted ambiguity-veto leaves behind must be refused, not trusted.

    The ambiguity handler runs in the parent during the resume scan, while the recompute it
    justifies happens later in a worker. So there is a window: if the run dies in between, the
    shard is left with tracking under the current suffix and no sibling. Were the checkpoint
    still standing, the next run would read that as (complete, unambiguous) and silently merge
    a file whose provenance is exactly what the veto declared unknowable -- and the runs that
    get interrupted are the same ones that produce both-encoding shards.

    So the handler removes `.ok` as well. This test pins the INVARIANT that makes that safe --
    a shard whose outputs are all present but whose checkpoint is gone must be recomputed --
    and deliberately not the handler itself: it deletes `.ok` directly rather than reaching
    that state through an ambiguity veto.

    That limit is worth stating: this test passes with the handler's `.ok` removal reverted, and
    so does the ambiguity test, whose log line can be emitted without the unlink happening. The
    handler's own mutation is pinned instead by
    test_invalidating_an_ambiguous_shard_removes_the_checkpoint_first, which calls it directly
    and asserts the removal AND its order -- end-to-end tests cannot see it, because a worker
    re-stamps `.ok` before the run finishes.

    The sentinel is content -- reuse would merge it verbatim.
    """
    gtf, genome = _write_inputs(tmp_path)
    bam = _write_bam(tmp_path / "r.bam", gtf, genome)
    base = [sys.executable, LRAA, "--quant_only",
            "--bam", str(bam), "--gtf", str(gtf), "--genome", str(genome),
            "--output_prefix", str(tmp_path / "out"),
            "--no_cleanup", "--no_parallelize_contigs"]

    r1 = subprocess.run(base, capture_output=True, text=True, cwd=str(tmp_path))
    assert r1.returncode == 0, r1.stderr[-2000:]

    tracks = [os.path.join(root, f)
              for root, _dirs, files in os.walk(tmp_path)
              for f in files if f.endswith(".quant.tracking.gz") and "contigtmp" in root]
    assert tracks, "no shard tracking produced; fixture is inert"

    SENTINEL = "SENTINEL_POST_CRASH_ROW"
    removed_ok = []
    for t in tracks:
        with gzip.open(t, "rt") as fh:
            lines = fh.readlines()
        with gzip.open(t, "wt") as fh:
            fh.writelines(lines[:1] + [f"{SENTINEL}\t0\t0\t0\t0\t0\t0\n"] + lines[1:])
        ok = t.replace(".quant.tracking.gz", ".ok")
        if os.path.exists(ok):
            os.unlink(ok)
            removed_ok.append(ok)
    assert removed_ok, "no .ok checkpoints found to remove; fixture assumption is wrong"

    r2 = subprocess.run(base, capture_output=True, text=True, cwd=str(tmp_path))
    assert r2.returncode == 0, r2.stderr[-2000:]

    merged = tmp_path / "out.LRAA.quant-only.quant.tracking.gz"
    body = "".join(_content_rows(merged))
    assert SENTINEL not in body, (
        "a shard with no .ok checkpoint was reused: the state left by an interrupted "
        "ambiguity veto would be trusted instead of recomputed"
    )


def _load_lraa_module():
    """Import the LRAA driver script as a module so its helpers can be called directly.

    Everything else in this file drives LRAA as a subprocess, which is right for behaviour but
    cannot observe a mutation that a later step undoes. The `.ok` invalidation is exactly that:
    the parent removes it, then a worker re-stamps it, so end-to-end the file is present both
    before and after and no subprocess test can see the removal happen.
    """
    import importlib.util

    spec = importlib.util.spec_from_loader(
        "lraa_driver", importlib.machinery.SourceFileLoader("lraa_driver", LRAA)
    )
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def test_invalidating_an_ambiguous_shard_removes_the_checkpoint_first(tmp_path):
    """Pins the one mutation the subprocess tests cannot see, including its ORDER.

    Order is load-bearing: removing the sibling before `.ok` leaves, at the crash point, a
    checkpoint standing with its ambiguity evidence gone -- which is the defect this exists to
    fix. Recording the calls proves the sequence rather than just the end state.
    """
    import importlib.machinery  # noqa: F401  (used via _load_lraa_module)

    mod = _load_lraa_module()

    shard = tmp_path / "chr1.+"
    ok = tmp_path / "chr1.+.ok"
    current = tmp_path / "chr1.+.quant.tracking.gz"
    stale = tmp_path / "chr1.+.quant.tracking"
    for f in (ok, current, stale):
        f.write_text("x")

    calls = []
    real_unlink = os.unlink

    def recording_unlink(p):
        calls.append(str(p))
        return real_unlink(p)

    os.unlink = recording_unlink
    try:
        mod._invalidate_ambiguous_shard(str(shard), str(stale))
    finally:
        os.unlink = real_unlink

    assert not ok.exists(), "the checkpoint must be removed, or a crash leaves it trusted"
    assert not stale.exists(), "the stale sibling must be removed so the state converges"
    assert current.exists(), "the current-suffix tracking must be left alone"

    assert len(calls) == 2, f"expected two unlinks, got {calls}"
    assert calls[0].endswith(".ok"), (
        f"the checkpoint must be removed FIRST; order was {calls}"
    )
    assert calls[1] == str(stale), f"the sibling must be removed second; order was {calls}"


def test_invalidation_tolerates_a_missing_file_but_not_a_permission_error(tmp_path):
    """Idempotent for an already-absent file; loud for anything else.

    Swallowing OSError broadly would let an EPERM on `.ok` fall through to deleting the
    sibling, leaving the checkpoint standing with its evidence gone -- the original defect by
    another route.
    """
    import importlib.machinery  # noqa: F401

    mod = _load_lraa_module()

    shard = tmp_path / "chr2.-"
    # neither file exists: must not raise
    mod._invalidate_ambiguous_shard(str(shard), str(tmp_path / "chr2.-.quant.tracking"))

    ok = tmp_path / "chr3.-.ok"
    stale = tmp_path / "chr3.-.quant.tracking"
    ok.write_text("x")
    stale.write_text("x")

    real_unlink = os.unlink

    def failing_unlink(p):
        if str(p).endswith(".ok"):
            raise PermissionError(13, "denied")
        return real_unlink(p)

    os.unlink = failing_unlink
    try:
        with pytest.raises(PermissionError):
            mod._invalidate_ambiguous_shard(str(tmp_path / "chr3.-"), str(stale))
    finally:
        os.unlink = real_unlink

    assert stale.exists(), (
        "a failure removing the checkpoint must not go on to remove the evidence: that "
        "leaves the shard reusable with nothing recording that it is untrustworthy"
    )
