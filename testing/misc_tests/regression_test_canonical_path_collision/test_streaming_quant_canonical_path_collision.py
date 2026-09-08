"""Regression test for the streaming quant-only canonical-path collision.

HISTORY
-------
LRAA v0.31.0 (0.31.0-03b9c9a) aborted a whole 25-contig single-cell quant-only
run on this locus:

    StreamingQuant.py:149, AssignmentTable.build
    RuntimeError: two multipaths map to canonical path chr21+|TSS:...,POLYA:...:
                  cannot build one row set for it

Two distinct read-level multipaths canonicalized to one simple path. The guard
raises rather than merges, by design -- its own comment notes the correct merge
is undefined because the two multipaths carry different weights and different
mp ids -- so a single locus took down the entire run:
`1 of 5 chunks failed; refusing to merge a partial result`.

It reproduced on BOTH the ref-guided and the de novo arm, at 48 min and 23 min
of wall-clock respectively, and was only reachable through the *streaming*
quant-only path (`--stream_reads`) against an UNCOLLAPSED annotation -- here the
raw `init_gtf` straight out of initial discovery, which has never been through
merge_LRAA_GTFs.py / collapse_LRAA_GTF_by_splice_pattern.py. The cluster-guided
final quant over the same reads and the same chunk geometry never hit it,
because it quantifies the collapsed final GTF.

That last point is what makes this worth a permanent test: the uncollapsed
init GTF is a first-class, documented reuse path (`precomputed_init_gtf`), and
nothing else in the suite exercises streaming quant against one.

FIXTURE
-------
`fixture/` is chr21:6,280,000-6,320,000 from that run, translated to offset 0
(subtract 6,279,999) so the contig is 40 kb instead of 6.7 Mb. It carries the
whole read-sharing component that collided, comp-1119 (originally
chr21:6,286,342-6,313,221), with ~13 kb of margin on each side; only reads
falling ENTIRELY inside the window were kept, so the component's splice graph is
the one the failure was built from. 200 KB, and it runs in ~2 s.

Translation moves the literal path string in the error message, so this test
asserts on the failure MODE, never on those coordinates.

STATUS
------
    v0.31.0-03b9c9a   exit 1, RuntimeError    (verified 2026-09-07)
    v0.33.0-7e37488   exit 0, 33 quant rows   (verified 2026-09-07)
"""
import gzip
import os
import shutil
import subprocess
from pathlib import Path

import pytest

HERE = Path(__file__).parent
FIXTURE = HERE / "fixture"

# The collision is unreachable without BOTH of these: streaming, and an
# uncollapsed annotation. Dropping either is what made every other quant path
# in the pipeline pass while this one failed.
REQUIRED_FLAGS = ["--quant_only", "--stream_reads", "--stream_reads_rescue_unassigned"]


def _lraa_cmd():
    """LRAA from $LRAA_HOME, else the repo checkout this test sits in, else PATH."""
    if os.environ.get("LRAA_HOME"):
        cand = Path(os.environ["LRAA_HOME"]) / "LRAA"
        if cand.exists():
            return [str(cand)]
    for up in HERE.parents:
        if (up / "LRAA").is_file():
            return [str(up / "LRAA")]
    found = shutil.which("LRAA")
    if found:
        return [found]
    pytest.skip("LRAA executable not found; set LRAA_HOME")


@pytest.fixture(scope="module")
def quant_run(tmp_path_factory):
    work = tmp_path_factory.mktemp("canonical_path_collision")
    for f in FIXTURE.iterdir():
        shutil.copy(f, work / f.name)
    (work / "cfg.json").write_text('{"HiFi": true, "cpu_budget": 2}\n')

    cmd = _lraa_cmd() + [
        "--genome", "locus.fa",
        "--bam", "locus.strand.+.bam",
        "--no_chunk",
        "--gtf", "locus.gtf",
        "--bam_for_sg", "locus.plus.norm.bam",
        "--no_norm",
        "--num_total_reads", "81523164",
        "--cpu_budget", "1",
        "--output_prefix", "q",
        "--min_mapping_quality", "0",
        "--min_mapping_quality_for_final_quant", "0",
        "--HiFi",
        "--config_update", "cfg.json",
    ] + REQUIRED_FLAGS
    proc = subprocess.run(cmd, cwd=work, capture_output=True, text=True, timeout=900)
    return proc, work


def test_no_canonical_path_collision(quant_run):
    """The exact v0.31.0 failure: two multipaths reaching one canonical path."""
    proc, _ = quant_run
    blob = proc.stdout + proc.stderr
    assert "two multipaths map to canonical path" not in blob, (
        "v0.31.0 canonical-path collision has returned:\n"
        + "\n".join(l for l in blob.splitlines() if "canonical path" in l)
    )


def test_quant_only_run_succeeds(quant_run):
    """A collision anywhere in a chunk fails the whole run, so exit status is
    the other half of the assertion -- StreamingQuant could raise something new."""
    proc, _ = quant_run
    assert proc.returncode == 0, f"LRAA exited {proc.returncode}\n{proc.stderr[-3000:]}"


def test_component_that_collided_is_quantified(quant_run):
    """Not merely 'did not crash': the component whose multipaths collided has
    to come out the far side with rows. A fix that dropped or skipped the
    offending path would satisfy the two assertions above and still be wrong."""
    _, work = quant_run
    expr = work / "q.LRAA.quant-only.quant.expr"
    assert expr.exists(), "no quant.expr produced"
    rows = [l for l in expr.read_text().splitlines()
            if l and not l.startswith("#") and not l.startswith("gene_id")]
    assert len(rows) >= 20, f"expected the locus's transcripts, got {len(rows)} rows"
    collided = [r for r in rows if "comp-1119" in r]
    assert collided, "the component that collided (comp-1119) is absent from quant.expr"


def test_tracking_is_readable_and_nonempty(quant_run):
    """AssignmentTable builds the tracking rows, so tracking is the output the
    guard sits directly upstream of."""
    _, work = quant_run
    tracking = work / "q.LRAA.quant-only.quant.tracking.gz"
    assert tracking.exists(), "no quant.tracking.gz produced"
    with gzip.open(tracking, "rt") as fh:
        lines = [l for l in fh if l.strip()]
    assert len(lines) > 1, "tracking has no data rows"
