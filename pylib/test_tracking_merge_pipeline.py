#!/usr/bin/env python3

"""The per-shard tracking merge in the WDL equals the python merge it replaced.

mergeQuantResults (LRAA.wdl) carries it. (LRAA_merge_trackings in
LRAA-cell_cluster_guided.wdl carried the same text until the per-cluster final quant was
replaced by one pooled quant, which leaves nothing to merge across clusters.) These tests
lift it out of the WDL and run it against a copy of the former python loop, on the shard shapes that break a naive pipeline: an empty first shard,
a plain-text shard, one input, more shards than cores, comment-only shards, and a schema
mismatch, which must still fail.
"""

import gzip
import shutil
import subprocess
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parents[1]
BLOCKS = {
    "mergeQuantResults": (REPO / "WDL" / "LRAA.wdl", "quantTrackingFiles", "~{outputFilePrefix}.quant.tracking.gz"),
}
HEADER = "gene_id\ttranscript_id\tread_name\tfrac\n"
MERGE_HEADER = "# LRAA version vX\n"

needs_tools = pytest.mark.skipif(shutil.which("pigz") is None, reason="pigz not installed")


def _block(task):
    wdl, files, out = BLOCKS[task]
    src = wdl.read_text()
    start = src.index("task {} {{".format(task))
    end = src.index(">>>", start)
    body = src[src.index("    set -o pipefail\n    mkdir parts", start) : end]
    return body, files, out


def _run(task, tmp_path, shards, cpu=4):
    body, files, out = _block(task)
    lst = tmp_path / "inputs.list"
    lst.write_text("".join(str(p) + "\n" for p in shards))
    (tmp_path / "merge_header.txt").write_text(MERGE_HEADER)
    script = (
        body.replace("~{write_lines(" + files + ")}", str(lst))
        .replace("~{merge_cpu}", str(cpu))
        .replace(out, "merged.tracking.gz")
    )
    assert "~{" not in script
    (tmp_path / "merge.sh").write_text("set -e\n" + script)
    return subprocess.run(["bash", "merge.sh"], cwd=tmp_path, capture_output=True, text=True)


def _old_python_merge(shards):
    """The loop the WDLs held through v0.47.0."""

    out = [MERGE_HEADER]
    wrote_header = False
    expected_header = None
    for tracking_file in shards:
        openf = gzip.open if str(tracking_file).split(".")[-1] == "gz" else open
        with openf(tracking_file, "rt") as fh:
            header = None
            for line in fh:
                if line.startswith("#"):
                    continue
                if header is None:
                    header = line
                    if expected_header is None:
                        expected_header = header
                    elif header != expected_header:
                        raise RuntimeError("schemas differ")
                    if not wrote_header:
                        out.append(header)
                        wrote_header = True
                    continue
                out.append(line)
    return "".join(out)


def _shard(path, text):
    if str(path).endswith(".gz"):
        with gzip.open(path, "wt") as fh:
            fh.write(text)
    else:
        path.write_text(text)
    return path


def _rows(tag, n):
    return "".join("g{0}\tt{0}\tread{1}\t0.5\n".format(tag, i) for i in range(n))


CASES = {
    "plain_case": lambda d: [
        _shard(d / "a.gz", "# c\n" + HEADER + _rows("a", 5)),
        _shard(d / "b.gz", "# c\n" + HEADER + _rows("b", 5)),
    ],
    "empty_first_shard": lambda d: [
        _shard(d / "a.gz", ""),
        _shard(d / "b.gz", HEADER + _rows("b", 3)),
        _shard(d / "c.gz", HEADER + _rows("c", 3)),
    ],
    "comment_only_first_shard": lambda d: [
        _shard(d / "a.gz", "# only a comment\n"),
        _shard(d / "b.gz", "# c\n" + HEADER + _rows("b", 3)),
    ],
    "header_only_shard": lambda d: [
        _shard(d / "a.gz", HEADER + _rows("a", 3)),
        _shard(d / "b.gz", HEADER),
        _shard(d / "c.gz", HEADER + _rows("c", 3)),
    ],
    "plain_text_shard": lambda d: [
        _shard(d / "a.gz", HEADER + _rows("a", 3)),
        _shard(d / "b.tracking", "# c\n" + HEADER + _rows("b", 3)),
    ],
    "one_input": lambda d: [_shard(d / "a.gz", "# c\n" + HEADER + _rows("a", 3))],
    "all_empty": lambda d: [_shard(d / "a.gz", ""), _shard(d / "b.gz", "# c\n")],
    "comment_mid_shard": lambda d: [
        _shard(d / "a.gz", HEADER + _rows("a", 2) + "# mid\n" + _rows("x", 2)),
    ],
    "more_shards_than_cores": lambda d: [
        _shard(d / "s{}.gz".format(i), "# c\n" + HEADER + _rows(i, 40)) for i in range(25)
    ],
}


@needs_tools
@pytest.mark.parametrize("task", sorted(BLOCKS))
@pytest.mark.parametrize("case", sorted(CASES))
def test_same_bytes_as_the_python_merge(task, case, tmp_path):
    shards = CASES[case](tmp_path)
    result = _run(task, tmp_path, shards)
    assert result.returncode == 0, result.stderr
    with gzip.open(tmp_path / "merged.tracking.gz", "rt") as fh:
        assert fh.read() == _old_python_merge(shards)


@needs_tools
@pytest.mark.parametrize("task", sorted(BLOCKS))
def test_a_schema_mismatch_still_fails(task, tmp_path):
    shards = [
        _shard(tmp_path / "a.gz", HEADER + _rows("a", 2)),
        _shard(tmp_path / "b.gz", "other\tcolumns\n" + _rows("b", 2)),
    ]
    result = _run(task, tmp_path, shards)
    assert result.returncode != 0
    assert "different schemas" in result.stderr


@needs_tools
@pytest.mark.parametrize("task", sorted(BLOCKS))
def test_a_corrupt_shard_fails(task, tmp_path):
    good = _shard(tmp_path / "a.gz", HEADER + _rows("a", 2000))
    data = good.read_bytes()
    bad = tmp_path / "b.gz"
    bad.write_bytes(data[: len(data) // 2])
    result = _run(task, tmp_path, [good, bad])
    assert result.returncode != 0

