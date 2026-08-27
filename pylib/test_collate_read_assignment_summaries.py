#!/usr/bin/env python3

"""Tests for util/misc/collate_read_assignment_summaries.py.

N per-cluster read-assignment summaries into ONE table with three nested levels:
``chunk`` (a cluster's own worker rows, each naming the cut its reads came from),
``cluster`` (that file's TOTAL, re-emitted) and one ``all_clusters`` aggregate.

THE HEADLINE TEST IS THE AVERAGING ONE. Every fraction on the ``all_clusters`` row
must be summed numerator over summed denominator. Averaging the per-cluster
fractions instead produces a number that is wrong and entirely plausible, and no
count in the table contradicts it. The fixture therefore gives the two clusters
DELIBERATELY UNEQUAL sizes -- 141 reads against 16,167, the spread these runs
really show -- because with equal sizes the mean and the ratio of sums coincide
and an averaging bug would pass every assertion. Here they differ by 58x:
0.008646 against 0.5.

``test_an_averaged_all_clusters_row_fails_that_test`` is the negative control,
kept as a test rather than run once by hand: it forges the averaged file the bug
would have written and asserts the check above REJECTS it. Delete the
recomputation in ``_fill_fractions`` and the headline test fails; weaken the
headline test and the control fails.

The inputs are built by LRAA's own summary writer and then by
``ChunkedRun.merge_read_assignment_summaries``, never hand-typed, so the schema
under test is the one the pipeline actually emits -- including the ``chunk_id``
column the merger fills in from each unit's chunk.
"""

import csv
import importlib.util
import os
import sys
from importlib.machinery import SourceFileLoader
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "pylib"))

import ChunkedRun  # noqa: E402


def _load(name, relative):
    path = REPO / relative
    loader = SourceFileLoader(name, str(path))
    spec = importlib.util.spec_from_loader(loader.name, loader)
    module = importlib.util.module_from_spec(spec)
    loader.exec_module(module)
    return module


collate_mod = _load(
    "collate_read_assignment_summaries_under_test",
    "util/misc/collate_read_assignment_summaries.py",
)
lraa = _load("lraa_driver_collate_fixture", "LRAA")


# Two clusters, deliberately 115x apart in read count. See the module docstring:
# equal sizes would let an averaging bug pass.
SMALL_CLUSTER = "cluster_03"
BIG_CLUSTER = "cluster_17"

CLUSTERS = [
    (
        SMALL_CLUSTER,
        [
            ("chrT_00", {"reads_total": 100, "reads_kept_genome": 100}),
            ("chrT_01", {"reads_total": 41, "reads_kept_genome": 41}),
        ],
    ),
    (
        BIG_CLUSTER,
        [
            (
                "chrT_00",
                {
                    "reads_total": 16000,
                    "reads_kept_genome": 0,
                    "reads_selected_tx_total": 8,
                    "reads_rescue_requested": 16000,
                    "reads_rescue_rescued": 8,
                    "reads_rescue_unrescued": 15992,
                    "rescue_alignment_rejections": "low_identity=3,short_overlap=1",
                },
            ),
            (
                "chrT_01",
                {
                    "reads_total": 167,
                    "reads_kept_genome": 0,
                    "reads_rescue_requested": 167,
                    "reads_rescue_unrescued": 167,
                    "rescue_alignment_rejections": "low_identity=2",
                },
            ),
        ],
    ),
]

SMALL_TOTAL = 141
BIG_TOTAL = 16167
ALL_TOTAL = SMALL_TOTAL + BIG_TOTAL


def _cluster_summary(root, cluster_id, chunks):
    """One cluster's merged summary, by the route a real run takes to it.

    LRAA's per-unit writer, LRAA's own single-unit merge (so each unit file holds a
    worker row AND its own TOTAL, which is what stage 6 has to skip), then
    ``ChunkedRun.merge_read_assignment_summaries`` over units carrying their
    ``chunk_id`` -- which is the only thing that puts a chunk id in the table.
    """

    cdir = root / cluster_id
    cdir.mkdir()
    units = []
    for chunk_id, stats in chunks:
        prefix = str(cdir / chunk_id)
        worker = prefix + ".read_assignment.summary.worker.tsv"
        lraa._write_read_assignment_summary(worker, "chrT", "+", stats)
        lraa._merge_read_assignment_summary_files(
            [worker], prefix + ".read_assignment.summary.tsv"
        )
        os.remove(worker)
        units.append(
            {
                # unit_id encodes the ORIENTATION as well as the cut, which is why
                # it is not the chunk label.
                "unit_id": "{}_plus".format(chunk_id),
                "chunk_id": chunk_id,
                "offset": 0,
                "quant_prefix": prefix,
            }
        )

    merged_dir = cdir / "merged"
    merged_dir.mkdir()
    merged = ChunkedRun.merge_read_assignment_summaries(str(merged_dir), units)
    assert merged is not None
    return merged["path"]


def _rows(path):
    with open(path, "rt", newline="") as fh:
        lines = [line for line in fh if not line.startswith("#")]
    return list(csv.DictReader(lines, delimiter="\t"))


def _header(path):
    with open(path, "rt", newline="") as fh:
        for line in fh:
            if not line.startswith("#"):
                # rstrip of BOTH, because the per-invocation summaries LRAA writes
                # are CRLF (csv.DictWriter's default) while the collated file is LF.
                return line.rstrip("\r\n").split("\t")
    raise AssertionError("no header in {}".format(path))


def _one(rows, level):
    picked = [r for r in rows if r["level"] == level]
    assert len(picked) == 1, picked
    return picked[0]


@pytest.fixture(scope="module")
def collated(tmp_path_factory):
    root = tmp_path_factory.mktemp("collate")
    pairs = [
        (cluster_id, _cluster_summary(root, cluster_id, chunks))
        for cluster_id, chunks in CLUSTERS
    ]
    out = root / "collated.read_assignment.summary.tsv"
    result = collate_mod.collate(pairs, str(out), population="cluster_thinned")
    assert result["clusters"] == 2
    return {"path": str(out), "pairs": pairs, "result": result, "root": root}


# -- the layout the contract names


def test_the_leading_columns_are_level_cluster_chunk(collated):
    """``level``, ``cluster_id``, ``chunk_id``, then the inherited schema.

    ``row_type`` is gone rather than kept beside ``level``: two columns saying what
    a row is, in two vocabularies, is one column too many. ``population`` is
    appended LAST so the inherited contig_acc-onward columns stay contiguous and in
    their original order.
    """

    header = _header(collated["path"])
    assert header[:3] == ["level", "cluster_id", "chunk_id"]
    assert header[-1] == "population"
    assert "row_type" not in header

    source = _header(collated["pairs"][0][1])
    tail = [f for f in source if f not in ("row_type", "chunk_id")]
    assert tail[0] == "contig_acc"
    assert header[3:-1] == tail


def test_chunk_rows_name_their_cluster_and_their_cut(collated):
    rows = _rows(collated["path"])
    chunk_rows = [r for r in rows if r["level"] == "chunk"]

    # Two chunks per cluster, each contributing one worker row here.
    assert len(chunk_rows) == 4
    assert {(r["cluster_id"], r["chunk_id"]) for r in chunk_rows} == {
        (SMALL_CLUSTER, "chrT_00"),
        (SMALL_CLUSTER, "chrT_01"),
        (BIG_CLUSTER, "chrT_00"),
        (BIG_CLUSTER, "chrT_01"),
    }
    # The per-chunk detail this table exists for: without chunk_id every one of
    # these rows reads "chrT +" and the outlier chunk is invisible.
    small = {
        r["chunk_id"]: int(r["reads_total"])
        for r in chunk_rows
        if r["cluster_id"] == SMALL_CLUSTER
    }
    assert small == {"chrT_00": 100, "chrT_01": 41}


def test_cluster_rows_carry_no_chunk_and_all_clusters_carries_no_cluster(collated):
    rows = _rows(collated["path"])
    cluster_rows = [r for r in rows if r["level"] == "cluster"]
    assert len(cluster_rows) == 2
    assert {r["cluster_id"] for r in cluster_rows} == {SMALL_CLUSTER, BIG_CLUSTER}
    assert [r["chunk_id"] for r in cluster_rows] == ["", ""]

    all_row = _one(rows, "all_clusters")
    assert all_row["cluster_id"] == ""
    assert all_row["chunk_id"] == ""

    # Last row, so the aggregate is read after what it aggregates.
    assert rows[-1]["level"] == "all_clusters"


def test_the_comment_lines_forbid_summing_and_forbid_the_pseudobulk_diff(collated):
    """The two traps this table can lead a reader into, both stated in the file.

    Prose rather than a schema, because a reader who has the file and not this
    repo has nowhere else to learn either one.
    """

    with open(collated["path"], "rt") as fh:
        comments = [line for line in fh if line.startswith("#")]
    assert comments, "the collated file must say what its levels mean"
    text = " ".join(comments)

    assert "MUST NOT BE SUMMED" in text
    assert "chunk | cluster | all_clusters" in text

    # The all_clusters row is a sum of independently-thinned per-cluster
    # populations, so it is not a library size and not the pseudobulk figure.
    assert "THINNED INDEPENDENTLY" in text
    assert "MUST NOT be compared to a pooled-normalized or pseudobulk" in text
    assert "population column" in text


def test_the_population_label_is_on_every_row_and_derived_on_the_aggregate(collated):
    """Machine-visible non-comparability, and an aggregate that cannot lie about it."""

    rows = _rows(collated["path"])
    for row in rows:
        if row["level"] == "all_clusters":
            # DERIVED from the row label, not taken separately.
            assert row["population"] == "sum_of_cluster_thinned"
        else:
            assert row["population"] == "cluster_thinned"


def test_the_population_label_defaults_to_unspecified_rather_than_guessing(
    collated, tmp_path
):
    """A wrong label is worse than an absent one.

    Defaulting to ``cluster_thinned`` would mislabel exactly the file the column
    exists to distinguish -- the pooled-normalized init/pseudobulk one.
    """

    out = tmp_path / "unlabelled.tsv"
    collate_mod.collate(collated["pairs"], str(out))
    rows = _rows(str(out))
    assert _one(rows, "all_clusters")["population"] == "sum_of_unspecified"
    assert {r["population"] for r in rows if r["level"] == "cluster"} == {
        "unspecified"
    }


# -- the arithmetic


def test_all_clusters_counts_equal_the_sum_of_the_cluster_rows(collated):
    """Every summable counter, not just reads_total."""

    rows = _rows(collated["path"])
    cluster_rows = [r for r in rows if r["level"] == "cluster"]
    all_row = _one(rows, "all_clusters")

    for key in collate_mod.SUMMABLE_FIELDS:
        expected = sum(int(r[key] or 0) for r in cluster_rows)
        assert int(all_row[key]) == expected, key

    assert int(all_row["reads_total"]) == ALL_TOTAL
    assert int(all_row["reads_kept_genome"]) == SMALL_TOTAL
    assert int(all_row["reads_rescue_requested"]) == BIG_TOTAL


def test_cluster_rows_are_the_input_totals_re_emitted(collated):
    """A cluster row is that file's own TOTAL, not a re-derivation of it."""

    rows = {
        r["cluster_id"]: r for r in _rows(collated["path"]) if r["level"] == "cluster"
    }
    for cluster_id, path in collated["pairs"]:
        source = [r for r in _rows(path) if r["row_type"] == "TOTAL"]
        assert len(source) == 1
        for field, value in source[0].items():
            if field in ("row_type", "chunk_id"):
                continue
            assert rows[cluster_id][field] == value, (cluster_id, field)


def _assert_fractions_are_recomputed(path):
    """Every ``all_clusters`` fraction is summed-over-summed.

    Shared with the negative control below, which forges the averaged row this
    rejects. Kept as a helper for exactly that reason: the control has to exercise
    the same assertions, or it proves nothing about them.
    """

    rows = _rows(path)
    cluster_rows = [r for r in rows if r["level"] == "cluster"]
    all_row = _one(rows, "all_clusters")

    totals = {
        key: sum(int(r[key] or 0) for r in cluster_rows)
        for key in collate_mod.SUMMABLE_FIELDS
    }
    denominator = totals["reads_total"]
    requested = totals["reads_rescue_requested"]

    for field in all_row:
        if not field.startswith("frac_"):
            continue
        if field == "frac_rescue_rescued_of_requested":
            expected = float(totals["reads_rescue_rescued"]) / float(requested)
        elif field == "frac_rescue_unrescued_of_requested":
            expected = float(totals["reads_rescue_unrescued"]) / float(requested)
        else:
            counter = "reads_" + field[len("frac_") :]
            assert counter in totals, field
            expected = float(totals[counter]) / float(denominator)

        got = float(all_row[field])
        assert got == pytest.approx(expected, abs=5e-7), field

        # and is NOT the mean of the per-cluster fractions. Checked per field
        # rather than once, because the two agree for any counter that happens to
        # be zero in every cluster and the point is the ones that are not.
        mean = sum(float(r[field]) for r in cluster_rows) / float(len(cluster_rows))
        if abs(mean - expected) > 1e-6:
            assert got != pytest.approx(mean, abs=1e-6), field


def test_all_clusters_fractions_are_recomputed_not_averaged(collated):
    _assert_fractions_are_recomputed(collated["path"])

    # The two numbers, named, so a future reader sees the gap the unequal cluster
    # sizes create: 141 kept of 16,308 against the mean of 1.0 and 0.0.
    all_row = _one(_rows(collated["path"]), "all_clusters")
    assert float(all_row["frac_kept_genome"]) == pytest.approx(
        SMALL_TOTAL / ALL_TOTAL, abs=5e-7
    )
    assert float(all_row["frac_kept_genome"]) != pytest.approx(0.5, abs=1e-3)
    # 8 rescued of 16,167 requested, against a mean of that fraction and a cluster
    # that requested nothing.
    assert float(all_row["frac_rescue_rescued_of_requested"]) == pytest.approx(
        8.0 / BIG_TOTAL, abs=5e-7
    )


def test_an_averaged_all_clusters_row_fails_that_test(collated, tmp_path):
    """THE NEGATIVE CONTROL. Swap recomputation for a mean and the check above bites.

    The averaged file is forged from the real one rather than produced by patching
    the util, so this pins the ASSERTIONS' sensitivity and does not depend on the
    shape of the code that would carry the bug.
    """

    rows = _rows(collated["path"])
    header = _header(collated["path"])
    cluster_rows = [r for r in rows if r["level"] == "cluster"]

    forged = tmp_path / "averaged.tsv"
    with open(forged, "wt", newline="") as ofh:
        writer = csv.DictWriter(
            ofh, fieldnames=header, delimiter="\t", lineterminator="\n"
        )
        writer.writeheader()
        for row in rows:
            out = dict(row)
            if row["level"] == "all_clusters":
                for field in header:
                    if field.startswith("frac_"):
                        out[field] = "{:.6f}".format(
                            sum(float(r[field]) for r in cluster_rows)
                            / float(len(cluster_rows))
                        )
            writer.writerow(out)

    with pytest.raises(AssertionError):
        _assert_fractions_are_recomputed(str(forged))


def test_rescue_alignment_rejections_sum_by_key(collated):
    """A ``key=count`` list, summed by key rather than concatenated."""

    all_row = _one(_rows(collated["path"]), "all_clusters")
    assert all_row["rescue_alignment_rejections"] == "low_identity=5,short_overlap=1"


# -- refusals


def test_a_missing_input_is_refused_naming_it(collated, tmp_path, capsys):
    pairs = list(collated["pairs"]) + [("cluster_99", str(tmp_path / "gone.tsv"))]
    out = tmp_path / "never_written.tsv"

    code = collate_mod.main(
        [arg for cid, path in pairs for arg in ("--cluster_id", cid, "--summary", path)]
        + ["--output", str(out)]
    )

    assert code == 2
    err = capsys.readouterr().err
    assert "cluster_99" in err
    assert "gone.tsv" in err
    # A partial file would be indistinguishable from a successful task.
    assert not out.exists()


def test_unpaired_arguments_are_refused(collated, tmp_path, capsys):
    out = tmp_path / "unpaired.tsv"
    code = collate_mod.main(
        [
            "--summary",
            collated["pairs"][0][1],
            "--summary",
            collated["pairs"][1][1],
            "--cluster_id",
            SMALL_CLUSTER,
            "--output",
            str(out),
        ]
    )
    assert code == 2
    assert "paired by position" in capsys.readouterr().err
    assert not out.exists()


def test_a_repeated_cluster_id_is_refused(collated, tmp_path):
    path = collated["pairs"][0][1]
    with pytest.raises(collate_mod.CollationError, match=SMALL_CLUSTER):
        collate_mod.collate(
            [(SMALL_CLUSTER, path), (SMALL_CLUSTER, path)],
            str(tmp_path / "dup.tsv"),
        )


def test_an_input_without_a_total_row_is_refused(collated, tmp_path):
    """Its TOTAL *is* the cluster row; deriving a replacement would hide the damage."""

    truncated = tmp_path / "no_total.tsv"
    rows = _rows(collated["pairs"][0][1])
    header = _header(collated["pairs"][0][1])
    with open(truncated, "wt", newline="") as ofh:
        writer = csv.DictWriter(
            ofh, fieldnames=header, delimiter="\t", lineterminator="\n"
        )
        writer.writeheader()
        for row in rows:
            if row["row_type"] != "TOTAL":
                writer.writerow(row)

    with pytest.raises(collate_mod.CollationError, match="TOTAL"):
        collate_mod.collate(
            [(SMALL_CLUSTER, str(truncated))], str(tmp_path / "out.tsv")
        )


def test_inputs_disagreeing_on_columns_are_refused(collated, tmp_path):
    narrowed = tmp_path / "narrow.tsv"
    rows = _rows(collated["pairs"][0][1])
    header = [f for f in _header(collated["pairs"][0][1]) if f != "frac_kept_genome"]
    with open(narrowed, "wt", newline="") as ofh:
        writer = csv.DictWriter(
            ofh, fieldnames=header, delimiter="\t", lineterminator="\n"
        )
        writer.writeheader()
        for row in rows:
            writer.writerow({f: row[f] for f in header})

    with pytest.raises(collate_mod.CollationError, match="cannot share a table"):
        collate_mod.collate(
            [
                (SMALL_CLUSTER, collated["pairs"][0][1]),
                (BIG_CLUSTER, str(narrowed)),
            ],
            str(tmp_path / "out.tsv"),
        )


# -- the pseudobulk shape


def test_no_all_clusters_omits_the_aggregate(collated, tmp_path):
    """The init/pseudobulk phase, reported on its own.

    One input, so an ``all_clusters`` row would only restate the single cluster
    row -- and the whole point of the separate file is that a pooled-normalized
    total is NOT commensurable with the per-cluster ones, which is also why this
    call labels its population differently.
    """

    out = tmp_path / "pseudobulk.tsv"
    code = collate_mod.main(
        [
            "--summary",
            collated["pairs"][1][1],
            "--cluster_id",
            "pseudobulk",
            "--output",
            str(out),
            "--no_all_clusters",
            "--population",
            "pooled_thinned",
        ]
    )
    assert code == 0

    rows = _rows(str(out))
    assert {r["level"] for r in rows} == {"chunk", "cluster"}
    assert _one(rows, "cluster")["cluster_id"] == "pseudobulk"
    assert int(_one(rows, "cluster")["reads_total"]) == BIG_TOTAL
    assert {r["population"] for r in rows} == {"pooled_thinned"}
