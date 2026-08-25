#!/usr/bin/env python3

"""Tests for util/misc/merge_chunk_outputs.py, the standalone stage-6 merger.

There is exactly ONE merge implementation: ``ChunkedRun.merge_and_translate``.
This script exists only to give it a file-based input contract (a JSON
manifest) and a CLI, so a workflow engine can invoke stage 6 as its own task.
``pylib/test_chunked_tracking_streaming_merge.py`` and
``pylib/test_chunked_merge_whole_contig_frame.py`` already prove the merge
itself is correct; these tests prove the WRAPPER does not diverge from it and
that the manifest/grouping/combine machinery this script adds is sound.

Four things are proven:

1. The CLI round-trips a manifest to EXACTLY what a direct
   ``ChunkedRun.merge_and_translate`` call produces, in both discovery modes.
2. A header mismatch still refuses through the CLI, and refuses before any
   output file exists (only the empty ``merged/`` directory is created).
3. Manifest and group validation name the offending unit/field/group rather
   than surfacing a bare ``KeyError`` deep inside the merge.
4. The two GROUPED-MERGE identities the module docstring claims:
   quant.tracking is byte-identical to a global merge once gathered with
   ``heapq.merge`` (not by concatenation), and quant.expr is byte-identical
   to a global merge once recombined with ``--combine-expr`` (not by
   concatenation either -- its TPM is scope-relative).
"""

import gzip
import heapq
import importlib.util
import json
import os
import subprocess
import sys
from importlib.machinery import SourceFileLoader
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "pylib"))

import ChunkedRun  # noqa: E402
import Util_funcs  # noqa: E402


def _load(name, relative):
    path = REPO / relative
    loader = SourceFileLoader(name, str(path))
    spec = importlib.util.spec_from_loader(loader.name, loader)
    module = importlib.util.module_from_spec(spec)
    loader.exec_module(module)
    return module


mco = _load("merge_chunk_outputs_under_test", "util/misc/merge_chunk_outputs.py")

# LRAA's own read-assignment summary writers. Stage 6 now REQUIRES a summary from
# every unit (``ChunkedRun.merge_read_assignment_summaries``), because every LRAA
# path writes one -- including the ``run_quant_only`` early returns that quantify
# nothing and report the reads they saw with zero assigned. Loaded rather than
# reimplemented: a hand-copied 28-column schema would be free to drift from the
# one LRAA writes, which is the disagreement the merge's field-name check exists
# to catch.
lraa = _load("lraa_driver_summary_fixture", "LRAA")


def _write_read_assignment_summary(prefix, reads_total):
    """The per-unit summary a real quant unit leaves beside its quant.expr.

    Written through LRAA's worker writer and then its own merger, because that is
    how a unit's file comes to hold a ``worker`` row AND its own ``TOTAL``: stage 6
    has to skip the latter or it counts every unit twice.
    """

    worker = prefix + ".read_assignment.summary.worker.tsv"
    lraa._write_read_assignment_summary(
        worker,
        "chrT",
        "+",
        {"reads_total": reads_total, "reads_kept_genome": reads_total},
    )
    lraa._merge_read_assignment_summary_files(
        [worker], prefix + ".read_assignment.summary.tsv"
    )
    os.remove(worker)


EXPR_HEADER = [
    "gene_id",
    "transcript_id",
    "uniq_reads",
    "all_reads",
    "isoform_fraction",
    "unique_gene_read_fraction",
    "TPM",
    "exons",
    "introns",
    "splice_hash_code",
    "RPM_total_reads",
]
TRACKING_HEADER = [
    "gene_id",
    "transcript_id",
    "transcript_splice_hash_code",
    "num_exons",
    "mp_id",
    "read_name",
    "frac_assigned",
    "read_weight",
]

# Two models per unit, distinct all_reads so the TPM rebase has something to
# actually rebase, and distinct enough between units that a merge that failed
# to isolate scopes (or a combine step reusing the wrong total) would show up
# as a row difference rather than being masked by symmetry.
MODELS = [
    ("g1", "t1", "chrT:(+)[[100, 200], [300, 400]]", "chrT:(+)[[201, 299]]"),
    ("g2", "t2", "chrT:(+)[[500, 600]]", ""),
]


def _write_unit(root, unit_id, offset, all_reads_by_model, read_tokens):
    """One quant unit's chunk-local ``quant.expr``, gzipped ``quant.tracking``,
    ``.gtf`` (``merge_discovery_gtf`` needs it) and ``read_assignment.summary.tsv``
    (stage 6 requires one from every unit).

    ``all_reads_by_model`` varies the counts per unit, deliberately, so the
    unrebased TPM differs from the correctly-rebased one -- a bug that used
    the wrong scope's total would produce a value that merely looks plausible
    rather than one that is obviously wrong.

    ``read_tokens`` are this unit's ``read_name`` values, chosen by the caller
    rather than derived from ``unit_id``: the tracking merge sorts on
    ``(read_name, transcript_id, gene_id)``, so if a read name embedded its
    own unit id, every unit's rows would already be contig-major and
    concatenating grouped output would accidentally match a global sort --
    masking exactly the identity the grouped-merge tests exist to check.
    """

    prefix = str(root / unit_id)
    reads = [all_reads_by_model[i] for i in range(len(MODELS))]
    total = float(sum(reads)) or 1.0

    with open(prefix + ".quant.expr", "wt") as fh:
        print("\t".join(EXPR_HEADER), file=fh)
        for (gene, tx, exons, introns), ar in zip(MODELS, reads):
            hash_code = Util_funcs.get_hash_code(introns) if introns else tx
            # This unit's own chunk-local TPM -- deliberately NOT the final
            # value, so a test that forgot to rebase would fail loudly.
            local_tpm = ar / total * 1e6
            print(
                "\t".join(
                    [
                        gene, tx, str(ar), "{:.1f}".format(ar),
                        "0.500", "1.000", "{:.3f}".format(local_tpm),
                        exons, introns, hash_code, "50.000",
                    ]
                ),
                file=fh,
            )

    with gzip.open(prefix + ".quant.tracking.gz", "wt") as fh:
        print("\t".join(TRACKING_HEADER), file=fh)
        for i, read_name in enumerate(read_tokens):
            gene, tx, _, introns = MODELS[i % len(MODELS)]
            hash_code = Util_funcs.get_hash_code(introns) if introns else tx
            print(
                "\t".join(
                    [
                        gene, tx, hash_code, "2", "mp{}".format(i),
                        read_name, "1.000", "1.000",
                    ]
                ),
                file=fh,
            )

    with open(prefix + ".gtf", "wt") as fh:
        print("# LRAA version test", file=fh)
        for gene, tx, exons, _ in MODELS:
            attrs = 'gene_id "{}"; transcript_id "{}";'.format(gene, tx)
            print(
                "\t".join(("chrT", "LRAA", "transcript", "1", "10", ".", "+", ".", attrs)),
                file=fh,
            )
            print(
                "\t".join(("chrT", "LRAA", "exon", "1", "10", ".", "+", ".", attrs)),
                file=fh,
            )

    # One read per tracking row, which is what a unit's summary reports and what
    # stage 6 now requires of every unit.
    _write_read_assignment_summary(prefix, len(read_tokens))

    return {"unit_id": unit_id, "offset": offset, "quant_prefix": prefix}


@pytest.fixture
def two_group_units(tmp_path):
    """Four units, two contigs x two orientations, nonzero offsets.

    Nonzero offsets exercise the coordinate-translation path (not the
    whole-contig skip); the contig split gives ``write_manifest``'s default
    grouping something real to key on; and the read-name tokens interleave
    chrA/chrB in the global sort order (m01 A, m02 A, m03 B, m04 B, m05 A,
    ...) so a global merge's row order and a concat-of-groups' row order
    provably differ -- the fixture would otherwise accidentally pass the
    grouped-identity tests for the wrong reason.
    """

    root = tmp_path / "chunks"
    root.mkdir()
    units = [
        _write_unit(root, "chrA_00_plus", 1000, [3, 7], ["m02", "m06", "m10"]),
        _write_unit(root, "chrA_00_minus", 1000, [5, 1], ["m01", "m05"]),
        _write_unit(root, "chrB_00_plus", 2000, [11, 2], ["m04", "m08"]),
        _write_unit(root, "chrB_00_minus", 2000, [4, 4], ["m03", "m07", "m09"]),
    ]
    return units


def _read_text(path):
    opener = gzip.open if path.endswith(".gz") else open
    with opener(path, "rt") as fh:
        return fh.read()


class TestManifestValidation:
    def test_missing_units_key_names_the_file(self, tmp_path):
        p = tmp_path / "m.json"
        p.write_text(json.dumps({"not_units": []}))
        with pytest.raises(ValueError, match=r"'units' array"):
            mco.load_manifest(str(p))

    def test_empty_units_is_rejected(self, tmp_path):
        p = tmp_path / "m.json"
        p.write_text(json.dumps({"units": []}))
        with pytest.raises(ValueError, match="no units"):
            mco.load_manifest(str(p))

    def test_missing_field_names_the_unit_and_field(self, tmp_path):
        p = tmp_path / "m.json"
        p.write_text(
            json.dumps({"units": [{"unit_id": "chrA_plus", "offset": 0}]})
        )
        with pytest.raises(ValueError, match=r"chrA_plus.*quant_prefix"):
            mco.load_manifest(str(p))

    def test_non_integer_offset_is_rejected(self, tmp_path):
        p = tmp_path / "m.json"
        p.write_text(
            json.dumps(
                {
                    "units": [
                        {
                            "unit_id": "chrA_plus",
                            "quant_prefix": "/x",
                            "offset": "0",
                        }
                    ]
                }
            )
        )
        with pytest.raises(ValueError, match="non-integer offset"):
            mco.load_manifest(str(p))

    def test_order_is_preserved(self, tmp_path):
        p = tmp_path / "m.json"
        ids = ["z_plus", "a_plus", "m_plus"]
        p.write_text(
            json.dumps(
                {
                    "units": [
                        {"unit_id": u, "quant_prefix": "/x/" + u, "offset": 0}
                        for u in ids
                    ]
                }
            )
        )
        units = mco.load_manifest(str(p))
        assert [u["unit_id"] for u in units] == ids


class TestGroupSelection:
    def test_unknown_group_lists_the_groups_present(self):
        units = [
            {"unit_id": "a", "group": "chrA"},
            {"unit_id": "b", "group": "chrB"},
        ]
        with pytest.raises(ValueError, match=r"chrA.*chrB|chrB.*chrA"):
            mco.select_group(units, "chrZ")

    def test_matching_group_keeps_order(self):
        units = [
            {"unit_id": "a", "group": "chrA"},
            {"unit_id": "b", "group": "chrB"},
            {"unit_id": "c", "group": "chrA"},
        ]
        kept = mco.select_group(units, "chrA")
        assert [u["unit_id"] for u in kept] == ["a", "c"]

    def test_none_group_is_a_no_op(self):
        units = [{"unit_id": "a", "group": "chrA"}]
        assert mco.select_group(units, None) is units


class TestScriptEntryPoint:
    """Drives the file as a SUBPROCESS, the way a workflow task actually would.

    Everything else in this suite calls ``mco.main(...)`` on the imported
    module object, which never executes the ``if __name__ == "__main__":``
    guard. That gap is exactly how the guard was once dropped by an edit that
    replaced ``main``'s body and swallowed the block after it -- every
    direct-import test still passed, and ``--help`` silently exited 0 having
    done nothing. These tests exist specifically to make that regress loudly.
    """

    SCRIPT = REPO / "util" / "misc" / "merge_chunk_outputs.py"

    def test_help_runs_and_describes_the_script(self):
        r = subprocess.run(
            [sys.executable, str(self.SCRIPT), "--help"],
            capture_output=True, text=True, timeout=30,
        )
        assert r.returncode == 0
        assert "Stage 6" in r.stdout
        assert "--manifest" in r.stdout

    def test_no_args_exits_nonzero_naming_what_is_missing(self):
        r = subprocess.run(
            [sys.executable, str(self.SCRIPT)],
            capture_output=True, text=True, timeout=30,
        )
        assert r.returncode != 0
        assert "--manifest and --outdir are required" in r.stderr

    def test_a_real_merge_runs_end_to_end_as_a_subprocess(
        self, two_group_units, tmp_path
    ):
        manifest = tmp_path / "manifest.json"
        # write_manifest is tested elsewhere; build the file directly here so
        # this test does not depend on the module under a DIFFERENT loaded
        # name (subprocess) matching the one imported in-process above.
        manifest.write_text(
            json.dumps(
                {
                    "units": [
                        {
                            "unit_id": u["unit_id"],
                            "quant_prefix": u["quant_prefix"],
                            "offset": u["offset"],
                        }
                        for u in two_group_units
                    ]
                }
            )
        )
        outdir = tmp_path / "out"
        r = subprocess.run(
            [
                sys.executable, str(self.SCRIPT),
                "--manifest", str(manifest),
                "--outdir", str(outdir),
            ],
            capture_output=True, text=True, timeout=30,
        )
        assert r.returncode == 0, r.stderr
        assert (outdir / "merged" / "chunked.quant.expr").exists()
        assert (outdir / "merged" / "chunked.quant.tracking.gz").exists()


class TestWriteManifest:
    def test_default_grouping_strips_orientation_and_index(self, tmp_path):
        units = [
            {"unit_id": "chrUn_KI270302v1_00_plus", "quant_prefix": "/x", "offset": 0},
            {"unit_id": "chr1_05_minus", "quant_prefix": "/y", "offset": 0},
        ]
        p = tmp_path / "m.json"
        mco.write_manifest(str(p), units)
        doc = json.loads(p.read_text())
        groups = {u["unit_id"]: u["group"] for u in doc["units"]}
        assert groups["chrUn_KI270302v1_00_plus"] == "chrUn_KI270302v1"
        assert groups["chr1_05_minus"] == "chr1"

    def test_round_trips_through_load_manifest(self, two_group_units, tmp_path):
        p = tmp_path / "m.json"
        mco.write_manifest(str(p), two_group_units)
        reloaded = mco.load_manifest(str(p))
        assert [u["unit_id"] for u in reloaded] == [
            u["unit_id"] for u in two_group_units
        ]
        assert [u["offset"] for u in reloaded] == [
            u["offset"] for u in two_group_units
        ]


class TestCliMatchesDirectCall:
    @pytest.mark.parametrize("discovery", [False, True])
    def test_cli_round_trip_byte_identical_to_direct_call(
        self, two_group_units, tmp_path, discovery
    ):
        manifest = tmp_path / "manifest.json"
        mco.write_manifest(str(manifest), two_group_units)

        cli_outdir = tmp_path / "via_cli"
        direct_outdir = tmp_path / "via_direct"

        rc = mco.main(
            [
                "--manifest", str(manifest),
                "--outdir", str(cli_outdir),
            ]
            + (["--discovery"] if discovery else [])
        )
        assert rc == 0

        direct = ChunkedRun.merge_and_translate(
            str(direct_outdir), two_group_units, discovery=discovery
        )

        cli_expr = cli_outdir / "merged" / "chunked.quant.expr"
        cli_track = cli_outdir / "merged" / "chunked.quant.tracking.gz"
        assert _read_text(str(cli_expr)) == _read_text(direct["quant_expr"])
        assert _read_text(str(cli_track)) == _read_text(direct["quant_tracking"])

    def test_discovery_namespaces_ids_and_plain_quant_does_not(
        self, two_group_units, tmp_path
    ):
        bulk = tmp_path / "bulk"
        denovo = tmp_path / "denovo"

        manifest = tmp_path / "manifest.json"
        mco.write_manifest(str(manifest), two_group_units)

        mco.main(["--manifest", str(manifest), "--outdir", str(bulk)])
        mco.main(
            ["--manifest", str(manifest), "--outdir", str(denovo), "--discovery"]
        )

        bulk_expr = _read_text(str(bulk / "merged" / "chunked.quant.expr"))
        denovo_expr = _read_text(str(denovo / "merged" / "chunked.quant.expr"))

        sep = ChunkedRun.NAMESPACE_SEP
        assert "\tt1\t" in bulk_expr and "{}t1\t".format(sep) not in bulk_expr
        assert "chrA_00_plus{}t1\t".format(sep) in denovo_expr

    def test_result_json_is_written_and_matches_the_direct_dict(
        self, two_group_units, tmp_path
    ):
        manifest = tmp_path / "manifest.json"
        mco.write_manifest(str(manifest), two_group_units)
        outdir = tmp_path / "out"
        result_path = tmp_path / "result.json"

        mco.main(
            [
                "--manifest", str(manifest),
                "--outdir", str(outdir),
                "--result", str(result_path),
            ]
        )
        result = json.loads(result_path.read_text())
        assert result["expr_rows"] == len(two_group_units) * len(MODELS)
        assert result["coordinates_translated"] is True  # nonzero offsets


class TestHeaderMismatchRefusesBeforeOutput:
    def test_expr_header_mismatch_raises_and_writes_nothing(
        self, two_group_units, tmp_path
    ):
        # Corrupt the second unit's expr header so it disagrees with the first.
        bad_prefix = two_group_units[1]["quant_prefix"]
        text = Path(bad_prefix + ".quant.expr").read_text()
        lines = text.splitlines()
        lines[0] = lines[0] + "\textra_column"
        Path(bad_prefix + ".quant.expr").write_text("\n".join(lines) + "\n")

        manifest = tmp_path / "manifest.json"
        mco.write_manifest(str(manifest), two_group_units)
        outdir = tmp_path / "out"

        with pytest.raises(ChunkedRun.PipelineError, match="header differs"):
            mco.main(["--manifest", str(manifest), "--outdir", str(outdir)])

        # merge_and_translate mkdirs merged/ before reading any unit, so the
        # directory exists -- but the refusal must land before either output
        # file inside it does.
        merged_dir = outdir / "merged"
        assert not (merged_dir / "chunked.quant.expr").exists()
        assert not (merged_dir / "chunked.quant.tracking.gz").exists()


class TestGroupedMergeIdentities:
    """The module docstring's two claims, exercised through the real CLI."""

    def _run_global(self, units, tmp_path):
        manifest = tmp_path / "global_manifest.json"
        mco.write_manifest(str(manifest), units)
        outdir = tmp_path / "global"
        mco.main(["--manifest", str(manifest), "--outdir", str(outdir)])
        return outdir / "merged"

    def _run_groups(self, units, tmp_path):
        manifest = tmp_path / "grouped_manifest.json"
        mco.write_manifest(str(manifest), units)
        groups = sorted({u["unit_id"].split("_")[0] for u in units})
        outputs = {}
        for g in groups:
            outdir = tmp_path / "group_{}".format(g)
            result = tmp_path / "result_{}.json".format(g)
            mco.main(
                [
                    "--manifest", str(manifest),
                    "--outdir", str(outdir),
                    "--group", g,
                    "--result", str(result),
                ]
            )
            outputs[g] = {
                "merged_dir": outdir / "merged",
                "result": result,
            }
        return outputs

    def test_tracking_is_byte_identical_via_heapq_merge_not_concatenation(
        self, two_group_units, tmp_path
    ):
        global_merged = self._run_global(two_group_units, tmp_path)
        grouped = self._run_groups(two_group_units, tmp_path)

        global_rows = _read_text(
            str(global_merged / "chunked.quant.tracking.gz")
        ).splitlines()
        header, global_body = global_rows[0], global_rows[1:]

        group_bodies = []
        for g in sorted(grouped):
            lines = _read_text(
                str(grouped[g]["merged_dir"] / "chunked.quant.tracking.gz")
            ).splitlines()
            assert lines[0] == header
            group_bodies.append(lines[1:])

        # Concatenation must NOT reproduce the global order (proves the test
        # fixture actually spans a real interleave, not a degenerate case).
        concatenated = [row for body in group_bodies for row in body]
        assert concatenated != global_body

        # heapq.merge over the pre-sorted group outputs MUST reproduce it.
        cols = header.split("\t")
        key_idx = [cols.index(c) for c in ("read_name", "transcript_id", "gene_id")]
        key = lambda line: [line.split("\t")[i] for i in key_idx]  # noqa: E731
        for body in group_bodies:
            assert body == sorted(body, key=key)
        merged_via_heapq = list(heapq.merge(*group_bodies, key=key))
        assert merged_via_heapq == global_body

    def test_expr_is_byte_identical_via_combine_expr_not_concatenation(
        self, two_group_units, tmp_path
    ):
        global_merged = self._run_global(two_group_units, tmp_path)
        grouped = self._run_groups(two_group_units, tmp_path)

        global_expr = _read_text(str(global_merged / "chunked.quant.expr"))

        combined_path = tmp_path / "combined.quant.expr"
        combined_result_path = tmp_path / "combined_result.json"
        rc = mco.main(
            [
                "--combine-expr",
                "--out", str(combined_path),
                "--result", str(combined_result_path),
            ]
            + [
                arg
                for g in sorted(grouped)
                for arg in ("--group-result", str(grouped[g]["result"]))
            ]
        )
        assert rc == 0

        # Naive concatenation of the grouped tables must NOT match: each
        # group's TPM was rebased over only that group's all_reads.
        naive_concat = "".join(
            _read_text(str(grouped[g]["merged_dir"] / "chunked.quant.expr"))
            for g in sorted(grouped)
        )
        assert naive_concat != global_expr

        combined = _read_text(str(combined_path))

        def rows_by_key(text):
            lines = text.splitlines()
            header = lines[0].split("\t")
            gi, ti, tpmi = header.index("gene_id"), header.index("transcript_id"), header.index("TPM")
            return {
                (r[gi], r[ti]): r[tpmi]
                for r in (line.split("\t") for line in lines[1:])
            }

        want = rows_by_key(global_expr)
        got = rows_by_key(combined)
        assert set(got) == set(want)
        assert got == want

        combined_result = json.loads(combined_result_path.read_text())
        assert combined_result["combined_from_groups"] == len(grouped)
        assert combined_result["expr_rows"] == len(two_group_units) * len(MODELS)

    def test_combine_expr_requires_at_least_one_group_result(self):
        with pytest.raises(SystemExit):
            mco.main(["--combine-expr", "--out", "/tmp/x.expr"])

    def test_combine_expr_requires_out(self, tmp_path):
        fake = tmp_path / "r.json"
        fake.write_text(json.dumps({"quant_expr": "/x", "merged_scope_total_all_reads": 1.0}))
        with pytest.raises(SystemExit):
            mco.main(["--combine-expr", "--group-result", str(fake)])


class TestRebaseTpm:
    """The arithmetic ``merge_and_translate`` and ``combine_grouped_expr`` share."""

    def test_matches_the_original_inline_formula(self):
        rows = [["10"], ["30"], ["0"]]
        header = ["all_reads", "TPM"]
        rows = [[r[0], "0.000"] for r in rows]
        total = ChunkedRun.rebase_tpm(rows, header)
        assert total == 40.0
        assert rows[0][1] == "{:.3f}".format(10 / 40 * 1e6)
        assert rows[1][1] == "{:.3f}".format(30 / 40 * 1e6)
        assert rows[2][1] == "0.000"

    def test_total_override_is_used_instead_of_summing(self):
        rows = [["10", "0.000"]]
        header = ["all_reads", "TPM"]
        total = ChunkedRun.rebase_tpm(rows, header, total_override=1000.0)
        assert total == 1000.0
        assert rows[0][1] == "{:.3f}".format(10 / 1000 * 1e6)

    def test_zero_total_yields_zero_tpm_not_a_division_error(self):
        rows = [["0", "0.000"]]
        header = ["all_reads", "TPM"]
        ChunkedRun.rebase_tpm(rows, header, total_override=0.0)
        assert rows[0][1] == "0.000"
