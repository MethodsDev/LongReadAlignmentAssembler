# Regression test: streaming quant-only canonical-path collision

A ~2 s, 200 KB test for the LRAA v0.31.0 bug that aborted both PBMC PacBio
cluster-guided quant-only runs on 2026-09-07.

    fixture/                                        40 kb locus, 200 KB total
    test_streaming_quant_canonical_path_collision.py  pytest, for the LRAA repo
    run_test.sh                                     no-pytest execution test, for an IMAGE
    build_fixture.py                                how fixture/ was derived

## Verified red/green

| | `run_test.sh` | pytest |
|---|---|---|
| `0.31.0-03b9c9a` (buggy) | FAIL rc=1 | 4 failed in 3.75 s |
| `0.33.0-7e37488` (fixed) | PASS, 30 rows | 4 passed in 4.13 s |

Both were run on 2026-09-07. A regression test that has never been observed to
fail is not a regression test, so the 0.31.0 leg is the point of the table.

## Running it

    ./run_test.sh /path/to/lraa-core_<tag>.sif        # or a docker:// uri, or a docker tag
    LRAA_HOME=/path/to/LRAA pytest test_streaming_quant_canonical_path_collision.py

`run_test.sh` uses `mktemp -d`, so export `TMPDIR` to somewhere exec-permitted if
`/tmp` is mounted `noexec` (it is on the Broad methods boxes, and Apptainer's
proot helpers exec out of the temp dir).

## What it covers

The failure needed BOTH conditions, which is why nothing else in the pipeline
caught it:

1. **streaming** quant (`--stream_reads`) -- `StreamingQuant.AssignmentTable.build`
   is the only caller of the guard that raised;
2. an **uncollapsed** annotation -- the raw `init_gtf` from initial discovery,
   never through `merge_LRAA_GTFs.py` /
   `collapse_LRAA_GTF_by_splice_pattern.py`.

The cluster-guided final quant runs the same streaming code over the same reads
and the same chunk geometry and passed, because it quantifies the *collapsed*
final GTF. The uncollapsed init GTF is nevertheless a documented, first-class
reuse path (`precomputed_init_gtf` in `LRAA-singlecell.wdl`), so it needs
coverage of its own. `REQUIRED_FLAGS` in the pytest file pins condition 1 so a
future edit cannot quietly drop the flags and leave a test that always passes.

Four assertions, in increasing strength:

- the v0.31.0 error string is absent;
- LRAA exits 0 (StreamingQuant could raise something new);
- `comp-1119`, the component whose multipaths collided, appears in
  `quant.expr` -- a "fix" that dropped the offending path would pass the first
  two and still be wrong;
- `quant.tracking.gz` is readable and has data rows -- `AssignmentTable` builds
  those rows, so tracking is the output the guard sits directly upstream of.

## Fixture provenance

`chr21:6,280,000-6,320,000` from `PBMC.BASIC.refguided`, chunk `chr21_04`, plus
strand, translated to offset 0 (subtract 6,279,999) so the contig is 40 kb
rather than 6.7 Mb. It holds the whole colliding component, `comp-1119`
(originally `chr21:6,286,342-6,313,221`), with ~13 kb of margin each side. Only
reads falling entirely inside the window were kept, so the component's splice
graph is the one the failure was built from.

Translation shifts the literal path string in the error message, so the test
asserts on the failure MODE, never on those coordinates.

`build_fixture.py` regenerates `fixture/` from the full 69 MB chunk, if that is
ever needed. It requires the original chunk, which is preserved at
`../repro/chr21_04_canonical_path_collision/`.
