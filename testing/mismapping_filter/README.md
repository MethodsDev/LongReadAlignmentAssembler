# Alignment-mismapping filter — minigenome regression tests

Execution regression tests for the whole-genome alignment-mismapping filter
added in v0.40.0 (`pylib/AlignmentMismappingFilter.py`; see
`docs/alignment_mismapping_filter.md`).

## What these are

Four self-contained **minigenomes** (a genome fasta + aligned reads + optional
truth GTF, coordinates rebased to a small local contig). Each contains a real,
highly-expressed **source gene** together with the alignment/strand-mismapping
**artifacts derived from it** that the filter is meant to remove:

| example | class | what it is |
|---|---|---|
| `s1_mirror`, `s2_mirror` | `s` | wrong-strand **near-mirrors**: opposite-strand models whose splice sites sit on the source gene's exon boundaries |
| `x1_runon`, `x2_runon` | `x` | **run-on / near-identical copies**: opposite-strand models whose cDNA is ≥99% identical over most of their length to the source gene |

Each minigenome was built to include the source transcript the artifacts were
mismapped from — the filter's sequence detector needs that higher-expressed
partner present to recognize the copy. (`x2_runon` has no annotated Sequin gene
in its window, so it ships without a truth GTF; its reconstructed source-strand
models are the reference, and the test is truth-independent regardless.)

## How the test works

`test_mismapping_filter.py` runs the checked-out `../../LRAA` de novo on each
minigenome TWICE — once with `--no_filter_mismappings` (artifacts present) and
once with the filter on — and asserts:

1. the artifact reproduces with the filter OFF (≥1 opposite-strand model);
2. the filter ON reduces the opposite-strand (artifact) model count;
3. the removal log removed at least the expected number of artifacts;
4. the real **source-strand** models are unchanged (the filter must only remove);
5. the filtered `quant.expr` TPM sums to 1e6 (row-drop + renormalize, no requant).

```
make test        # runs LRAA 8x (~15-20 min); this is a SLOW target
make clean
```

## De-identification

The reads are derived from an ONT cDNA Sequin spike-in sample but are
**de-identified** so they cannot be traced back to the source run:

- every read QNAME is renamed to an anonymous sequential id (`r1`, `r2`, …);
- the BAM header carries only `@HD` and `@SQ` — no `@PG`/`@RG`/`@CO` (which had
  carried the run accession and file paths);
- only the `ts` (splice-motif strand) and `NM` (edit distance) aux tags are
  retained — `ts` is functionally required to reproduce the strand-mismapping
  behavior, and neither tag carries sample metadata; all other tags are stripped;
- read sequences and base qualities are preserved unchanged.

The genome contig (`chrIS`) and any gene ids are the public Sequin reference, not
sample-identifying.
