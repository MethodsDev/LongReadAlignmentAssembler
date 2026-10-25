# --pasa_mode and --preassembled_input

## --pasa_mode: chaining fragments into longer transcripts

LRAA normally only collapses reads into the transcript they fully represent; a read
never extends another. That suits full-length long reads, but not inputs that are
fragments of longer transcripts: ESTs, partial cDNAs, or staggered reads that each
cover part of a transcript. `--pasa_mode` chains such fragments into longer models,
as PASA's alignment assembler does:

- Two spliced alignments are chained when they are compatible, overlap, and **share
  at least one intron**. Fragments that overlap only within an exon are not chained,
  and unspliced alignments never chain (they are assembled in LRAA's separate
  single-exon pass, as usual).
- A terminal exon running at most 20 bp into an intron, past a splice site the
  terminal exon also covers, is trimmed back to that splice site (PASA's fuzz
  distance), so short overhangs do not break the intron chain.

Everything else runs as in a standard LRAA run: TSS and polyA inference, EM
quantification, isoform filtering, transcript ends (trimmed from the read end
distribution) and the outputs (GTF, tracking, quant.expr). Not every input
alignment has to end up in an assembled transcript. The settings are
`LRAA_Globals.PASA_MODE_CONFIG`; explicit flags and `--config_update` override them.

TSS inference stays on: stacked 5' read ends are taken as evidence for a TSS, and a
TSS ends a chain there. When the 5' ends of fragments are not expected to mark TSSs,
add `--no_infer_TSS`.

## --preassembled_input: one alignment is enough

By default an intron needs 2 supporting reads to enter the splice graph
(`min_intron_support`), and a novel isoform needs 2 reads (`--min_reads_novel`).
That is the right default for long-read RNA-seq, where a single read can carry an
alignment error.

When the inputs are preassembled transcripts (Trinity or StringTie assemblies, for
example), each alignment already summarizes many reads, and one is evidence enough.
`--preassembled_input` sets both minimums to 1. It works with or without
`--pasa_mode`. An explicit `--min_reads_novel` or `--config_update` still wins.

**Do not use `--preassembled_input` for long-read RNA-seq.**

| input | flags |
|---|---|
| full-length long reads | (none) |
| fragmented / staggered long reads | `--pasa_mode` |
| preassembled transcripts, ESTs, cDNAs | `--pasa_mode --preassembled_input` |

## Chunked mode

Both switches are config presets. In chunked mode (the default), the driver writes the
resolved settings to `worker_config_overrides.<digest>.json` in the chunk work
directory, and every chunk worker applies that file last, so the workers assemble
with the same settings as an unchunked run.

## Tests

- `pylib/test_chain_shared_intron.py`, `pylib/test_preassembled_input.py`
- `testing/pasa_mode`: synthetic staggered reads (`build_fixture.py`), checked for the
  expected intron chains (`make test`).
