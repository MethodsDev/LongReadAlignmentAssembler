# Chunked vs Whole-Contig Parity Evaluation

## Overview

Chunked quantification cuts a contig into segments, quantifies each independently, and
merges. The question that decides whether that is sound is whether it produces the same
numbers as quantifying the whole contig at once. The pipeline answers it directly:
`--arm both` runs the chunked arm and a whole-contig control over the same substrate, and
the two outputs can be compared transcript by transcript.

The orchestration lives in `pylib/ChunkedRun.py` as of v0.21.0. Two things reach it:
`LRAA --chunk`, which is how users run it, and `util/misc/run_chunked_quant_pipeline.py`,
which is the front end for this evaluation and the only route offering `--arm baseline`.
A whole-contig control is a measurement device, not a mode LRAA has reason to offer, which
is why the split falls there. Both routes produce byte-identical merged output, and
`testing/single_contig` asserts that on every run.

**The full evaluation is not part of routine testing** and is not expected to become part
of it: it needs a real chromosome with real annotation to be worth running, and `testing/`
has no room for one. What follows is enough to reconstruct it when it is needed: which
corpora make it meaningful, the exact invocations, the comparison, and the numbers a
healthy run produced on 2026-08-15 so a future run has something to differ from.

What *is* routine, since v0.21.0, is `make test` in `testing/single_contig`, which runs
both routes over the bundled minigenome with the cut size forced down to produce 14 chunks
and diffs their merged output. It executes every stage as a subprocess, so unlike the unit
tests it can catch a defect in the wiring between stages. That is the class this evaluation
found twice, and the class that let the v0.20.0 tracking-name break survive two releases
with nothing running chunking at all.

Unit coverage of the pieces lives in `util/misc/test_chunked_pipeline_{budget,checkpoints,parity}.py`,
`pylib/test_select_contig_cut_points.py`, `pylib/test_extract_contig_region_inputs.py`,
`pylib/test_gtf_tabix_index.py` and `pylib/test_chunked_entry_point.py`. None of those
executes a stage subprocess.

## What has to be true for the comparison to mean anything

A cut severs the reads that span it, and both neighbouring chunks drop them. That loss is
accepted; it is what makes chunking possible. The whole-contig control has no cut to lose them
at, so unless it drops them too it quantifies a record set the chunked arm never saw, and every
difference between the arms is confounded by exactly those reads.

`run_baseline` therefore prunes them, using the union of the cut selector's
`<prefix>.dropped_reads.txt` files, and quantifies `baseline/whole.parity.bam`. Confirm this
happened before trusting any comparison:

```
timing.json  ->  "baseline_excluded_severed_reads": N
baseline/whole.parity.bam exists when N > 0
```

With `N == 0` the pruned and unpruned inputs are the same file set and the control is the plain
merge, which is correct — there is nothing to subtract.

## Choosing a corpus

Two properties matter independently, and no single corpus available in-tree has both.

| property | why it matters | `testing/single_contig` | real chr20 |
|---|---|---|---|
| multi-gene read-sharing components | theta is normalized over a component, not a gene; splitting one changes every isoform fraction in it | **absent** | **present** |
| nonzero severing at the chosen cut geometry | exercises the parity subtraction | present when forced | absent at default geometry |

`testing/single_contig/pacbio.PBLR.bam` is 50 gene loci concatenated into one 3.4 Mb
pseudo-contig with a uniform 4,000 bp spacer and **no overlapping gene neighbours**, so no read
can be compatible with transcripts of two genes and every component is a single gene. It is
useful for exercising the plumbing and nothing else. Reads there that appear to span 60 kb are
minimap2 splicing across the spacers, not biology.

A real chromosome is required for the component half. The corpus used here:

```
/home/unix/bhaas/projects/LRAA_PAPER_Analyses/__LRAA_local_runs/indel_calibration/
    ont_chr20.bam        120,370 records, ONT
    chr20.fa             64 Mb
    chr20.annot.gtf      1,480 genes on chr20
```

## Running it

Real chromosome, default 10 Mb cut geometry:

```bash
cd /home/unix/bhaas/projects/LRAA_PAPER_Analyses/__LRAA_local_runs/indel_calibration
python3 <repo>/util/misc/run_chunked_quant_pipeline.py \
    --bam ont_chr20.bam --genome_fa chr20.fa --gtf chr20.annot.gtf \
    --contig chr20 --output_dir /tmp/chr20p --arm both --cpu_budget 10
```

Small synthetic contig, cuts forced because 3.4 Mb yields none at the default:

```bash
cd <repo>/testing/single_contig
python3 ../../util/misc/run_chunked_quant_pipeline.py \
    --bam pacbio.PBLR.bam --genome_fa ref_genome.fa --gtf ref_annot.gtf \
    --output_dir /tmp/parity --arm both --cpu_budget 6 --HiFi \
    --approx_MB_per_cut 0.5 --approx_MB_per_cut_wiggle_window 0.1
```

Outputs to compare are named in `<outdir>/outputs.json`:

```
baseline["quant_expr"]      baseline/baseline_quant.LRAA.quant-only.quant.expr
chunked["quant_expr"]       merged/chunked.quant.expr
   ... and the matching quant.tracking pair
```

## The comparison

Three levels, in increasing strength. `quant.expr` is TSV with `#` comment lines above the
header row; key on `transcript_id`.

**1. Inputs.** The arms must consume the same records before anything else is meaningful.

```bash
samtools view -c /tmp/chr20p/baseline/whole.parity.bam     # or whole.primary.bam when N==0
for b in /tmp/chr20p/chunks/*/chunk.bam; do samtools view -c "$b"; done | paste -sd+ | bc
```

Read names are the stronger form: extract `cut -f1 | sort -u` from the baseline input and from
the union of the chunk BAMs and `diff` them. Chunk BAMs carry rebased coordinates, so compare
names, not positions.

**2. Per-transcript numbers.** Load both `quant.expr` files into dicts keyed on
`transcript_id` and compare `uniq_reads`, `all_reads`, `isoform_fraction`,
`unique_gene_read_fraction`, `TPM` and `RPM_total_reads`, plus the structural columns
`exons`, `introns`, `splice_hash_code` and `gene_id`. Report the maximum absolute delta per
column rather than a pass/fail, so a small real divergence is distinguishable from a large one.

**3. Per-read assignments.** Key the two `quant.tracking` files on `(read_name, transcript_id)`
and compare the key sets and `frac_assigned`. This catches compensating errors that cancel in
the per-transcript totals.

**Confirm the hard case was present.** Before drawing any conclusion, count reads assigned to
transcripts of more than one `gene_id` in the baseline tracking. If that count is zero the run
says nothing about components, whatever the deltas look like.

## Results on 2026-08-15

Repository at `2cdc0da`, LRAA v0.19.2.

Real chr20, default geometry — 5 cuts per strand, **0 reads severed**, 5,822 transcripts:

| column | transcripts differing | max abs delta |
|---|---|---|
| `uniq_reads` | 0 | 0 |
| `all_reads` | 0 | 0 |
| `isoform_fraction` | 0 | 0 |
| `unique_gene_read_fraction` | 0 | 0 |
| `TPM` | 0 | 0 |
| `RPM_total_reads` | 0 | 0 |

Total `all_reads` 93,204.90 in both arms. Tracking identical: same `(read, transcript)` key set,
`frac_assigned` delta 0.

The hard case was present: **2,593 reads (2.8% of 93,212) compatible with transcripts of ≥2
genes, in 25 distinct multi-gene groupings**, the largest being
`MKKS^ENSG00000125863.20 + MKKS^ENSG00000285508.1` with 1,624 reads. The read-to-gene map was
identical in both arms.

Synthetic 3.4 Mb contig, cuts forced to 0.5 Mb — 6 cuts per strand, 5 reads severed and
excluded from the control, 555 transcripts: all columns equal, 541 tracking rows identical,
408 distinct reads assigned in both arms. No multi-gene components exist in this corpus, so
this run demonstrates the subtraction and the plumbing only.

## Why the arms agree, and the case neither run covers

Cut placement is annotation-aware: a boundary may not fall inside an annotated isoform, and
overlapping loci are collapsed into one indivisible island. Two genes that share a read overlap
or abut, so they occupy one island and no cut can be placed between them. Measured on the chr20
run: **0 of 25 multi-gene groupings had a cut inside their combined span**, with the four
largest at 10,401,009-10,434,222, 50,081,124-50,153,637, 5,950,949-6,040,053 and
646,615-675,800 against cuts at 9,985,800, 20,002,700, 30,000,000, 40,000,000 and 50,000,000.

So a cut splitting a read-sharing component appears to be unreachable by construction rather
than merely unobserved. Neither run tests it directly, and the selector actively prevents the
configuration from arising, so constructing it would mean defeating the annotation blocking on
purpose. Should that ever be attempted, the expected behaviour is still parity: a read bridging
two genes spans any cut placed between them, so it is severed and dropped from both arms
equally.

## Prior defects this evaluation found

Recorded because both were invisible to the unit tests and would be the first things to
re-check if a future run diverges.

- **The severed-alignment BAM was always empty** (`86de2f3`). Cut scoring and severed-read
  emission used predicates that disagreed on what a falsy strand means, and the CLI coerced
  argparse's `None` into `""`. Every run reported reads dropped at its cuts and wrote an empty
  BAM beside the report. Symptom to watch for: `spanning_alignments_dropped` nonzero in
  `cuts.json` while `samtools view -c` on the severed BAM returns 0.
- **The control kept the severed reads** (`2cdc0da`). `run_baseline` merged both strand BAMs
  and quantified all of them while its docstring claimed the merge was "exactly the union of
  every chunk's mini bam". On the synthetic corpus that was 2,266 records against the chunked
  arm's 2,261. Symptom to watch for: baseline record count exceeding the chunk sum by exactly
  the number of distinct names in `cuts/*.dropped_reads.txt`.
