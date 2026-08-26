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

## The sibling comparison: strand-first against strandless chunking

Everything above compares chunked against whole-contig. `--strandless_chunks` asks a
different question with the same shape — does moving the orientation split from a whole-BAM
serial phase into the per-chunk parallel phase change the answer — and it is answered by
`pylib/StrandlessParity.py`, driven from `testing/strandless_parity/`. It reuses this
evaluation's machinery rather than restating it: `ChunkedRun.severed_read_names`,
`ChunkedRun.read_tsv` and `ChunkedRun.resolve_tracking` are what it reads the arms with, so
the two comparisons cannot drift into disagreeing about what an arm consumed.

Two things differ from the chunked-vs-whole comparison and both change what the bar can be.

**There is a cheaper and stronger check available first.** Strand assignment is per-read
local — `split_bam_by_strand` decides from the read alone — so a read's strand must be
identical computed over the whole contig or inside the chunk holding it. That is checkable
without quantifying anything (`strand-invariant`), and if it fails no downstream number can
be made to agree. Measured on chr1 (HG002 PacBio Kinnex, 1,909,780 records, 1,064,804
retained, 24 union cuts at the 10 Mb default), in both modes:

| mode | records compared | agree | flag rewritten by the split |
|---|---|---|---|
| aligner flag (what the pipeline runs today) | 1,064,804 | 1,064,804 | 0 |
| `--infer_read_orient` | 1,064,804 | 1,064,804 | 12,471 whole-genome, the same 12,471 in-chunk |

**Quote the second row, not the first.** Under plain aligner flags the split is a filter and a
partition on `read.is_forward`, so a record's strand is a property of the record and cannot
disagree with itself no matter where the split ran: that row is a wiring check, and 100%
agreement in it is not evidence for the design. Inference is the only regime in which a read's
strand can differ from its alignment flag at all, and therefore the only regime in which the
invariant could have been false. There the split rewrote 12,471 flags genome-wide, rewrote the
same 12,471 inside the chunks, and no read changed strand. Nothing was severed at that
geometry, so the two arms' read sets are identical too, and every disagreement the check could
have reported would have been a real one.

**The arms may legitimately consume different reads.** Cuts are chosen from different
evidence — one strand's coverage against both strands' union — so their POSITIONS can differ
and different positions sever different reads. On the minigenome at 0.5 Mb the strandless arm
severs 6 where strand-first severs 5, and the extra read is named rather than counted:
`m160629_191727_..._s1_p0/69114/ccs`, reverse-oriented, spanning 928,251-1,020,755, severed
by the shared cut at 1,014,000 that the minus strand's own cut at 1,046,300 misses. It reaches
neither arm's `quant.tracking`, so the difference costs the comparison nothing and the merged
tables come out identical anyway. The gate therefore does not demand byte identity
unconditionally: it demands that every differing row be attributable to a read in the
symmetric difference of the two severed sets, and names any row that is not.

`quant.tracking` differs between the arms in `mp_id` alone, which is a per-process MultiPath
counter rather than data (`ChunkedRun.merge_and_translate` says the same). The gate tests row
ORDER separately from row CONTENT so that an ordering difference is proven — equal row
multisets, which no content difference can satisfy — rather than asserted; on this corpus the
rows are in the same order, so "ordering" would have been the wrong description and the gate
says `mp_id` instead.

**The ordering mistake, measured rather than warned about.** `retained_for_extraction`
filters on the RAW alignment flag and `split_bam_by_strand` rewrites that flag, so extracting
`chr1+:...` from a raw BAM keeps every read the aligner called forward and then treats all of
them as forward — including those whose splice dinucleotides say otherwise, which no later
stage revisits. `ordering-cost` commits the mistake on purpose and names what it costs:
56 of 245 emitted records on `minigenome+:1014001-1502500`, and 21 of 25,757 on
`chr1+:9919601-20000000`. The two rates differ by 300x and the reason is the substrate, not
the geometry: the minigenome's PacBio BAM is not oriented with transcription (1,076 of 2,266
records are flipped by the split, 47.5%) while HG002 Kinnex largely is (12,471 of 1,064,804,
1.2%). Quoting either figure alone would misstate the risk — the minigenome overstates how
visible the mistake is, chr1 understates it — so both are here, and neither is zero.

Unit coverage is `pylib/test_strandless_parity.py`, whose fixture plants splice
dinucleotides that CONTRADICT the aligner flag: without a read the split rewrites, the
invariant would pass while testing nothing. It also runs each check against a deliberately
broken arrangement, because an instrument that cannot register a fault is worthless.

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

Strandless chunking, `--strandless_chunks`: cut the RAW bam, extract one chunk per
interval holding both orientations, and run the orientation split inside each chunk.
The chunked arm then does not run stage 1 at all, so it has no retained record count
to default `-N` to and requires one -- use `stage1_retained_records.total` from a
strand-first run over the same bam, which is what makes the two arms' RPM column
comparable. `--arm both` still runs stage 1, because the control IS the strand-split
whole bam.

```bash
cd <repo>/testing/single_contig
python3 ../../util/misc/run_chunked_quant_pipeline.py \
    --bam pacbio.PBLR.bam --genome_fa ref_genome.fa --gtf ref_annot.gtf \
    --output_dir /tmp/parity_strandless --arm chunked --cpu_budget 6 --HiFi \
    --strandless_chunks -N 2266 \
    --approx_MB_per_cut 0.5 --approx_MB_per_cut_wiggle_window 0.1
```

An output directory serves one mode or the other: `merged/`, `timing.json` and
`outputs.json` are at fixed paths and stage 6 is not checkpointed, so the second mode
would overwrite the first's results while every per-stage sentinel reported new work.
The run refuses instead. `--dry_run` stops after cut selection -- which is the plan --
and prints the units it would run: 25 intervals and 50 quant units on chr1 at 10 Mb,
against 50 chunks and 50 quant units strand-first.

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

Synthetic 3.4 Mb contig, cuts forced to 0.5 Mb — 6 cuts per strand, 6 reads severed and
excluded from the control, 555 transcripts: all columns equal, 544 tracking rows identical,
409 distinct reads assigned in both arms. No multi-gene components exist in this corpus, so
this run demonstrates the subtraction and the plumbing only.

Two of those figures were re-measured after the 2026-08-15 run, and neither is a parity
finding — both arms moved together each time, with rows identical in order, every numeric
column delta 0.0, and severed symmetric difference 0:

- **severed reads 5 → 6**, re-pinned in `669459f4`, which verified against an unmodified
  `devel` worktree that the figure had been stale since before that branch existed.
- **tracking rows 541 → 544, distinct reads 408 → 409**, bisected to `c940609b`, "Compute
  percent identity over the columns the mismatch tag counts". That commit repaired an
  identity calculation which counted gaps in the numerator but not the denominator, leaving
  the value unbounded (it measured −13900%) and discarding alignments at true identities of
  99.0–99.9%. This corpus is PacBio and the gate runs `--HiFi`, i.e. a high identity floor,
  which is exactly where those reads were being discarded: one more read is now assigned,
  carrying three more rows. The `+3` is this corpus's measured outcome rather than a general
  rule — admitting a read changes EM and can change the survivor set, so row counts need not
  move monotonically with retention.

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

## The single-cell three-input shape (2026-08-26)

Everything above concerns the bulk shape, where one BAM serves every role. The
cluster-guided single-cell path splits that into THREE inputs, and the parity
comparison has to match all three or it measures the wrong thing.

    --bam              this cluster's full reads        pass-2 assignment
    --bam_for_sg       shared merged normalized BAM     splice graph only
    --bam_for_priors   this cluster's normalized reads  pass-1 theta

### What chunking removes from each

Cut selection severs reads that straddle a chunk boundary, and each input loses its
own set. MEASURED on the by_chunk qonly fixture, cluster 00, chr19 2 Mb, plan of 10
chunks emitted over the pre-partition superset -- comparing each per-chunk slice
union against the whole input:

    --bam              5 of 27,142
    priors slice       5 of 13,320   (the same 5 ids)
    sg slice         156 of 95,064

### Result

Chunked merged output against an unchunked run, three matchings of increasing
strictness:

    unchunked, only --bam pruned          0/654 rows differ    INVALID
    unchunked, --bam + priors pruned      0/654 rows differ
    unchunked, all three matched          0/654 rows differ    sum 24,083.3 both

Exact on reads AND TPM, 654 transcripts, no set difference.

SCOPE, and it is narrow: ONE cluster, ONE 2 Mb contig, 10 chunks, quant-only. This
does NOT establish that by_chunk is generally quantification-neutral, and it does
NOT substitute for a whole-genome correctness check. Broader parity is UNMEASURED.

### Why matching only `--bam` is invalid

Pruning only `--bam` reports exact parity while being invalid: the unchunked arm's
graph still sees 156 records and its theta 5 records the chunked arm never had. It
gave the right answer for the wrong reason. The code tracks
`sg_dropped_read_names` separately from `dropped_read_names` for exactly this
shared-cuts case -- see `verify_severed_accounting` in `pylib/ChunkedRun.py`. The
reliable construction is to restrict `--bam_for_sg` and `--bam_for_priors` to the
UNION of read names present across the per-chunk `chunk.sg.bam` / `chunk.priors.bam`
slices, which is by definition what the chunked arm saw.

`--arm both` cannot produce this comparison. Its baseline self-normalizes the whole
contig, so with `--bam_for_sg`/`--bam_for_priors` supplied it would build the graph
and estimate theta from evidence the chunked arm never saw; the pipeline refuses the
combination rather than reporting an incomparable pair. The subtraction must be
applied to the externally supplied files instead.
