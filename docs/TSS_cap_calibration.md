# TSS cap signature: library calibration and the cap rules

A read whose 5' end is a capped RNA end carries one more untemplated G than the library's
uncapped reads: reverse transcriptase adds a C opposite the cap during template switching.
A 5' truncation or degradation product starts at a pile-up that looks like a TSS but lacks
that extra G. LRAA uses the signature to prune such sites (`prune_uncapped_alt_TSS`, on by
default) and, in the single-cell workflow's cluster runs, to admit new TSSs only with it
(`require_cap_for_new_TSS`). Both need to know first whether the library carries the
signature at all, and its baseline G count.

## Read-level definitions

- **G run**: the run of G's in a read's 5' soft clip that touches the alignment, recorded by
  `Pretty_alignment` before the untemplated-G strip (a reverse-strand read: the C run at the
  head of its right clip).
- **Trusted 5' end**: the read's start is taken as the first base of its molecule only when
  its 5' clip is accounted for: nothing but 1-3 G's (Kinnex, Iso-Seq: the template-switch oligo
  already trimmed), or a run of at least 3 G's next to the alignment (ONT PCS109: adapter +
  the strand-switching primer, whose 3' GGG anneals to the C's reverse transcriptase adds).
  The >= 3 G rule marks the primer, not the cap. A trusted end is a real molecule end, not
  necessarily a TSS: away from peaks most are truncated or degraded molecules. Reads with an
  unexplained clip are skipped.
- **Informative**: the genomic base just 5' of the read start is not G (otherwise a capped
  read's extra G aligns to it and the read shows nothing).
- **Baseline**: the most common G run among trusted, informative, spliced read starts away
  from every TSS peak, pooled over the library (ties to the smaller run); an int
  `TSS_cap_G_baseline` in the config is used as given. Measured: 0 for sc Kinnex 10x, 1 for
  bulk PacBio Kinnex (likely an extra G in the bulk Kinnex PCR adapters; the library's uncapped
  SIRV spike-ins carry it too), 3 for ONT PCS109 (the primer's GGG). A read is **capped** when
  its G run exceeds the baseline: the cap adds one G (reads with two extra are rarer).

## The library decision (pre-pass, `pylib/TSS_cap_calibration.py`)

Before any chunking, the bam is scanned once, one worker per contig, strands kept apart. TSS
peaks are called as a graph build calls them (`aggregate_sites_within_window`,
`filter_non_peaky_positions`), and each peak's capped fraction is taken over the spliced
informative read starts within 50 nt of it. Peaks with at least
`TSS_cap_prepass_min_peak_reads` (20) such reads are judged.

With the signature, the judged peaks' capped fractions are **bimodal**: an uncapped mode at
the transcript-body rate (truncation and degradation pile-ups) and a capped mode (promoters).
Without it (SIRVs, simulated reads, 5'-trimmed reads) there is one mode at 0. The fractions are
split in two where they separate best (Otsu's threshold: the cut maximising the between-group
variance), and:

| decision | when | effect |
|---|---|---|
| on | the capped group holds >= `TSS_cap_prepass_min_capped_mode_frac` (10%) of judged peaks and >= `TSS_cap_prepass_min_capped_mode_peaks` (50) peaks, and its median exceeds the uncapped group's by >= `TSS_cap_prepass_min_mode_separation` (0.20) | `TSS_cap_signature` "on", `TSS_cap_G_baseline` fixed |
| off | otherwise | `TSS_cap_signature` "off": no cap rules |
| undetermined | fewer than `TSS_cap_prepass_min_peaks` (50) judged peaks, or a clear but thin capped mode | nothing set: `TSS_cap_signature` stays "auto", each graph build decides |

Measured on ten libraries (LRAA paper analyses, `__TSS_5prime_Completeness/cap_calibration`):
uncapped / capped mode medians FTD nuclei 1.4 / 48%, PBMCs 1.5 / 65%, BT474 and K562 PacBio
2 / 77%, MCF7 and A549 ONT 9 / 39-44%; 91-97% of capped-mode peaks at a FANTOM5 CAGE peak
against 15% of uncapped-mode ones; MORFs, SIRVs and an Arabidopsis simulation off or
undetermined.

## Where it runs

- **Inside LRAA**: at the start of a discovery run whose `TSS_cap_signature` is "auto" (with
  `TSS_cap_prepass`), after every config source and before chunking; the result reaches every
  chunk worker with the forwarded config. Not in quant-only runs.
- **Standalone**: `util/calibrate_TSS_cap_signature.py --bam --genome --output_prefix`
  (or `--scan_table` to re-decide from a saved scan). It writes
  `<prefix>.TSS_cap_calibration.config.json` for `--config_update`: the decision's keys and
  `TSS_cap_prepass` false, so a run given it never calibrates again on its own reads.
- **WDL**: `LRAA.wdl` runs the standalone step once on the whole bam before any sharding
  (task `calibrate_TSS_cap_signature`, input `run_TSS_cap_calibration`, default true) and
  merges the result into every shard's `config_update_json`; skipped under `quant_only` or
  when the caller's config sets `TSS_cap_signature`. Outputs `capCalibrationJson` (the keys)
  and `capCalibrationReport`. `LRAA-singlecell.wdl` hands the init run's calibration to every
  per-cluster discovery run (`internal_cap_calibration_json`), so clusters inherit the
  library's decision; it surfaces `cap_calibration_json` and `cap_calibration_report`.

Outputs, either way: `<prefix>.TSS_cap_calibration.json` (decision, config keys, statistics),
`.hist.tsv` (per-peak capped-fraction histogram), `.svg` (its plot, no plotting library
needed) and `.scan.tsv.gz` (the counts behind the decision). Execution test:
`testing/TSS_cap_calibration`.

## The per-site rules (unchanged by the calibration)

A site is judged with at least `min_TSS_cap_informative_reads` (5) informative reads. In
discovery, a TSS downstream of a capped site on any splice-graph path is pruned when under
`min_TSS_cap_frac` (10%) of its reads are capped. With `require_cap_for_new_TSS` (cluster
runs), a TSS not within 50 nt of a guide 5' end with `TSS "True"` needs at least
`min_new_TSS_capped_reads` (2) capped reads and 10%. The 10% threshold sits at or below the
valley between the two modes in every library measured (valleys 0.15-0.35).

## Possible extension: a fitted mixture model and per-site posteriors

Not implemented; recorded for later (2026-10-29). The current rule is deliberately simple.

The per-peak counts (informative reads n, capped reads k) could be fitted with a
**two-component beta-binomial mixture** by EM: the binomial part accounts for each peak's
depth (20 reads at 15% is weak evidence, 500 reads is not), the beta part for real spread
within a mode (promoters run 30-80% capped). It would give:

1. **A model-based library decision**, replacing the Otsu split: two components preferred
   over one (BIC), with effect-size criteria on the fitted parameters -- the capped
   component's weight >= 10%, its mean >= 0.20 above the uncapped one's, and separation such
   as Ashman's D >= 2. BIC alone is not enough: with thousands of peaks any heterogeneity
   favours two components (the MORFs' sparse capped-looking tail would pass). A parametric
   bootstrap would put confidence intervals on the weights for the report.
2. **Per-site posteriors**, replacing the fixed 10% threshold: for a site with n informative
   reads and k capped, the posterior probability that it belongs to the uncapped mode. A site
   with 5 reads and none capped is then uncertain and not pruned on thin evidence; one with
   200 reads at 4% is confidently uncapped; and the cut adapts to each library (ONT's uncapped
   mode reaches past 10%, so a fixed 10% keeps some ONT truncation sites).

The cost: per-site decisions would depend on a fitted model, needing its own validation (FTD,
PBMC and the benchmark chromosomes against the 10% rule) and a benchmark rerun. Step one would
be to fit the mixture to the ten libraries' scan tables and check it agrees with the current
decisions.
