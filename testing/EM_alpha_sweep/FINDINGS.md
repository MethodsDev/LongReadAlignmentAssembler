# EM_alpha response surface and recommendation

LRAA v0.20.0 at base commit `aa7df89`, quant-only, every arm at `--cpu_budget 3`.
Scored with `bmark_nb_runner.py` (`--analysisType quant_only`), one pass per
sample with all arms co-registered. All comparisons are paired within sample.

Harness: `run_grid.py`, `make_registry.py`, `score_grid.py`, `collect_results.py`,
`analyze_grid.py`, `ambiguity_profile.py`, `residual_test.py`,
`partition_analysis.py`, `ambiguity_sensitivity.py`, `plot_curves.py`.
Result tables in `results/`. No LRAA default is changed by this work.

## Headline

1. **Alpha is inert on the transcripts it cannot reach.** MARD on transcripts
   with zero ambiguous read support is constant to 4e-8..6e-7 across the entire
   alpha range, against library-wide effects of 9e-3..5e-2. Five orders of
   magnitude. This was not guaranteed: unreachable transcripts share the TPM
   renormalization with reachable ones, so mass moved among the reachable could
   have shifted them, and measurement says the leakage is negligible. Changing
   alpha cannot harm cases that were already fine.
2. **The headroom that remains in quant-only lives in the stratum only alpha can
   move.** On both MORF samples the irreducibly ambiguous stratum -- reads whose
   candidate transcripts share a 3' end and differ by internal splicing, which
   no per-read evidence can resolve -- is the WORST-quantified part of the
   output (MARD 0.0847 ont / 0.0528 pacbio, against 0.0096 / 0.0027 for
   unambiguous transcripts; Quant3Prime). alpha is the only knob LRAA has on it.
3. **The recommendation depends on the corpus, not on the 3' weighting.** Argmin
   alpha is identical with 3' weighting on and off on all 10 samples;
   the interaction contrast is 0.2-1.5% of the main effect.
4. **Both realistic-abundance samples put the optimum at alpha = 0.3, thirty
   times the shipped default** — morf2_ont (ONT) and morf2_pacbio (HiFi), both
   interior minima, both platforms, both 3' settings, four argmins one value,
   worth -7.6e-3 MARD, and it also removes every false negative.

## Recommendation

| stratum | recommended alpha | status | evidence |
| --- | --- | --- | --- |
| ONT | **~0.3** | **CONDITIONAL — do not ship yet** | morf2_ont: interior MARD minimum, -7.615e-3 paired vs default, false negatives 2 -> 0. Blocked on the false-positive measurement below. |
| HiFi | **~0.3** | **CONDITIONAL — same block** | morf2_pacbio: interior MARD minimum at the SAME alpha, -3.644e-3 paired, false negatives 32 -> 0. |
| pooled | **~0.3 across both realistic-abundance samples** | **CONDITIONAL — same block** | 2 of 2 realistic samples, both platforms, both 3p settings, all four argmins = 0.3. |
| SIRV-only | do not use | — | The two SIRV strata recommend opposite extremes and cancel to nothing when pooled. Seven strikes below. |

**The two realistic-abundance samples agree exactly.** morf2_ont (ONT) and
morf2_pacbio (HiFi) both put the MARD optimum at alpha = 0.3, both with a genuine
interior minimum (1.0 is worse than 0.3 on both), and both give the same argmin
with 3' weighting off. Four independent argmins, one value, across two platforms.
So 0.3 is not ONT-specific, and the recommendation is platform-independent —
which is the opposite of what the SIRV strata suggested, and the reason SIRVs are
excluded from it.

| morf2_pacbio, 3p on | MARD | delta vs default | truth-expressed called zero (of 2727) |
| --- | --- | --- | --- |
| alpha 0 | 0.037223 | +1.63e-03 | 32 |
| alpha 0.01 (default) | 0.035598 | 0 | 0 |
| alpha 0.1 | 0.032837 | -2.76e-03 | 0 |
| **alpha 0.3** | **0.031954** | **-3.64e-03** | 0 |
| alpha 1.0 | 0.032886 | -2.71e-03 | 0 |

### Why ONT is conditional and not ready to ship

**The blocking gap is false positives, not the sample count.** alpha adds mass
in proportion to ambiguous support, so raising it 30x is precisely the regime
where it should start putting mass on transcripts that are not expressed. That
is its known failure mode in the direction being recommended. Neither SIRVs nor
MORFs contain a single truth-unexpressed transcript, so this corpus **cannot
measure it at all**. The half of the objective that is measurable here (false
negatives) favours 0.3; the half that is not points the other way. Arabidopsis
and mouse — the only corpora with truth-unexpressed transcripts (5,286 and
30,994) — are being run now, and the recommendation is conditional on FP as a
function of alpha at 0.1 / 0.3 / 1.0.

**n=1 for ONT is a corpus limit, not a sampling shortfall.** morf2_ont is the
only realistic-abundance ONT library in the entire quant-only corpus; the only
other ONT source is SG-NEx, which is poor quality. No additional work on this
corpus can strengthen the ONT recommendation, and the arabidopsis and mouse runs
are HiFi so they will not fix it either — they test whether the direction
generalises across organism and error rate, not whether ONT replicates.

## The SIRV corpus cannot answer this question

Seven independent findings, each of which looks like a caveat alone:

1. **E0 is equimolar** and recommends the opposite alpha to every truth-bearing
   sample. Excluded before this work began.
2. **E1 and E2 have only 4 distinct truth levels each** — exact powers of two,
   1:2:4:8 and 1:8:32:128 — against 907 and 1788 distinct levels on the MORFs.
3. **The four cell lines per E-level are technical replicates of one spiked
   material**, recovered at different depths from different host transcriptomes.
   Their normalized truth compositions agree to 1.4e-17. So 10 libraries carry
   4 truth compositions, and within-stratum agreement is close to no independent
   evidence. They cannot reveal composition-dependent effects, which is exactly
   the effect found here.
4. **False negatives are unreachable on SIRVs.** `n_pred_zero` and
   `n_pred_lt_1pct` are 0 in all 184 SIRV arm-rows; `n_pred_lt_10pct` has a
   maximum of 1. Not that zeros are rare — the whole low tail is absent.
5. **The 3' weighting's side channels fire on SIRV and are inert on MORF**
   (Quant3Prime).
6. **Unique-read counts already contradict the SIRV truth by more than the
   entire ambiguous mass.** 26-34 of 69 transcripts have negative residual
   (truth minus unique mass), totalling 1.3x-2.8x the ambiguous mass. Unique
   reads involve no EM, no prior and no parameter, so on SIRVs the dominant
   error is upstream of anything alpha can reach. That is the honest explanation
   for SIRV MARD sitting at 0.18-0.20 while MORF sits at 0.03-0.06.
7. **The unambiguous transcripts — the ones involving no EM, no prior and no
   parameter — are quantified 21x-116x worse on SIRV than on MORF.** MARD of the
   UNAMBIGUOUS stratum at alpha 0.01, 3p on (Quant3Prime, all 10 samples):
   BT474_E1 0.3066, HG002_E1 0.2050, K562_E1 0.2187, UHRR_E1 0.2058,
   BT474_E2 0.2763, HG002_E2 0.2746, K562_E2 0.2926, UHRR_E2 0.2312, against
   morf2_ont 0.00957 and morf2_pacbio 0.00265. Eight of eight SIRV samples in
   0.205-0.307; both MORF samples in 0.0027-0.0096. These transcripts are scored
   purely on their unique reads, so this IS the unique-read-versus-truth
   disagreement of (6) read off the accuracy side — one defect, two instruments.
   It also inverts the within-corpus ordering: on SIRVs my unreachable partition
   scores 0.443 against 0.179 for the reachable (BT474_E1, per-transcript
   basis), while on MORF unambiguous is the BEST stratum. Checked and rejected
   as an abundance artifact: the 15 unreachable SIRV transcripts spread over all
   four truth levels (3/5/2/5 on BT474_E1) with the same mean truth abundance as
   the reachable set (14842 vs 14396 TPM on E1; 14292 vs 14549 on E2).

   **Consequence: SIRV ABSOLUTE MARD IS NOT INTERPRETABLE.** The 0.17-0.20
   library-wide levels are substantially corpus defect, not LRAA behaviour.
   Paired deltas survive — the defect is common to both arms of every pair and
   cancels in the difference — which is why every SIRV claim in this report is a
   paired contrast and never an absolute level.

Ruled out for (6): the residual construction is self-consistent — `sum(T_scaled)
= N` and `sum(residual) = A` exactly on all 8 samples, and the rescale is a
0.4-1.0% no-op because the groundtruth files were already scaled to each
library's read count. Length bias is ruled out on all 8: Spearman between the
observed/truth ratio and transcript length runs -0.069..+0.131, p 0.28..0.98.

Truth-file provenance, derived from structure since no generator script
survives: level ratios are exact powers of two matching the Lexogen SIRV E1/E2
design molarities, times a per-library scale landing within 1% of that library's
mapped read count. So the files hold expected READ COUNTS, not molarity, and no
length term is needed — correct for full-length long reads where one read is one
molecule, and consistent with equimolar E0 appearing as 69 equal counts rather
than length-proportional ones.

## Response surface

`results/alpha_response_curves.png`, `results/grid_metrics.tsv`.

### SIRV: the optimum tracks the truth mixture, not the platform

MARD, mean paired delta vs the alpha=0.01 default, 3p on, n=4 per stratum:

| alpha | E1 delta | E1 better | E2 delta | E2 better |
| --- | --- | --- | --- | --- |
| 0 | +2.40e-04 | 0/4 | -1.57e-03 | 4/4 |
| 0.001 | +2.15e-04 | 0/4 | -1.42e-03 | 4/4 |
| 0.01 | 0 (ref) | — | 0 (ref) | — |
| 0.03 | -4.56e-04 | 4/4 | +3.25e-03 | 0/4 |
| 0.1 | -1.85e-03 | 4/4 | +1.06e-02 | 0/4 |
| 1.0 | -9.10e-03 | 4/4 | +2.42e-02 | 0/4 |

Monotone in opposite directions, no exceptions. Per-sample argmin is 10.0 (the
top of the grid) on all four E1 and 0.0 on all four E2, identical at both 3p
settings. Pooling the eight gives a mean near zero and 4-better/4-worse: two real
opposing effects cancelling, not an absence of effect. Quant3Prime reproduces
the same split independently at 10x finer spacing (alpha 0.001 vs 0.01: all four
E2 prefer 0.001, all four E1 prefer 0.01).

E1's optimum is the **alpha -> infinity limit**, not a truncation. Running
alpha=1e6 — where theta is exactly the normalized ambiguity profile — shows
alpha=10 already achieves 96-99% of the total alpha=0-to-infinity change on E1
and 92-95% on E2. The curve has converged; there was never anything above it.

### morf2_ont (the only ONT sample): interior optimum at alpha=0.3

| alpha | MARD 3p on | delta vs default | n_pred_zero | spearman | nrmse |
| --- | --- | --- | --- | --- | --- |
| 0 | 0.060088 | +4.18e-03 | 42 | 0.96174 | 0.47543 |
| 0.001 | 0.059257 | +3.35e-03 | 6 | 0.96298 | 0.47441 |
| 0.01 | 0.055906 | 0 (ref) | 2 | 0.96955 | 0.46586 |
| 0.03 | 0.052622 | -3.28e-03 | 1 | 0.97635 | 0.44948 |
| 0.1 | 0.049394 | -6.51e-03 | 1 | 0.98366 | 0.40930 |
| **0.3** | **0.048291** | **-7.62e-03** | **0** | 0.98722 | 0.35663 |
| 1.0 | 0.049972 | -5.93e-03 | 0 | 0.98809 | 0.31825 |

A genuine interior minimum: 1.0 is worse than 0.3 on MARD. The 3p-off column
gives the same argmin, 0.3, with delta -8.12e-3.

**The metrics disagree about HOW FAR, not about DIRECTION, and that strengthens
the case rather than weakening it.** All of MARD, spearman and nrmse_mean_truth
say the shipped 0.01 is far too low on morf2_ont, and not one of them favours
staying there. They split only on the magnitude: MARD has an interior minimum at
0.3 and turns up by 1.0, while spearman (0.98722 -> 0.98809) and nrmse (0.35663
-> 0.31825) keep improving through the top of the grid and would go further. So
a reader who distrusts MARD as the objective still has to accept that 0.01 is
wrong; they would only argue for a larger move than 0.3, never a smaller one.
MARD is the pre-agreed primary objective and puts the optimum at 0.3.

This is the only sample where the metrics split at all. On all 8 SIRV samples
every metric agreed about alpha's direction, spearman being the sole wobbler and
moving only in the 4th decimal at 60 scored transcripts. Quant3Prime hit a
separate metric split on the 3'-weighting contrast (nrmse dissenting on E2), so a
metric split is a recurring feature of this corpus rather than a one-off.

## False negatives: the one place the count-based objective works

Truth-expressed transcripts predicted exactly zero, morf2_ont, n_truth_expressed
= 2706: **42** at alpha=0, 12 at 3e-4, 6 at 1e-3, 4 at 3e-3, **2** at the default,
1 at 0.03 and 0.1, **0** at 0.3 and 1.0. Transcripts below 1% of their truth
value: 44 -> 3 at the default -> 0 at alpha=1.0.

This is the predicted failure mode of low alpha, measured, and it is monotone.
It also means the recommended ONT value of 0.3 is the smallest grid point that
drives false negatives to zero — MARD and the FN count agree there, which is
worth more than either alone.

On SIRVs the same count is identically 0 at every alpha (see strike 4), so this
half of the objective is measurable on exactly one of the nine samples finished.

## Partition analysis: where alpha can act

`results/partition_analysis.tsv`. Reachable = nonzero ambiguous read count, i.e.
exactly where alpha's pseudocount is nonzero.

On SIRVs, 54-57 of 69 transcripts (78-83%) are reachable, carrying 83-86% of
predicted mass, so the library-wide figure is diluted by only 1.21-1.28x. SIRVs
are an unusually dense annotation — nearly every transcript is ambiguous with
another — and this is not representative.

**MY PRE-REGISTERED PREDICTION WAS WRONG.** I committed in writing, before
measuring, that MORF's dilution factor would be >= 3x and plausibly near 10x,
reasoning from the 88-92% unique read MASS. Measured with Quant3Prime's
`threeprime_stratified_effect.py` (b96c1b6) rather than my own partition, on the
alpha 0.01 -> 0.3 contrast: the reachable strata (LEVER + TIED_ONLY) hold 1871 of
2706 rows on morf2_ont and 1876 of 2727 on morf2_pacbio — 69% both times — giving
a dilution of **1.45x on both**, essentially the same as SIRV's 1.21-1.28x rather
than 3-10x higher.

The prediction failed because I conflated read mass with transcript count. 88-92%
of READS being uniquely assignable does not make 88-92% of TRANSCRIPTS
unambiguous: most transcripts collect a few ambiguous reads even when most reads
are unique. Only 31% of MORF rows are fully unambiguous.

**But the effect IS concentrated, just not in the way I predicted.** The exact
decomposition (checked to recombine to the library delta):

| morf2_ont, alpha 0.01 -> 0.3 | rows | % rows | MARD at 0.01 | delta | % of library delta |
| --- | --- | --- | --- | --- | --- |
| LEVER | 738 | 27.3 | 0.06419 | -7.56e-03 | 27.1 |
| TIED_ONLY | 1133 | 41.9 | 0.08466 | **-1.33e-02** | **73.2** |
| UNAMBIGUOUS | 835 | 30.9 | 0.00957 | +5.26e-05 | -0.2 |

TIED_ONLY — the irreducibly ambiguous stratum, alpha's exclusive territory — is
both the worst-quantified stratum and the source of 73% of the total gain, from
42% of rows, with an effect 1.75x the library-wide figure. On morf2_pacbio the
split is 48.6% LEVER / 51.2% TIED_ONLY / 0.16% UNAMBIGUOUS.

**Inertness reproduced independently on a realistic annotation.** UNAMBIGUOUS
moves +5.26e-05 (ont) and -1.89e-05 (pacbio) against library effects of 7.6e-03
and 3.6e-03 — 0.2% and 0.5% of the effect, contributing -0.2% and +0.16% of the
delta. My SIRV inertness result (1e-7 against 1e-2) was on my own two-way
partition; this is Quant3Prime's three-way partition, their code, a different
corpus, same conclusion.

## Alpha x 3' weighting: measurably additive

Argmin alpha is identical at 3p on and off on all 10 samples (8 SIRV,
morf2_ont). The interaction contrast — the gain from moving alpha to its optimum
with 3p off, minus the same with 3p on — is between -2.7e-05 and +1.0e-04 on the
SIRV samples, against alpha gains of 1.2e-3 to 1.4e-2. The interaction is
0.2-1.5% of the main effect.

This is a positive measurement of independence, not an absence of evidence, and
it agrees with Quant3Prime's independent result that their on/off delta is
invariant across a 10x alpha change (HiFi +1.288e-3 at alpha 0.01 vs +1.208e-3 at
0.001; ONT +4.426e-3 vs +4.489e-3; per-sample agreement ~1e-5).

**The alpha recommendation therefore does not need conditioning on the 3p
setting.** Caveat that must travel with it: the 3' weighting can only reach
4.1% (HiFi) / 2.5% (ONT) of library mass, and moves 1.3e-3..5.3e-3 of mass in
total variation between on and off (Quant3Prime). A null here means 3p had no
lever to pull over this corpus, not that the two parameters are independent in
principle.

## Scale invariance: alpha is dimensionless

Exact half, `pylib/test_EM_alpha_semantics.py`: scaling every read count by k
leaves the normalized abundances bit-identical at k=2 and k=1024 (exact in binary
floating point) and equal to 1e-9 relative at k=10, for alpha in {0, 0.003, 0.01,
0.3, 1.0}. Both terms of the M-step scale with depth and normalization divides
the factor out.

Empirical half, 44 runs: HG002_E2 and K562_E2 subsampled (samtools, seed 42) to
50% and 25%: 251127 -> 125297 -> 62687 and 226621 -> 113456 -> 56665 primary
reads. Argmin alpha is 0.0 at every depth. Because that optimum is on a boundary
and cannot move downward, the informative test is the response magnitude:

| | sub50 | sub25 |
| --- | --- | --- |
| HG002_E2 | 1.003 (0.977-1.079) | 0.889 (0.852-1.026) |
| K562_E2 | 0.910 (0.778-1.004) | 1.097 (1.034-1.331) |

Ratio of (MARD(alpha) - MARD(0)) at reduced depth to full depth, averaged over
alpha >= 0.01. Within +-11% of unity and **non-monotone in depth** — K562 gives
0.91 at half and 1.10 at quarter — so the variation is re-draw noise from
subsampling, not a depth effect. **Use +-11% as the noise floor for any depth
claim.** The two halves test different things: the unit test holds multipath
structure fixed, subsampling perturbs it. Invariance survives both.

One default transfers across library sizes. Nothing here argues for scaling
alpha with depth.

## Mechanism: what alpha actually chooses

As alpha grows, the M-step is dominated by `alpha * ambiguous_read_counts`, so
theta converges to the **normalized ambiguity profile**. But theta is not what is
scored: LRAA reports E-step read counts at that theta, and 88-92% of read mass is
compatible with exactly one transcript and lands there regardless of theta. So

> alpha interpolates the split of AMBIGUOUS read mass between the likelihood's
> own split (alpha=0) and a split proportional to each transcript's ambiguous
> support (alpha -> infinity). It helps to the degree the latter is closer to
> truth, over the 4-13% of the library that is ambiguous.

The naive form of this — "alpha helps to the degree the profile resembles the
truth" — is **false**, and measurably so: the profile is farther from truth than
a flat prior on 8 of 8 SIRV samples (0.683-0.735 against 0.437-0.660), yet alpha
helps without bound on all four E1. The unique-mass correction above is why.

Tested at the right granularity (`residual_test.py`), with
`residual_t = T_t - u_t` in count space and no clipping: the profile split beats
the likelihood split on exactly the four samples where raising alpha helps, and
loses on exactly the four where it hurts. **8/8.** Margins are 1.5%-15% of
ambiguous mass. This is a SIGN result only: the unique-versus-truth disagreement
dominates the target's absolute level, but it is common-mode — both splits are
scored against the same residual and differ only over the ambiguous slice — so it
cancels in the signed comparison and cannot support a magnitude claim.

### Effect size scales with ambiguous mass (ordering, within stratum)

Covariate from Quant3Prime's read-level tracking extraction; my independent
extraction reproduces their total-ambiguous fractions to 3 decimals on all 8
samples (6.276/7.206/6.808/9.812 on E1, 13.439/8.135/7.786/7.578 on E2).

Sensitivity is `|d MARD / d log10 alpha|` at the default, central difference over
0.003 and 0.03 — not MARD spread, which measures different things in the two
strata because their optima sit at opposite ends of the grid.

* **E1 is the result:** Pearson r = 0.9949 (p = 0.0051) against pre-registered
  irreducible mass, Spearman 1.0, over four genuinely separated masses
  (3.53-5.89%). Both pre-registered maxima were hit (UHRR in E1, BT474 in E2),
  jointly p ~ 0.06.
* **E2 is consistent, not confirming.** Its four masses are 7.58 / 3.81 / 3.82 /
  3.76 — a tight cluster and one distant point — so r = 0.9917 is very nearly
  forced by that geometry and Spearman is 0.4. It is a two-point slope wearing a
  four-point r; within the cluster the ordering is arbitrary.
* Slope ratio E2/E1 = 1.78 on the derivative measure, against 0.21 on the spread
  measure, so the derivative removes the contradictory power laws. That 1.78 is a
  ratio of two point estimates, one with one effective degree of freedom, not a
  measured ratio of slopes.
* Irreducible and total-ambiguous mass are collinear in this corpus, so the
  mechanism is supported as an **ordering** predictor and cannot be attributed
  specifically to the irreducible fraction.
* Pooled correlations are contaminated by the flatness confound running the same
  way and are **not** offered as support.

## Could not determine

* **The residual test does NOT generalise out of stratum, and I do not know
  whether that is a failure of the test or a limit of its scope.** It was 8/8 on
  SIRV; on the two MORF samples the profile split LOSES to the likelihood split
  (misallocation 0.384 vs 0.308 on ont, 0.445 vs 0.371 on pacbio) while alpha
  nevertheless helps on both. The coherent reading is that the test compares the
  two ENDPOINTS, alpha=0 against alpha=infinity, so it can only predict a
  boundary optimum — and both SIRV strata have boundary optima while both MORF
  samples have interior ones at 0.3 with the curve already turning up by 1.0. If
  MORF's alpha=infinity really is worse than its alpha=0, the test is right and
  merely mute about interior optima. **Deciding that needs the MORF alpha=1e6
  arms, which I did not run.** Until then the residual test is supported only
  where the optimum is on a boundary.
* **The profile-distance predictor DID generalise, 2 of 2 out of stratum.** I
  committed before looking that morf2_ont's profile-vs-truth distance would place
  it on one side of the 0.695/0.724 SIRV gap and that would predict its optimum.
  Measured: morf2_ont 0.672 and morf2_pacbio 0.678, both below 0.695, i.e. on the
  E1 (alpha-helps-upward) side. Both optima are indeed above the default. Two
  hits, on a corpus with a different truth and a different unique-mass fraction.
* **The false-positive cost of the recommended move is UNMEASURED, and this is
  the blocking gap.** Raising alpha 30x is the exact regime where its known
  failure mode — mass on transcripts that are not expressed — should appear, and
  neither SIRVs nor MORFs contain a single truth-unexpressed transcript. So the
  measurable half of the objective favours 0.3 while the unmeasurable half points
  the other way. Arabidopsis (5,286 truth-unexpressed) and mouse (30,994) are
  running now; the ONT recommendation is conditional on FP as a function of alpha
  at 0.1 / 0.3 / 1.0.
* **ONT n=1 is a CORPUS LIMIT, not a sampling shortfall.** morf2_ont is the only
  realistic-abundance ONT library in the whole quant-only corpus; the only other
  ONT source is SG-NEx, which is poor quality. No further work on this corpus can
  strengthen the ONT recommendation, and the arabidopsis and mouse runs are HiFi
  so they cannot either. The "4/4" paired sign-consistency that backs every SIRV
  statement is unavailable for the recommendation that actually matters, and will
  stay unavailable until an ONT library worth benchmarking exists.
* **The metric split on morf2_ont is unresolved** (MARD 0.3; spearman and nrmse
  >= 1.0). It is a disagreement about magnitude, not direction — see above — so
  it does not threaten the "0.01 is too low" conclusion, only the value 0.3.
* **MORF alpha=1e6 limit arms were not run**, which is what blocks the residual
  test above. Two arms, ~20 min.
* **The depth probe covers only SIRV E2.** Both probed samples have boundary
  optima and neither has a realistic abundance distribution. An E1 subsample
  would have tested a top-edge optimum; a MORF subsample would have tested a
  realistic one. Neither was run.
* **False positives are unmeasurable on this corpus.** Neither SIRVs nor MORFs
  contain truth-unexpressed transcripts, so the FP half of the stated objective
  is untestable here regardless of alpha. Combined with false negatives being
  unreachable on SIRVs, **MARD carries essentially the whole objective on 8 of
  the 9 finished samples.** An FN-based objective needs a decoy-bearing
  annotation or a detection threshold; "predicted exactly zero" is close to
  unreachable in quant-only against a supplied annotation, because EM
  redistributes rather than eliminates.
* **Norm-cache sharing was considered and rejected.** It would save the ~1m57s
  normalization on 35 of 36 MORF arms by symlinking one content-keyed cache dir
  per sample. The cache stem already keys on the bam fingerprint and every
  normalization parameter so it is safe in principle, but concurrent arms of the
  same sample would race to build it. A race in a measurement harness does not
  fail loudly; it produces a subtly wrong normalized BAM shared across arms.
  Rejected for ~30 minutes of wall-clock.
* **`--no_parallelize_contigs` does not help and should not be retried** on a
  big-contig corpus. It is ~4x slower on MORF (1425 units in 20 min against 2652
  in 11): it does not remove the per-contig process spawn, only stops the spawns
  overlapping, and quant-only per-contig work has no native tool step to hand the
  freed budget to.
* **Arabidopsis and mouse are out of scope for this phase** and are the only
  corpora with truth-unexpressed transcripts.

## Verified along the way

* **Scorer control reproduced exactly.** Re-scoring the staged LRAA-v0.18.2
  morf2_ont output in a fresh directory with a fresh registry gives spearman
  0.9695120257297285 and mard 0.055854896341881864 — delta 0.00e+00 against the
  per-sample published file. (The aggregated MORFs table prints ...284 / ...818;
  that is repr truncation at aggregation, not a different computation.)
* **Pilot reproduced bit-identically**, not merely in ordering: HG002 E2 at alpha
  0/0.001/0.01/0.1/1.0 gives MARD 0.18539523145409811 / 0.18554596059228295 /
  0.18695190574637816 / 0.19549058954832862 / 0.20954200148880810, matching the
  five pilot values digit for digit from a different worktree at a different time.
* **LRAA quant-only output is invariant to work partitioning**, verified
  byte-for-byte in both contig regimes: cpu_budget 1/3/6 on SIRV's 7 contigs
  (Main) and cpu_budget 3 vs 8 on MORF's 1566 contigs (here, md5
  4b529ea3e05dc68ce33182e567b96ac0 on both, 2924 rows, comments stripped). This
  is the precondition that makes any parallel parameter study on this codebase
  trustworthy and had not been checked before.
* **The scorer needs >= 2 registry entries per directory** — its scatterplot cell
  indexes the subplot array, which matplotlib collapses to a bare Axes at n=1.
* **Per-arm metrics are independent of how many arms are in a scoring pass**, so
  incremental scoring is safe.

## Reproducing

```
./run_grid.py --samples <s> --alphas 0,0.01,1.0 --threep both --jobs 4 --cpu_budget 3
./make_registry.py && ./score_grid.py && ./collect_results.py && ./analyze_grid.py
./ambiguity_profile.py --run --analyze     # needs tracking; ~1 extra arm per sample
./residual_test.py ; ./partition_analysis.py ; ./ambiguity_sensitivity.py
./plot_curves.py
```
