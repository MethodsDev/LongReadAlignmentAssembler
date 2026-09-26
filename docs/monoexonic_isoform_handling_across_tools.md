# Monoexonic isoform handling: IsoQuant, bambu, StringTie 3 (and LRAA for contrast)

**Purpose.** Hand-off brief for an agent working on monoexonic (single-exon / unspliced) transcript
discovery and filtering. Answers: how each tool represents monoexonic models, whether novel
monoexonic discovery is on by default, and which thresholds apply to monoexonic models that do not
apply to multiexonic ones.

**Provenance.** Compiled 2026-09-19 from direct upstream source inspection plus publications.
Versions inspected:

| tool | version | commit / tag |
|---|---|---|
| IsoQuant | 4.0.0 | `7d8268918a770b8d0c6925e8cea99d6c969b4eae` (tag `v4.0.0`) |
| bambu | 3.13.1 devel | `066b92ff222df50b9c74c9eb8ae306d2947ffc8c`; decisive defaults re-checked in tag `v3.12.1` (`8991474`) |
| StringTie | 3.0.3 | `3436ad6dfd0ffc806a94086cf747ac6ff2b0dc19` (tag `v3.0.3`, 2025-11-13) |
| LRAA | working tree, branch `devel` | this repository, as of the date above |

Where code and publication disagree, **code is authoritative** and the disagreement is called out.
Claims that could not be established from source are listed in the final section; do not fill those
gaps from general knowledge of the tools.

---

## 0. Read this first: two distinctions that are easy to collapse

1. **Discovery eligibility is not quantification.** bambu's `min.txScore.singleExon = 1` and
   IsoQuant's platform-dependent `--report_novel_unspliced` are *discovery* policies: they decide
   whether a novel monoexonic model is added to the annotation. Both tools quantify monoexonic
   transcripts that are already in the supplied annotation regardless, and in bambu the read
   classes rejected for discovery still participate in quantification. A shared, exon-count-agnostic
   quantification path therefore does **not** imply a shared output gate, and "novel single-exon
   discovery is off" does **not** mean "no single-exon rows in the output".
2. **Annotated vs novel.** Every tool here treats a monoexonic transcript that the supplied
   annotation already asserts differently from one inferred de novo. Statements below are about
   *novel* monoexonic models unless explicitly marked "known/annotated".

---

## 1. Summary table

| | novel monoexon default | separate construction path? | decisive monoexon-specific gate | genomic internal-priming check? |
|---|---|---|---|---|
| **IsoQuant 4.0.0** | **ON** for PacBio/assembly, **OFF** for ONT | Yes — bypasses the intron graph | external polyA/polyT tail is **mandatory, unconditional** | No (read sequence only) |
| **bambu 3.13.x** | **OFF** for all data | Yes — separate unspliced RC path + separate XGBoost model | `min.txScore.singleExon = 1` applied as strict `>` to a probability ⇒ unsatisfiable | No (10 bp genomic A/T counts as ML features only) |
| **StringTie 3.0.3** | **ON**; no mono-only switch exists | No — same splice graph, degenerate `gno == 3` case | aligned terminal A/T discard for `exons.Count() <= 2` | No (read sequence only) |
| **LRAA** | ON; under `--HiFi` requires a terminal TSS/PolyA node | Yes — separate SE splice graph | `min_monoexonic_TPM`, read-span coherence, adjusted-TPM ratio | **Yes** — `looks_internally_primed()` against contig sequence |

---

## 2. IsoQuant 4.0.0 — polyA evidence is the currency

### Representation
A monoexonic model is an ordinary `TranscriptModel` with a one-element `exon_blocks` and
`intron_path = ()`; there is no distinct class. The *construction route* is separate:
`GraphBasedConstructor.construct_fl_isoforms` skips empty intron paths outright
(`if not intron_path: continue`, `isoquant_lib/model_construction/fl_graph_constructor.py`), so novel
monoexons are never vertices or zero-intron paths in the intron graph. They are built by
`AssignmentBasedConstructor.construct_monoexon_novel`
(`isoquant_lib/model_construction/assignment_based_constructor.py`), which clusters candidate reads
**solely by 3' polyA/polyT position** within `apa_delta = 50 bp`. The opposite (5') boundary is
simply the most extreme observed read start/end in the cluster. Models are typed
`novel_not_in_catalog` with a new gene ID, then rejoin the shared `ModelFilter` / end-refinement /
counting pipeline.

### Default on/off
Controlled by `--report_novel_unspliced` (`-u`). Parser default is `None`; the decisive line is
`if args.report_novel_unspliced is None: args.report_novel_unspliced = strategy.novel_monoexonic`
(`isoquant.py`, `set_model_construction_options`).

| data type | preset | novel monoexon |
|---|---|---|
| `pacbio_ccs` / `pacbio` | `default_pacbio` | **ON** |
| `pacbio` + `--fl_data` | `fl_pacbio` | **ON** |
| `nanopore` / `ont` | `default_ont` | **OFF** |
| `assembly` / `transcripts` | `assembly` | **ON** |

The gate wraps only `construct_monoexon_novel`; **known** monoexon recovery is not disabled by it.
`--model_construction_strategy sensitive_ont` or `all` also yields ON; `reliable` yields OFF.

### Mandatory boundary evidence
Candidate reads must be single-exon, `polyA_found`, non-multimapping, and unassigned-or-inconsistent.
Clustering then accepts **only** `external_polya_pos` or `external_polyt_pos` — an internal-only tail
call never enters a cluster. This requirement is hard-coded and survives
`--polya_requirement never`, which modifies only `requires_polya_for_construction`,
`require_monointronic_polya`, and `require_monoexonic_polya` (the latter is the *known*-monoexon
switch). Docs concur: "polyA tails are always required for reporting novel unspliced isoforms".
`--polya_trimmed stranded|all` synthesizes external positions so trimmed data can satisfy the path.

No TSS is required. A trained TSS model refines 5' ends only under `--fl_data`, and
`TranscriptEndProcessor.correct_novel_transcript_ends` "never creates or drops a model".

### Monoexon-specific thresholds
Defaults shown PacBio / ONT / assembly.

| parameter | default | gates |
|---|---|---|
| `--report_novel_unspliced` | true / false / true | whether the novel-mono constructor runs |
| `min_novel_count` (preset-internal, no flag) | 2 / 3 / 1 | reads per 3' cluster; also the absolute floor in `ModelFilter` |
| `min_mono_count_rel` (preset-internal) | 0.005 / 0.02 / 0.01 | final floor = `max(min_novel_count, rel * overlapping_component_max_coverage)`; applies to monoexons **and one-intron models** |
| `simple_models_mapq_cutoff` (internal) | 30 | mean supporting-read MAPQ for ≤2-exon models |
| `--simple_alignments_mapq_cutoff` | 1 | per-read MAPQ/multimapper prefilter for ≤2-exon alignments |
| `apa_delta` (from `minor_exon_extension`) | 50 bp | 3'-position clustering window |
| `min_mono_exon_coverage` (internal) | 0.75 | **known monoexons only**: breadth (binary union) of annotated exon covered |
| `require_monoexonic_polya` | PB true / ONT true / assembly false; `auto` flips to true when sample polyA fraction ≥ 0.7 | **known monoexons only** |
| candidate-overlap predicate (not configurable) | reject if any exon of an already-built model overlaps by `> candidate_length/2` | **strand-blind** |

### Multiexon-only thresholds monoexons escape
All intron-graph parameters: `min_novel_intron_count`, `graph_clustering_distance`,
`graph_clustering_ratio`, `min_novel_isolated_intron_abs`, `singleton_adjacent_cov`,
`terminal_position_abs` / `_rel`, `terminal_internal_position_rel`, `min_novel_count_rel` (≥2-intron
models only), `require_monointronic_polya`, `report_canonical_strategy`, and the `--use_replicas`
multi-file requirement. Conversely, ≥3-exon models escape both MAPQ gates.

### Context filters
- Purely intronic reads are `noninformative` / `genic_intron`, which counts as unassigned, so they
  **are** eligible candidates.
- `is_internal_monoexonic_read` suppresses a candidate whose 3' end matches a spliced terminal exon
  within `apa_delta` and which lies internally within it.
- Surviving same-strand candidates are later merged into the host gene by `TranscriptToGeneJoiner`
  (shared junction, or ≥ `STRONG_OVERLAP_FRACTION = 0.5` overlap of the shorter gene interval).
- Strand comes from tail orientation: polyA only ⇒ `+`, polyT only ⇒ `-`, else `.`.

### Quantification
**Not EM in v4.0.0.** `ReadWeightCounter` gives unique reads weight 1 and splits ambiguous reads
equally (`1.0 / feature_count`) only when the strategy admits them; transcript-level default is
`unique_only`, so ambiguous reads contribute zero. TPM is `count * 1e6 / total`, with no effective
length correction. Two monoexon-aware details: `fsm_only` accepts `mono_exon_match` as a full match,
and `confirms_feature` lets an unspliced read confirm a zero-intron target but not a spliced one.
Ambiguity resolution is *more* aggressive for monoexons: under the default presets
`resolve_ambiguous = monoexon_and_fsm`, so every ambiguous unspliced match is score-resolved
(Jaccard − flanking fraction) while spliced matches are resolved only when an FSM is present.

---

## 3. bambu 3.13.x — off by an unsatisfiable threshold

### The off-switch is a number, not a boolean
`setIsoreParameters()` (`R/bambu_utilityFunctions.R`):

```r
min.txScore.multiExon  = 0,
min.txScore.singleExon = 1,
```

The monoexon predicate is strict: `NSampleTxScore = sum(txScore > min.txScore.singleExon)` in
`makeUnsplicedTibble()` (`R/bambu-extendAnnotations-utilityCombine.R`). `txScore` is the output of an
XGBoost `binary:logistic` model, i.e. a probability in [0,1], so `> 1` is unsatisfiable: **no novel
unspliced candidate can pass by default.** The vignette states it and gives the enabling call:
`opt.discovery = list(min.txScore.singleExon = 0)`. A non-zero user value acts as a genuine strict
TPS gate, not merely an on/off.

**NDR does not override this.** The TPS / read-count / gene-fraction / sample gates run in
`filterTranscripts()` *before* `calculateNDROnTranscripts()`, so even `NDR = 1` reports no novel
monoexons while the default stands.

### Separate machinery throughout
- `isore.constructReadClasses()` splits on `elementNROWS(readGrgList) == 1` and routes one-exon
  alignments to `constructUnsplicedReadClasses()`. Spliced RCs are keyed by corrected junction chain;
  unspliced RCs are formed by containment / interval reduction.
- Reads contained in a reference exon become `unsplicedWithin` with coordinates
  `start = max(start), end = min(end)` over hits, and are removed from novel consideration.
  Remaining reads are `reduce()`d into union components ⇒ `unsplicedNew`.
- Boundary estimation differs: spliced RCs use the 20th-percentile start / 80th-percentile end;
  unspliced use interval unions or reference intersections.
- `getTranscriptScore()` predicts with `transcriptModelME`, then **overwrites** `numExons == 1` rows
  with `transcriptModelSE`. Nine features: scaled/log read count, gene read proportion, start SD,
  end SD, genomic A and T counts at each end, strand bias.
- Cross-sample aggregation differs: spliced uses `pmax`; unspliced uses
  `weighted.mean(txScore, readCount_tmp)` — despite the result being named `maxTxScore`.

### Extra gates once enabled
| gate | default | note |
|---|---|---|
| hard-coded `readCount > 1` in `extractNewUnsplicedRanges()` | effectively 2 reads | `min.readCount = 1` does **not** bypass it |
| `min.exonOverlap` veto in `addNewUnsplicedReadClasses()` | 10 bp | any overlap ≥ this with an annotated transcript exon **discards** the candidate; only `is.na(overlap)` survives |
| `stranded` | FALSE | default unstranded reduction/containment ignores strand and stores `*` |

Common gates still apply: `min.readCount = 2`, `min.readFractionByGene = 0.05`,
`min.sampleNumber = 1`, then the pooled NDR stage.

### Multiexon-only thresholds monoexons escape
`min.txScore.multiExon` (0), the `highConfidenceJunctionReads` requirement, `min.exonDistance` (35)
and `min.primarySecondaryDist*` (5) in the discovery path, and `remove.subsetTx = TRUE` — whose
string predicate targets `compatible` spliced classifications and therefore never matches
`unsplicedNew`. The monoexon annotation-overlap veto is the stricter substitute.

### Boundary / internal priming
No TSS, TES, poly(A) tail, or PAS is required. Consistency is preferred *probabilistically*:
`startSD` / `endSD` are features, and `countPolyATerminals()` counts A and T in **10 genomic bases**
at each terminus (`numAstart`, `numAend`, `numTstart`, `numTend`) as four more features. These are
genomic-context features, not observed read tails, and there is no hard internal-priming rule.

### Quantification
Shared path, **no `numExons` branch**. Discovery-rejected read classes still participate. But
compatibility for a one-exon RC has no intron constraint in `myGaps()`, so its compatible set is
intrinsically larger; since `uniqueAval = aval * (!multi_align)`, monoexons mechanically accrue fewer
unique reads with no compensating rule. `degradationBias = TRUE` applies uniformly; there is no
monoexon-specific effective-length correction.

---

## 4. StringTie 3.0.3 — same graph, and the manual is wrong in two places

### Representation
No separate monoexonic subsystem. `CTransfrag` stores every model as a node vector plus sparse
`GBitVec` pattern; a pure unspliced locus is the degenerate source–node–sink graph, identified in
`get_trf_long()` as `gno == 3` and literally named `singleExonGene`. `-L` and `--mix` use different
*mode-level* routines (`get_trf_long()` / `get_trf_long_mix()`), not mono-specific ones.

The one clear mono-specific support gate is node-coverage eligibility in `get_trf_long`
(`rlink.cpp`):

```c
gno==3 && (polyStartUnaligned>10 || polyEndUnaligned>10 || abundance >= CHI_WIN/2)   // CHI_WIN = 100
```

i.e. unaligned-tail evidence at either end **or** abundance ≥ 50, as alternatives.

### Default on/off
**ON**, and there is no `--report-novel-unspliced` equivalent in the v3.0.3 option set. `-e -G ref.gtf`
disables *all* de novo discovery, not monoexons specifically.

### Two documentation/code contradictions — do not trust the manual here
1. **`-L` does not set `-s 1.5`.** Help and README say `-L` "enforces `-s 1.5 -g 0`". The source sets
   `bundledist = 0`; the `singlethr = 1.5` assignment is **commented out**, so the compiled default
   **4.75** stands unless `-s` is passed explicitly (`stringtie.cpp`, `processOptions`).
2. **`-s` is not an exon-count-aware threshold.** Both active consumers of `singlethr` test
   `guided && !...->guide && abundance < singlethr` inside loops over already-kept long clusters and
   **never inspect exon count** (`rlink.cpp`, `process_transfrags`). Symmetrically, `-c`
   (documented multi-exon-only) is applied to rescue transfrags without an exon-count check.

### Thresholds
| parameter | default | note |
|---|---|---|
| `-s` / `singlethr` | **4.75** compiled | see contradiction 2; not a verified mono-only output filter |
| `-m` | 200 nt | min assembled length, all novel models; guides bypass |
| `-c` / `readthr` | 1 read/bp | used as `guide \|\| abundance>=readthr` without exon-count test |
| `-f` / `isofrac` | 0.01; floored to **0.10** in `--mix` (`ERROR_PERC`) | shared isoform-fraction gate |
| terminal aligned A/T screen | ≤20 bp window; ≥5 A/T, ≥0.80 fraction, or 5-base anchored run | discards novel `exons.Count() <= 2` alignments |
| whole terminal-exon A/T | ≥0.80 | removes an A-rich last / T-rich first exon on ≥2-exon alignments |
| `POLY_TAIL_STOP_COUNT` | 8 | aggregated unaligned tails promote `hardstart`/`hardend` |

Multiexon-only gates monoexons escape: `-a` (10 bp junction anchor), `-j` (junction coverage),
`-E` (25 bp splice-site correction window), the 25-bp anchor for introns > 100 kb, and the
consecutive-splice-edge long-read witness check.

### Boundary / internal priming
No TSS or PAS requirement: an unabsorbed candidate is retained if it is a guide or has both
`(longstart || hardstart)` and `(longend || hardend)` — ordinary alignment endpoints suffice.
Optional `--ptf` point features supply `GPFT_TSS` / `GPFT_CPAS` hard endpoints; guide ends are always
hard. Internal-priming detection is read-level and aggressive: any aligned terminal T (left) or A
(right) discards a novel ≤2-exon alignment outright, with the proposed unaligned-tail exemption still
commented out. Guide overlap (≥5 bp, `BundleData::evalReadAln`) exempts a read from both the ≤2-exon
discard and terminal-exon removal.

### Quantification
No post-assembly EM: reconstruction and quantification are simultaneous via generalized maximum flow
(`push_max_flow()` / `long_max_flow()`). Coverage is base-weighted and divided by assembled exonic
length — ordinary normalization, not an effective-length model. Multimappers are `NH`-divided
identically for all exon counts. Mono and multi share compatibility and flow allocation; the only
difference is the `gno == 3` eligibility predicate and the absence of splice-edge constraints.

---

## 5. LRAA, for contrast

Grounded in this repository (`pylib/`, `LRAA`):

- **Separate SE graph.** `ME_graph_emits_monoexonic_models = False` by default; single-exon
  reconstruction is routed to a separately built `Splice_graph(restrict_splice_type="SE")` via
  `build_SE_transcripts`. `SE_subtract_intronic_background = True` subtracts the local intronic
  coverage floor before SE segmentation.
- **Monoexon-specific filters** (all in `TranscriptFiltering`, defaults in `LRAA_Globals.config`):
  `min_monoexonic_TPM` 1.0, `min_monoexonic_read_span_peak_frac` 0.5 (do supporting reads stack on a
  common interval, or merely tile?), `min_monoexonic_adjusted_TPM_ratio` 0.20,
  `min_monoexonic_supporting_cells` 5 (single-cell only).
- **3'-end evidence requirement** (replaces the former HiFi-only boundary rule). In
  `filter_monoexonic_isoforms_by_terminal_evidence`, gated by
  `require_terminal_evidence_for_monoexonic` (default True) and run after PAS annotation, a
  single-exon model is retained if **any** of three channels vouches for its 3' end: an inferred
  PolyA site, a canonical PAS hexamer upstream of its terminus, or >=1 assigned read whose own soft
  clip there is >=80% A/T over the 20 bases nearest the alignment
  (`max_dist_tail_evidence_to_terminus` 100). Reference-containing models bypass via
  `reference_model_reprieved()`; multi-exonic models are untouched.
  The rule it replaces (`require_terminal_feature_for_monoexonic`, demanding
  `has_TSS() or has_PolyA()`) measured whether terminal-feature *inference* succeeded rather than
  whether a model is real: HiFi models satisfy it at ~88% while ONT models satisfy it at ~5%
  *including multi-exonic ones*, so it deleted 99.98% of monoexonic models on A549 ONT and 55% of
  all output on MCF7. The disjunction is what makes one rule portable, because the channels are
  near-complementary across platforms -- measured on MCF7 ONT the retained models are carried by
  PAS 65.2% / read tail 34.4% / inferred PolyA 0.3%, and on BT474 PacBio HiFi by inferred PolyA
  89.7% / PAS 10.3% / read tail 0.0%, the last because Kinnex strips tails before alignment.
  End-to-end on MCF7 ONT it retains 88.6% of models matching an annotated single-exon transcript
  (against 3.0% for the old rule) while removing 80.4% of the antisense class; on BT474 HiFi it
  retains 89.7% of monoexons against the old rule's 86.1%, for 0.16% of total output.
- **Proximal-window tail scoring.** Tail *evidence* is scored over the 20 bases of the soft clip
  nearest the alignment (`polyA_tail_proximal_window`, `min_proximal_tail_base_frac`), separately
  from the whole-clip test that decides whether to *strip* a clip and move a boundary. An ONT cDNA
  clip is tail + adapter + barcode, so whole-clip averaging fires on 1.1% of reads at GENCODE 3'
  ends where the proximal window fires on 65.2% (400 sites, 44,605 reads, SGNex MCF7). The windows
  are asymmetric: a reverse alignment stores the reverse complement, so its 3' end is the LEFT clip
  and the alignment-proximal end is that clip's tail. Scoring both ends with the same `[:20]` would
  report 0.9% instead of 54.6% on reverse reads.
- **Annotation-overlap veto: deliberately absent.** bambu's `min.exonOverlap` predicate was
  measured on 4,786 restored MCF7 monoexons and is anti-correlated with quality here: it vetoes
  129/129 models matching an annotated single-exon transcript while keeping 70% of the antisense
  class. It answers "is this redundant with the annotation" for novel-discovery gating, not "is
  this real", and it cannot fire at all in de novo mode.
- **Genomic internal priming.** `Util_funcs.looks_internally_primed(contig_seq_str, position, strand)`
  is applied at PolyA-site identification in `Splice_graph._incorporate_PolyA_objects`, with
  `restrict_internal_priming_filter_to_monoexonic = True`, a reference-3'-end reprieve, and a
  `reject_internally_primed_polyA_sites` policy defaulting to `"spliced_only"`.
- **Quantification is EM** (`Quantify.py`), unlike all three tools above.

### Where LRAA is distinctive
1. **Only tool of the four with a reference-genome internal-priming check.** IsoQuant, bambu and
   StringTie all detect A/T from the *read* sequence (or, in bambu, use a 10 bp genomic window purely
   as an ML feature). None performs a downstream-genome A-richness veto or a canonical PAS-motif test
   in the inspected paths.
2. **Only tool with read-span coherence tests.** `min_monoexonic_read_span_peak_frac` and
   `min_monoexonic_adjusted_TPM_ratio` ask whether supporting reads describe one molecule or a
   covered region. No counterpart in any of the three.
3. **Now the only tool requiring 3'-end evidence as a disjunction rather than a fixed channel.**
   IsoQuant demands a per-read polyA tail for novel monoexons; LRAA accepts a tail, an inferred
   PolyA site, or a genomic PAS motif, whichever the input can supply. This matters because the
   channels are not interchangeable across platforms: IsoQuant's tail requirement is enabled
   specifically for PacBio, yet Kinnex/Iso-Seq strips tails before alignment (98.4% of BT474 reads
   carry a 3' clip <7 bp), so on that input the rule is unsatisfiable rather than merely strict.
4. **Conversely**, LRAA still has no 5'-end counterpart to any of this. All three of its channels
   interrogate the 3' terminus, so a 3'-anchored degradation fragment passes.

---

## 6. Version drift: published values that no longer match code

- **IsoQuant.** The 2023 paper's novel-model support ("at least five FSM reads, three for PacBio") no
  longer matches v4.0.0 (`min_novel_count` = 3 ONT / 2 PacBio). The paper's known-monoexon rule
  (unique read + confirmed polyA) omits the ≥75% exon-breadth requirement now present, and
  assembly/PB-FL presets do not require polyA unless the `auto` rule (≥70% polyA reads) engages.
- **bambu.** The paper defines multi-sample TPS as the **maximum** across samples; code does that for
  spliced candidates but a read-count-weighted **mean** for unspliced ones. The paper describes
  terminal A/T over 20 bp; code uses `width = 10` at each end. Paper describes `min.exonOverlap` as a
  merge tolerance; code uses it as a discard veto. Vignette says de novo default NDR < 0.1; code sets
  `NDR = 0.5` when no annotation is supplied.
- **StringTie.** The StringTie3 paper says the aligned-terminal discard rule applies to single-exon
  alignments; code applies it at `exons.Count() <= 2`. The paper's partial-terminal-segment
  "shorten to 3 bp" rule has no active implementation (`g_longread_shortened` is declared and printed
  but never incremented). The paper's 5-bp/20-read CPAS anchoring has definitions
  (`CPAS_POS_BIN = 5`, `CPAS_MIN_SUPPORT = 20`, `cluster_positions_with_counts()`,
  `add_cpas_trimpoint()`) but no located call sites; active hardening occurs at 8 aggregated
  unaligned-tail reads instead. `-M` is parsed but has no located consumer.

---

## 7. NOT ESTABLISHED — do not assert these without further work

- **bambu:** whether terminal genomic A-richness raises or lowers `txScore`, and its effect size. The
  pretrained single-exon model is serialized in `R/sysdata.rda` / `inst/extdata/defaultModels.rds`.
  Also: `lmNDR.SE` is created by `trainBambu()` but has no located consumer in `recommendNDR()`, so a
  distinct single-exon NDR calibration is **not** implemented in the inspected filter.
- **IsoQuant:** whether any minimum novel-monoexon *length* exists (none in the constructor or
  `ModelFilter`, but an indirect upstream limit was not ruled out); whether
  `min_novel_isolated_intron_rel` has any consumer (appears vestigial in v4.0.0).
- **StringTie:** whether the CPAS 5-bp/20-read algorithm executes at all; whether `-s` ever acts as an
  unconditional monoexonic output threshold (source argues against); the origin/update of
  `CGroup::neg_prop`, hence the exact fallback strand split for unstranded mono reads.
- **All three:** no reference-flank internal-priming or PAS-motif predicate was found in the inspected
  paths. This is an absence-of-evidence claim scoped to the files read, listed above per tool.

---

## 8. Sources

Source trees read directly at the commits in the provenance table; file paths and symbol names are
given inline throughout so each claim can be re-checked.

1. Prjibelski AD, et al. "Accurate isoform discovery with IsoQuant using long reads."
   *Nature Biotechnology* 41, 915–918 (2023). doi:10.1038/s41587-022-01565-y —
   https://pmc.ncbi.nlm.nih.gov/articles/PMC10344776/
2. Chen Y, et al. "Context-aware transcript quantification from long-read RNA-seq data with Bambu."
   *Nature Methods* 20, 1187–1195 (2023). doi:10.1038/s41592-023-01908-w —
   https://pmc.ncbi.nlm.nih.gov/articles/PMC10448944/
3. Shinder I, et al. StringTie3 (2026); methods identify v3.0.3 —
   https://pmc.ncbi.nlm.nih.gov/articles/PMC13250942/
4. Shumate A, et al. "Improved transcriptome assembly using a hybrid of long and short reads with
   StringTie." *PLoS Computational Biology* (2022) —
   https://pmc.ncbi.nlm.nih.gov/articles/PMC9191730/
5. Kovaka S, et al. "Transcriptome assembly from long-read RNA-seq alignments with StringTie2."
   *Genome Biology* 20, 278 (2019) — https://pmc.ncbi.nlm.nih.gov/articles/PMC6912988/
