# Why an FSM isoform can lack FSM (or unique-FSM) read support

A reconstructed isoform can be categorized **FSM** (its intron chain reproduces a reference
transcript's exactly) yet report **`has_FSM_read = 0`** and/or **`uniq_FSM_reads = 0`**. This is
not a bug and not a contradiction — it follows from how the two quant columns are defined and how
reads are assigned to isoforms. This note explains the definitions, the two mechanisms that
produce it, and two worked examples traced to the read.

## The two quant columns, and the per-read flags behind them

For each `(read, isoform)` pair LRAA computes two booleans (`Quantify.py`):

- `mp_is_FSM` — the read's intron chain **equals** the isoform's intron chain (reproduces it
  exactly; not a prefix, not a subset, not merely compatible).
- `mp_is_unique` — `n_compat_mp == 1`, i.e. the read is compatible with **exactly one** isoform.

and reports two isoform-level columns in `*.quant.expr`:

- `has_FSM_read` — TRUE if **any** assigned read is FSM to the isoform, exclusivity aside.
- `uniq_FSM_reads` — count of reads that are **both** unique **and** FSM
  (increment guarded by `if mp_is_unique: if mp_is_FSM:`).

The per-read `is_unique` and `is_FSM` flags are also written to `*.quant.tracking.gz`, one row per
`(read, isoform)` assignment, so both columns are recomputable from the tracking.

**Crucial:** these flags are only computed for reads that are **assigned** to the isoform. A read
that fails read→isoform assignment never scores against it, whatever its chain.

Downstream analysis (e.g. the SQANTI-like read-support plots) keys "FSM read support" on
`uniq_FSM_reads >= 1` — the strict *unique AND FSM* reading. `has_FSM_read` is the looser reading
and is the column that distinguishes the two mechanisms below. Its comment in `Quantify.py` says
so: `uniq_FSM_reads` "cannot distinguish 'no read traverses this chain' from 'reads do, but each
also fits another model'."

## Two mechanisms for `uniq_FSM_reads == 0`

### A. `has_FSM_read == 1` — FSM reads exist but none is unique

Reads reproduce the chain, but each is FSM to **more than one** isoform, so none is unique and
`uniq_FSM_reads` stays 0. This is the common case, and it is a direct consequence of
**splice-pattern collapse**: several isoforms can share one intron chain while differing only at
the TSS/PolyA, and a read reproducing that chain is FSM to all of them at once.

Collapsing by splice pattern pools those siblings into one feature, so the same reads become
unique to the single collapsed chain and `uniq_FSM_reads` recovers. Collapse therefore *reduces*
the count of unsupported FSM features rather than causing it.

**Worked example — GLIPR1 (`ENST00000266659`).** LRAA emitted three FSM models sharing one
splice hashcode (identical intron chain, differing only in 3′ end). One of them had 1,545 reads
that reproduce the chain but **all** with `is_unique = 0`; each was also assigned to the two
siblings (a single read was `is_FSM = 1` for all three at once). Its unique reads were the
partial, non-FSM ones — so `has_FSM_read = 1`, `uniq_FSM_reads = 0`. Summing the shared pattern
over the three models restored positive unique-FSM support.

### B. `has_FSM_read == 0` — the full-length read exists but isn't assigned to the isoform

Here a molecule reproducing the exact chain **does exist**, but it is never **assigned** to the
isoform, so it never scores as FSM. The dominant cause found so far is a **terminal read-through**
interacting with the assignment overlap gate.

Read→isoform assignment requires at least `fraction_read_align_overlap` of the read's aligned
length to overlap the isoform's structure (`LRAA_Globals.py`, default **0.75**), applied on every
test of the compatibility cascade. An isoform's TSS/PolyA termini are set to the **majority**
reads' termini. A **minority** molecule that alone carries the full intron chain but reads through
one of those boundaries — extending a terminal exon past the isoform's inferred TSS or PolyA into
space no isoform models — can have >25% of its length outside the isoform, failing the 0.75 gate.
It is compatible in introns but rejected from assignment, so it is left unassigned and never
written to the tracking; the isoform reports `has_FSM_read = 0`.

**Worked example — OGA (`ENST00000370094`, 9 introns; chr10).** Exactly one read carried the
isoform's exact 9-intron chain (99.93% identity, MAPQ 60, from a retained cell). It survived
strand-splitting and pre-assembly normalization (present in the normalized assembly BAM), so it
was available when the isoform was built — but it is absent from every `quant.tracking`. Its exon
blocks vs the isoform:

```
                3' terminal exon        5' terminal exon        9 introns
read            13,048 - 14,155 (1108)  33,825 - 33,893 (69)    identical
isoform iso-4   13,936 - 14,155 ( 220)  33,825 - 34,443 (619)   identical
```

The read reads through **~888 bp past the isoform's inferred 3′/PolyA terminus** (into a region
no isoform models). That overhang is ~32% of the read's ~2,787 bp aligned length, so its overlap
with the isoform is ~68% — below the 0.75 gate. It is therefore not assigned to the isoform (nor
to anything, since nothing models its 3′ extension), so `has_FSM_read = 0` even though a molecule
carrying the exact chain exists. This was reproduced in a controlled single-gene minigenome run at
the same LRAA version.

Stated generally: an FSM isoform is real and correctly reconstructed, but its **own** defining
full-length molecule can be excluded from its quantification when that molecule reads through a
terminus the isoform did not adopt.

## Note on assembly mode

These observations are under **collapse mode** (`restrict_asm_to_collapse = True`, the default),
which disables chaining of overlapping reads — `MultiPathGraph.py:340` draws an overlap edge
between two compatible read-path nodes only when that flag is False. So a reconstructed isoform is
some read's structure, not a path stitched from several partial reads; situation B is therefore
about read **assignment**, not about the chain being assembled from fragments.

## Reproducing / investigating a specific case

Use the single-gene mechanism deep-dive (skill `lraa-reconstruction-eval`): BAM intron-chain scan
→ `quant.tracking` lookup → minigenome extraction (`util/misc/extract_contig_region_inputs.py`) +
`LRAA --debug` at the data's version, tracing the read across `*.chunked_work/` intermediates
(`chunk.strand.<±>.bam` → `chunk.<±>.norm.bam` → chunk `quant.tracking.gz`) to find the first
stage where it disappears.
