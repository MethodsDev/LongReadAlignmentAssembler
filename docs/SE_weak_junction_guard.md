# Single-exon genes hidden behind a weak junction

## The defect (v0.47.0)

Some single-exon genes are reported only as a spliced 2-exon model, with no single-exon
model, although nearly every read is unspliced and full length. TUSCO single-exon genes
in the bulk PacBio cell lines, LRAA de novo:

| sample | gene | reads | reads with the model's intron (J) | unspliced reads through the intron (T) | T/J |
|---|---|---|---|---|---|
| BT474 | RNF26 | 1,119 | 2 | 959 | 480 |
| BT474 | SEPHS2 | 1,446 | 2 | 1,201 (1,217 weighted) | 608 |
| BT474 | C6orf47 | 388 | 2 | 300 | 150 |
| BT474 | NCBP2AS2 | 1,241 | 2 | 1,216 | 614 |
| HG002 | SEPHS2 | 810 | 2 | 423 | 211 |
| HG002 | NCBP2AS2 | 889 | 5 | 807 | 161 |
| K562 | MARS2 | 155 | 2 | 123 | 61 |
| UHRR | NCBP2AS2 | 672 | 2 | 631 | 315 |

Two SE-pass steps both treat the spliced model as authoritative, whatever its support:

1. **Encapsulation mask** (`min_SE_read_ME_exon_overlap_pct`, 50). An unspliced read
   is withheld from the SE graph when half of it lies in one exon of a spliced model.
   When one exon of the 2-exon model covers most of the gene (SEPHS2, C6orf47, MARS2,
   and HG002 NCBP2AS2 with its 51 bp intron), every unspliced read is withheld.
2. **Intronic background subtraction** (`SE_subtract_intronic_background`, v0.35.0). The
   median SE coverage over the spliced model's intron is subtracted as a floor. At RNF26
   and NCBP2AS2 the unspliced reads survive the mask, so the "floor" is the read-through
   pile itself (956 and 1,217) and subtracting it zeroes the intron, which splits the gene.

Minigenome toggles (v0.47.0, benchmark flags): disabling only the subtraction fixes RNF26
and BT474/UHRR NCBP2AS2; disabling only the mask leaves truncated SE fragments (the
subtraction then eats the unmasked pile); disabling both fixes all eight.

Before v0.35.0 the ME graph itself emitted monoexonic fragments cut at the weak junction's
splice sites (v0.34.0 RNF26: 119334527-119336055), and those joined the mask, hiding the
full-length reads. The ME monoexonic suppression removed that path and the v0.35.0
subtraction took over as the blocker.

## The guards

`build_ME_transcripts` records each model intron's spliced-read support J (the ME graph's
`Intron` count, normalization-weighted) on the model (`Transcript.set_intron_read_support`).
In the SE pass:

- `SE_mask_max_readthrough_to_junction_ratio` (default 20): T is the
  normalization-weighted count of unspliced reads whose span contains the whole intron.
  When T > ratio x J, the two exons flanking that intron do not mask unspliced reads.
  Other exons of the same model still mask.
- `SE_intronic_background_max_floor_to_junction_ratio` (default 20): an intron whose
  floor (the median it would subtract) exceeds ratio x J is not subtracted.

None disables a guard (v0.47.0 behaviour); unknown J never triggers one.
`SE_weak_junction_stats_path` appends a per-intron diagnostics table (J, T, floor) for tuning.

## Choosing the default (BT474 chr19, guards off)

Model introns 11,257 (T > 0 at 663). T/J > 1: 79; > 10: 21; > 20: 13; > 50: 3. No intron
with T/J > 10 is a GENCODE v39 intron. The highest annotated introns are ZNF625 (9.3),
FOSB's retained introns (3.0-3.5) and CACTIN (3.5). The TUSCO failures sit at 61-614.
Floor/J over the 510 introns with a non-zero floor: p99 3, max 11, and no annotated intron
above 8 (BLOC1S3). 20 leaves a 2x margin above every annotated intron and a 3x margin
below the failures.

## Effect

- Minigenomes: all eight cases above give the full-length single-exon model (RNF26 1,057
  reads; SEPHS2 1,432; C6orf47 372; NCBP2AS2 1,217). The weak spliced model is dropped
  at four of them and kept with 2-5 reads at MARS2 and HG002 NCBP2AS2. Nothing else in
  the minigenome regions changes.
- BT474 chr19, defaults (20): the mask guard exempts 15 model-intron pairs and the floor
  guard skips 1 intron; every v0.47.0 model is kept unchanged and one single-exon model is
  gained, at the annotated single-exon gene ENSG00000261526 (87 reads). At ratio 10: 4
  single-exon models gained, 0 lost. gffcompare against GENCODE is unchanged at the base,
  exon, intron and intron-chain levels (transcript precision 31.0 -> 31.1).

## Caveat

Tuned on bulk PacBio (mostly cytoplasmic RNA). Nuclear libraries carry much more
read-through, and the subtraction was designed on that kind of data, so measure the
ratio distribution on one nuclear chromosome before relying on the default there.
