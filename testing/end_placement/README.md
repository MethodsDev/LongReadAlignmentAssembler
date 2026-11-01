# end_placement

Execution test for where LRAA puts isoform ends (v0.48.0). Five minigenomes cut from the BT474
PacBio Kinnex benchmark library (`util/misc/extract_contig_region_inputs.py`, GRCh38, GENCODE v39
loci kept whole). Coordinates in the checks are local: genomic = local + region start - 1.

| locus | region (GRCh38) | what it guards | v0.47.x result |
|---|---|---|---|
| `TIGD7` | chr16:3,262,000-3,320,000 | a 5' end extended past its exon stays at the reads' position (43431) inside the node it reaches; TSS-site form at 43728 kept | dominant form at 43820, site form moved to 43975 |
| `TIGD7_merge` | same; `data/TIGD7.merge_input.c{1,2}.gtf` are v0.47.0 runs on the two halves of the reads | merged models keep their input models' 5' ends (merge uses the normal end logic, one end per input model) | ends moved to node edges (43975, 43975, 43212) |
| `SPACER_revert` | chr19:14,071,536-14,207,187 | an extension onto an exon segment the model's 5' exon has no edge to is reverted (EXTENSION_REVERTED), and the model is kept | no extension (5' end at the exon edge) |
| `SEPHS2` | chr16:30,438,631-30,450,874 | single-exon TUSCO gene beside a 2-read junction is reported as one exon with its ~1,400 reads (SE weak-junction guard) | spliced model only, 213 reads |
| `COMT` | chr22:19,874,517-20,071,569 | MB-COMT and the internal-promoter S-COMT (an ISM of it) both called with a TSS | same (regression guard) |

`data/cap_calib.json` sets BT474's TSS cap calibration (signature on, baseline 1): the
minigenomes are too small to calibrate.

    make test     # runs and checks all five, about a minute
