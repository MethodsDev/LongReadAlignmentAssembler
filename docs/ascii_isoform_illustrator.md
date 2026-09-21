# ASCII isoform structure views

`pylib/Ascii_isoform_illustrator.py` draws transcript structures as one line of text
each, so several models stack and their differences line up column-wise.
`util/ascii_isoform_view.py` is the command-line front end over GTFs.

Use it when the question is structural — did LRAA rebuild the annotated model, which
junction is novel, where does a read's multipath stop short — and reading exon
coordinate lists is the slow way to answer.

Drawing conventions follow `PASApipeline/PerlLib/Ascii_genome_illustrator.pm`:
right-aligned coordinate label, single-character exon glyph, arrowhead at the 3' end
(`>` at the high coordinate for `+`, `<` at the low for `-`), name to the right.

## Command line

```bash
# reference and LRAA output overlaid at one gene
util/ascii_isoform_view.py -g ref.gtf=# -g LRAA.gtf -G SIRV5

# compare every row against one model: each gets a verdict, and its introns that
# the reference does not have are drawn '*'
util/ascii_isoform_view.py -g ref.gtf=# -g LRAA.gtf -G SIRV5 --ref-transcript SIRV508
```

```
SIRV5:784-12,084   [intron-compressed, 100 cols]
  '*' marks an intron absent from SIRV508
                                   2,488      3,643      5,627          6,723     7,307   ...
                       |             |          |          |              |         |
           1009-10991+     #----#####-----#----##----#----#----#----#----#----#----##--->  SIRV508      (reference)
           1009-10991+     =----=====-----=----==----=----=----=----=*********=----==--->  ...iso-3     novel_junctions introns shared=14 query_only=1 ref_only=2
           1009-10991+     =----=====-----=----==----=----=----=----=----=----=----==--->  ...iso-4     identical_splice_pattern introns shared=16 query_only=0 ref_only=0
            1009-2398+     #----#####>                                                     SIRV506      contained_subchain introns shared=1 query_only=0 ref_only=15
```

Useful flags:

| flag | effect |
|---|---|
| `-g PATH[=GLYPH]` | GTF to draw, repeatable; glyphs auto-assigned and never collide |
| `-r chr:lend-rend` | window; selects every transcript **overlapping** it, drawing clipped |
| `-G`, `-t` | restrict by gene_id/gene_name, or by transcript_id (substring ok) |
| `--ref-transcript ID` | comparison anchor; adds verdicts and `*` marks |
| `--sqanti REF.gtf` | adds the category from `pylib/SQANTI_like_annotator.py` (FSM/ISM/NIC/NNIC/...) |
| `--annot-key KEY` | show a GTF attribute per row (`TPM`, `TSS`, `PolyA_read_count`, ...) |
| `--mode proportional` | true bp scale; default `compressed` gives each intron 4 columns |
| `--unstranded-as +` | needed for `__mpgns.*.gtf`, whose records carry strand `?` |
| `--mark POS` | annotate a genomic position under the drawing |

`compressed` is the default because at true scale a 100 bp exon inside a 100 kb locus
is under one column: the picture becomes a row of dashes. Switch to `proportional`
when relative spacing is the thing being shown.

A window keeps every **overlapping** transcript and clips the drawing, unlike
`GTF_contig_to_transcripts.parse_GTF_to_Transcripts(lend_restrict=...)`, which keeps
only fully contained ones. Containment is right for chunked runs, where a truncated
isoform would be indistinguishable from a real one; it is wrong for a viewer, where it
would hide the long models you are looking for. A model whose 3' end lies outside the
window is drawn without its arrowhead, so a clipped end is never mistaken for a real
terminus.

One drawing is one contig. A selection spanning contigs is refused rather than
overlaid, because a single coordinate axis cannot honestly hold two loci: the rows
would line up on the page while sharing nothing in the genome. Narrow with
`--region`, or with a gene/transcript unique to one contig.

Exon structures come from `GTF_contig_to_transcripts.parse_GTF_to_Transcripts`, the
same reader the pipeline uses, so what you see is what LRAA sees — including its
`read_aln_gap_merge_int` rule, which fuses exons separated by a gap of 10 bp or less
(`Transcript.__init__`). A sub-10 bp gap in the input GTF therefore does not appear
as a junction here, and it does not appear as one anywhere else in LRAA either.

## From Python

```python
import sys; sys.path.insert(0, "pylib")
from Ascii_isoform_illustrator import illustrate, transcripts_from_gtf

models = transcripts_from_gtf("LRAA.gtf", contig_acc="SIRV5")
print(illustrate(models, reference="t:SIRV5:+:comp-1:iso-4", contig_acc="SIRV5"))
```

`illustrate()` accepts LRAA `Transcript` objects, `Track` objects, or
`(name, exons, strand[, annot])` tuples, so ad-hoc structures mix with parsed ones.

For finer control build a `View` and add rows yourself:

```python
from Ascii_isoform_illustrator import View, exons_of

view = View("SIRV5", 1000, 11000, width=100, mode="compressed")
view.rule("annotation")
view.add_transcript("SIRV508", ref_exons, "+", glyph="#")
view.rule("reconstructed")
view.add_lraa_transcript(transcript_obj, glyph="=", annot="TPM=12.4")
view.mark(7675053, "^", "novel donor")
print(view.render())
```

## Simple paths through the splice graph

An LRAA simple path is a list of splice-graph node ids (`E:12`, `I:7`, `TSS:3`,
`POLYA:1`, `SPACER`), so it can only be resolved with the `Splice_graph` that issued
those ids — in-process, during a run or in a debugger:

```python
from Ascii_isoform_illustrator import track_from_simple_path

view.add(track_from_simple_path(splice_graph, mpgn.get_simple_path(), "mp4x"))
```

Adjacent exon nodes are fused into one drawn exon (the splice graph splits an exon at
every internal splice site, and drawing those splits would invent junctions). A
`SPACER` — a gap the path never resolved — is drawn `?` rather than as a confident
intron. `TSS`/`PolyA` nodes appear as `T`/`A` marks on the row.

Outside a live run, the multipath debug dumps are already GTF: draw
`__mpgns.*.gtf` with `--unstranded-as`.

## Structural verdicts

`compare_structures(query_exons, ref_exons)` reports junction-level facts only:

| relationship | meaning |
|---|---|
| `identical_splice_pattern` | same intron chain |
| `contained_subchain` | query introns are a contiguous run of the reference's |
| `contains_reference_subchain` | the reference's introns are a contiguous run of the query's |
| `known_junctions_recombined` | every query intron is in the reference, but not as a run |
| `novel_junctions` | at least one query intron is absent from the reference |
| `no_shared_junctions` | both spliced, no intron in common — a junction statement, not a locus one: a single exon skip can land here while overlapping heavily, so read `overlap_bp` for position |
| `monoexonic_*` | one or both sides unspliced; reports overlap bp instead |

Introns use the `Transcript.get_introns()` convention, `(prev_exon_rend + 1,
next_exon_lend - 1)`, which is the same identity `SQANTI_like_annotator` keys on.
This is deliberately not a second SQANTI taxonomy — pass `--sqanti` to get categories
from the annotator that owns them.

Tests: `pylib/test_ascii_isoform_illustrator.py`.
