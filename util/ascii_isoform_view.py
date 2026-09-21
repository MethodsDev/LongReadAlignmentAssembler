#!/usr/bin/env python3
# encoding: utf-8

"""Draw isoform structures from one or more GTFs as stacked ASCII rows.

Made for eyeballing structural questions that coordinates alone answer badly:
did LRAA reconstruct the annotated model or a variant of it, which junction is
novel, where does a multipath read path stop short of the full-length isoform.

Examples
--------
  # every model at a gene, reference and LRAA output overlaid
  util/ascii_isoform_view.py -g ref.gtf=# -g LRAA.gtf -G SIRV5

  # compare against one named model: rows get a structural verdict and their
  # non-reference introns are drawn with '*'
  util/ascii_isoform_view.py -g ref.gtf=# -g LRAA.gtf -G SIRV5 --ref-transcript SIRV503

  # a genomic window at true scale, with two positions marked
  util/ascii_isoform_view.py -g LRAA.gtf -r chr17:7,668,000-7,688,300 \\
      --mode proportional --mark 7675053 --mark 7674290

  # LRAA multipath debug dump; those records carry no strand
  util/ascii_isoform_view.py -g __mpgns.pre.gtf --unstranded-as + -r contig:900-4200

  # add SQANTI-like categories, assigned by pylib/SQANTI_like_annotator.py
  util/ascii_isoform_view.py -g LRAA.gtf -G SIRV5 --sqanti ref.gtf
"""

import argparse
import logging
import os
import re
import sys

PYLIB_DIR = os.path.join(os.path.dirname(os.path.realpath(__file__)), "../pylib")
if PYLIB_DIR not in sys.path:
    sys.path.insert(0, PYLIB_DIR)

from Ascii_isoform_illustrator import (  # noqa: E402
    Track,
    View,
    compare_structures,
    exons_of,
    transcripts_from_gtf,
)

logger = logging.getLogger(__name__)

GLYPH_CYCLE = ["=", "#", "o", "+", "~", "x"]


def parse_region(region_str):
    match = re.match(r"^(\S+):([\d,]+)-([\d,]+)$", region_str)
    if match is None:
        raise argparse.ArgumentTypeError(
            "region must look like chr:lend-rend, got {}".format(region_str)
        )
    contig_acc = match.group(1)
    lend = int(match.group(2).replace(",", ""))
    rend = int(match.group(3).replace(",", ""))
    if rend < lend:
        raise argparse.ArgumentTypeError("region rend < lend: {}".format(region_str))
    return contig_acc, lend, rend


def parse_gtf_spec(spec):
    """'path' or 'path=GLYPH' -> (path, glyph or None)."""
    if "=" in spec:
        path, glyph = spec.rsplit("=", 1)
        if len(glyph) != 1:
            raise argparse.ArgumentTypeError(
                "glyph must be one character: {}".format(spec)
            )
        return path, glyph
    return spec, None


def main():
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument(
        "-g",
        "--gtf",
        action="append",
        required=True,
        metavar="PATH[=GLYPH]",
        help="GTF to draw; repeatable. Optional '=GLYPH' sets the exon character.",
    )
    parser.add_argument(
        "-r", "--region", help="restrict to chr:lend-rend (transcripts must be contained)"
    )
    parser.add_argument("-G", "--gene", help="restrict to a gene_id / gene_name (substring ok)")
    parser.add_argument(
        "-t",
        "--transcript",
        action="append",
        help="restrict to a transcript_id (substring ok); repeatable",
    )
    parser.add_argument(
        "--ref-transcript",
        help="transcript_id to compare every other row against (substring ok)",
    )
    parser.add_argument(
        "--sqanti",
        metavar="REF_GTF",
        help="annotate each row with a SQANTI-like category from pylib/SQANTI_like_annotator.py",
    )
    parser.add_argument(
        "--annot-key",
        action="append",
        default=[],
        metavar="KEY",
        help="GTF attribute to show beside each row (e.g. TPM); repeatable",
    )
    parser.add_argument("-w", "--width", type=int, default=100, help="drawing width in columns")
    parser.add_argument(
        "--mode",
        choices=["compressed", "proportional"],
        default="compressed",
        help="compressed gives introns a fixed 4 columns; proportional is true bp scale",
    )
    parser.add_argument("--max-ticks", type=int, default=8, help="max coordinate ticks")
    parser.add_argument("--color", action="store_true", help="ANSI color")
    parser.add_argument(
        "--mark", action="append", type=int, default=[], help="mark a genomic position; repeatable"
    )
    parser.add_argument(
        "--max-transcripts",
        type=int,
        default=40,
        help="stop after this many rows (0 = no limit)",
    )
    parser.add_argument(
        "--unstranded-as",
        choices=["+", "-"],
        help="strand to assume for records whose GTF strand is not +/- (e.g. __mpgns.*.gtf)",
    )
    parser.add_argument("--quiet", action="store_true", help="suppress parser logging")

    args = parser.parse_args()

    logging.basicConfig(
        level=logging.ERROR if args.quiet else logging.WARNING,
        format="%(levelname)s %(module)s: %(message)s",
    )

    contig_acc = lend = rend = None
    if args.region:
        contig_acc, lend, rend = parse_region(args.region)

    gtf_specs = [parse_gtf_spec(spec) for spec in args.gtf]
    # an explicitly requested glyph wins; auto-assignment then avoids it, so two
    # GTFs never end up indistinguishable on the page
    taken = {glyph for _, glyph in gtf_specs if glyph is not None}
    available = [g for g in GLYPH_CYCLE if g not in taken]
    for i, (gtf_path, glyph) in enumerate(gtf_specs):
        if glyph is None:
            gtf_specs[i] = (gtf_path, available.pop(0) if available else GLYPH_CYCLE[-1])

    entries = []  # (Track, Transcript)
    seen_names = set()
    for gtf_path, glyph in gtf_specs:
        transcripts = transcripts_from_gtf(
            gtf_path,
            contig_acc=contig_acc,
            lend=lend,
            rend=rend,
            gene=args.gene,
            transcript_ids=args.transcript,
            unstranded_as=args.unstranded_as,
        )
        if not transcripts:
            logger.warning("no transcripts selected from %s", gtf_path)
        source_tag = os.path.basename(gtf_path)
        for transcript in transcripts:
            name = transcript.get_transcript_id()
            if name in seen_names:
                name = "{}@{}".format(name, source_tag)
            seen_names.add(name)
            meta = transcript.get_meta()
            annot_bits = [
                "{}={}".format(key, meta[key]) for key in args.annot_key if key in meta
            ]
            entries.append(
                (
                    Track(
                        name,
                        exons_of(transcript),
                        transcript.get_strand(),
                        glyph=glyph,
                        annot=" ".join(annot_bits),
                    ),
                    transcript,
                )
            )

    if not entries:
        sys.exit("Error, no transcripts selected. Check --region / --gene / --transcript.")

    # One drawing is one coordinate axis.  Without this, a gene name or transcript
    # substring matching on two contigs silently overlays unrelated loci, and the
    # rows line up on the page while sharing nothing in the genome.
    selected_contigs = sorted({transcript.get_contig_acc() for _, transcript in entries})
    if len(selected_contigs) > 1:
        sys.exit(
            "Error, selection spans {} contigs ({}); one drawing is one coordinate "
            "axis. Narrow it with --region chr:lend-rend, or a gene/transcript "
            "unique to one contig.".format(
                len(selected_contigs), ", ".join(selected_contigs[:5])
            )
        )

    entries.sort(key=lambda pair: (pair[0].lend, pair[0].rend, pair[0].name))
    truncated = 0
    if args.max_transcripts and len(entries) > args.max_transcripts:
        truncated = len(entries) - args.max_transcripts
        entries = entries[: args.max_transcripts]

    if args.sqanti:
        _annotate_sqanti(entries, args.sqanti)

    ref_track = None
    if args.ref_transcript:
        candidates = [t for t, _ in entries if t.name == args.ref_transcript]
        if not candidates:
            candidates = [t for t, _ in entries if args.ref_transcript in t.name]
        if not candidates:
            sys.exit(
                "Error, --ref-transcript {} not among the selected transcripts".format(
                    args.ref_transcript
                )
            )
        ref_track = candidates[0]

    tracks = [track for track, _ in entries]
    if ref_track is not None:
        for track in tracks:
            if track is ref_track:
                track.annot = (track.annot + " " if track.annot else "") + "(reference)"
                continue
            comparison = compare_structures(track.exons, ref_track.exons)
            track.annot = (track.annot + " " if track.annot else "") + comparison.summary()
            for intron in comparison.query_only_introns:
                track.junc_chars[intron] = "*"

    view_contig = contig_acc or entries[0][1].get_contig_acc()
    if args.region:
        view_lend, view_rend = lend, rend
    else:
        view_lend = min(t.lend for t in tracks)
        view_rend = max(t.rend for t in tracks)
        pad = max(1, int((view_rend - view_lend + 1) * 0.02))
        view_lend, view_rend = view_lend - pad, view_rend + pad

    view = View(
        view_contig,
        view_lend,
        view_rend,
        width=args.width,
        mode=args.mode,
        color=args.color,
        max_ticks=args.max_ticks,
    )
    if ref_track is not None:
        view.note("'*' marks an intron absent from {}".format(ref_track.name))
    for track in tracks:
        view.add(track)
    for pos in args.mark:
        view.mark(pos, "^")  # the legend line prints the coordinate itself
    if truncated:
        view.note("{} further transcript(s) not shown (--max-transcripts)".format(truncated))

    print(view.render())

    return 0


def _annotate_sqanti(entries, ref_gtf):
    """Append the SQANTI-like category assigned by the repo's own annotator."""
    from SQANTI_like_annotator import SQANTI_like_annotator

    annotator = SQANTI_like_annotator(ref_gtf)
    for track, transcript in entries:
        class_info = annotator.classify_alignment_or_isoform(
            transcript.get_contig_acc(),
            transcript.get_strand(),
            transcript.get_transcript_id(),
            transcript,
        )
        category = class_info.get("sqanti_cat", "?")
        matching = class_info.get("matching_isoforms", "")
        bit = category if not matching else "{}:{}".format(category, matching)
        track.annot = (track.annot + " " if track.annot else "") + bit


if __name__ == "__main__":
    sys.exit(main())
