#!/usr/bin/env python3
# encoding: utf-8

"""
Ascii_isoform_illustrator - one-row ASCII drawings of transcript structures.

Comparing isoform structures (a reference GTF model against an LRAA-reconstructed
one, or against the multipath simple path a read followed through the splice graph)
is coordinate arithmetic that is tedious to do in your head and immediate to see.
This renders each structure as a single line of text, so structures stack and the
differences line up column-wise.

Drawing conventions follow PASApipeline/PerlLib/Ascii_genome_illustrator.pm:
right-aligned coordinate label, single-character glyphs, a transcription-direction
arrowhead at the 3' end ('>' at the high coordinate for '+', '<' at the low for
'-'), and the feature name to the right of the drawing.

Coordinates
-----------
1-based inclusive genomic coordinates throughout, matching Transcript and GTF.
Introns are (prev_exon_rend + 1, next_exon_lend - 1), i.e. exactly what
Transcript.get_introns() returns, so intron identity here is the same identity
SQANTI_like_annotator uses.

Scaling
-------
mode="compressed" (default) gives exonic sequence proportional space and every
intron a fixed small number of columns.  At locus scale that is the difference
between seeing the exon structure and seeing a row of dashes: a 100 bp exon in a
100 kb locus is under one column at true scale.  Use mode="proportional" when
true spacing is the point (e.g. showing that an intron is enormous).

Structural comparison
---------------------
compare_structures() reports junction-level facts only: which introns are shared,
which are unique to each model, and how the two chains nest.  It deliberately does
NOT assign SQANTI categories -- SQANTI_like_annotator.py owns that taxonomy, and a
second one would drift from it.  util/ascii_isoform_view.py --sqanti calls the real
annotator when categories are wanted.
"""

from dataclasses import dataclass, field
from typing import Optional

from LRAA_Globals import SPACER
from Transcript import Transcript, GTF_contig_to_transcripts

__all__ = [
    "Track",
    "View",
    "StructureComparison",
    "compare_structures",
    "intron_chain",
    "exons_of",
    "transcripts_from_gtf",
    "exons_from_simple_path",
    "track_from_simple_path",
    "illustrate",
]


SPACER_CHAR = "?"  # intron glyph for an unresolved (SPACER) gap in a simple path

_ANSI = {
    "green": "32",
    "red": "31",
    "yellow": "33",
    "blue": "34",
    "magenta": "35",
    "cyan": "36",
    "grey": "90",
}


def _paint(text, color):
    code = _ANSI.get(color)
    return text if code is None else "\033[{}m{}\033[0m".format(code, text)


##############
# track model
##############


@dataclass
class Track:
    """One drawable row: exon blocks joined by introns."""

    name: str
    exons: list  # [(lend, rend), ...] 1-based inclusive
    strand: str = "."  # '+', '-', or '.' (no arrowhead)
    glyph: str = "="  # exon character
    intron: str = "-"  # default intron character
    annot: str = ""  # free text, right of the name
    span: Optional[tuple] = None  # envelope drawn as '.' outside the exons
    color: Optional[str] = None  # ANSI color name, honored when View(color=True)
    junc_chars: dict = field(default_factory=dict)  # (intron_lend, intron_rend) -> char
    site_marks: tuple = ()  # ((pos, char), ...) painted over the row

    def __post_init__(self):
        if len(self.glyph) != 1 or self.glyph.isspace():
            raise ValueError("glyph must be a single non-whitespace character")
        self.exons = sorted((min(a, b), max(a, b)) for a, b in self.exons)
        if not self.exons:
            raise ValueError("track {} has no exons".format(repr(self.name)))

    @property
    def lend(self):
        return self.exons[0][0]

    @property
    def rend(self):
        return self.exons[-1][1]

    def introns(self):
        return intron_chain(self.exons)


def intron_chain(exons):
    """Exon blocks -> intron coordinate pairs, Transcript.get_introns() convention."""
    exons = sorted(exons)
    return [
        (exons[i][1] + 1, exons[i + 1][0] - 1)
        for i in range(len(exons) - 1)
        if exons[i + 1][0] - 1 >= exons[i][1] + 1
    ]


def exons_of(transcript):
    """Exon blocks of an LRAA Transcript, sorted."""
    return sorted(tuple(seg) for seg in transcript.get_exon_segments())


####################
# coordinate mapping
####################


class _Scale:
    """Maps a genomic position to a column. Modes: proportional, compressed."""

    def __init__(self, lend, rend, width, mode="compressed", exonic=(), intron_cols=4):
        if rend < lend:
            raise ValueError("rend < lend")
        if width < 10:
            raise ValueError("width must be >= 10")
        self.lend, self.rend, self.width, self.mode = lend, rend, width, mode

        if mode == "proportional":
            self.segments = [(lend, rend, 0, width)]
            return
        if mode != "compressed":
            raise ValueError("mode must be 'proportional' or 'compressed'")

        # merge the exonic intervals of every track, clipped to the window
        merged = []
        for a, b in sorted((max(a, lend), min(b, rend)) for a, b in exonic):
            if b < a:
                continue
            if merged and a <= merged[-1][1] + 1:
                merged[-1][1] = max(merged[-1][1], b)
            else:
                merged.append([a, b])
        if not merged:
            self.mode = "proportional"
            self.segments = [(lend, rend, 0, width)]
            return

        # alternating blocks: [gap] exon [gap] exon ... [gap]
        blocks = []  # (start, end, is_exonic)
        cur = lend
        for a, b in merged:
            if a > cur:
                blocks.append((cur, a - 1, False))
            blocks.append((a, b, True))
            cur = b + 1
        if cur <= rend:
            blocks.append((cur, rend, False))

        n_gaps = sum(1 for _, _, ex in blocks if not ex)
        exon_bp = sum(b - a + 1 for a, b, ex in blocks if ex)
        gap_budget = min(n_gaps * intron_cols, max(0, width - 2 * len(merged)))
        exon_budget = width - gap_budget
        if exon_budget < len(merged):  # too many blocks for the width; give up on it
            self.mode = "proportional"
            self.segments = [(lend, rend, 0, width)]
            return

        self.segments = []
        col = 0
        gap_each = gap_budget // n_gaps if n_gaps else 0
        for a, b, ex in blocks:
            if ex:
                w = max(1, round((b - a + 1) / exon_bp * exon_budget))
            else:
                w = gap_each
            w = min(w, width - col)
            if w > 0:
                self.segments.append((a, b, col, col + w))
            col += w
        if self.segments:  # absorb rounding into the last block
            a, b, c0, _ = self.segments[-1]
            self.segments[-1] = (a, b, c0, width)

    def col(self, pos):
        """Column for a genomic position, or None when outside the window."""
        if pos < self.lend or pos > self.rend:
            return None
        for a, b, c0, c1 in self.segments:
            if a <= pos <= b:
                if b == a:
                    return c0
                frac = (pos - a) / (b - a + 1)
                return min(c1 - 1, c0 + int(frac * (c1 - c0)))
        return None

    def clamp_col(self, pos):
        c = self.col(min(max(pos, self.lend), self.rend))
        return 0 if c is None else c


########
# view
########


class View:
    """Accumulate tracks, then render an ASCII picture of them."""

    def __init__(
        self,
        contig_acc,
        lend,
        rend,
        width=100,
        mode="compressed",
        intron_cols=4,
        label_width=21,
        name_width=24,
        color=False,
        max_ticks=8,
    ):
        self.contig_acc, self.lend, self.rend = contig_acc, lend, rend
        self.width, self.mode, self.intron_cols = width, mode, intron_cols
        self.label_width, self.name_width = label_width, name_width
        self.color, self.max_ticks = color, max_ticks
        self.tracks = []
        self.marks = []  # (pos, char, label)
        self.rules = []  # (track index, label)
        self.notes = []  # free text under the title

    ## building

    def add(self, track):
        self.tracks.append(track)
        return self

    def add_transcript(
        self,
        name,
        exons,
        strand=".",
        glyph="=",
        annot="",
        span=None,
        color=None,
        highlight=(),
        highlight_char="*",
        site_marks=(),
    ):
        """Add one structure.  `highlight` introns get `highlight_char` instead of '-'."""
        junc_chars = {tuple(intron): highlight_char for intron in highlight}
        return self.add(
            Track(
                name,
                list(exons),
                strand,
                glyph,
                annot=annot,
                span=span,
                color=color,
                junc_chars=junc_chars,
                site_marks=tuple(site_marks),
            )
        )

    def add_lraa_transcript(
        self, transcript, name=None, glyph="=", annot="", **kwargs
    ):
        """Add an LRAA Transcript object."""
        return self.add_transcript(
            name if name is not None else transcript.get_transcript_id(),
            exons_of(transcript),
            transcript.get_strand(),
            glyph=glyph,
            annot=annot,
            **kwargs
        )

    def note(self, text):
        self.notes.append(text)
        return self

    def mark(self, pos, char="^", label=""):
        self.marks.append((pos, char, label))
        return self

    def rule(self, text=""):
        """Section separator inserted before the next track added."""
        self.rules.append((len(self.tracks), text))
        return self

    ## rendering

    def _scale(self):
        exonic = [e for t in self.tracks for e in t.exons]
        return _Scale(
            self.lend, self.rend, self.width, self.mode, exonic, self.intron_cols
        )

    def _row(self, track, sc):
        cells = [" "] * self.width

        if track.span:
            a, b = track.span
            if b >= self.lend and a <= self.rend:
                for i in range(sc.clamp_col(a), sc.clamp_col(b) + 1):
                    cells[i] = "."

        # introns first, so exon glyphs overwrite them at shared columns
        for intron in track.introns():
            if intron[1] < self.lend or intron[0] > self.rend:
                continue
            ch = track.junc_chars.get(intron, track.intron)
            c0 = sc.clamp_col(max(intron[0], self.lend))
            c1 = sc.clamp_col(min(intron[1], self.rend))
            for i in range(c0, c1 + 1):
                if cells[i] in (" ", "."):
                    cells[i] = ch

        for a, b in track.exons:
            if b < self.lend or a > self.rend:
                continue
            for i in range(sc.clamp_col(max(a, self.lend)), sc.clamp_col(min(b, self.rend)) + 1):
                cells[i] = track.glyph

        for pos, ch in track.site_marks:
            c = sc.col(pos)
            if c is not None:
                cells[c] = ch

        # Arrowhead at the 3' end, but only when that end is actually in view: a
        # clipped model that stops at the window edge must not be drawn as if it
        # terminated there.  Its true span stays readable in the coordinate label.
        if track.strand == "+" and self.lend <= track.rend <= self.rend:
            cells[sc.clamp_col(track.rend)] = ">"
        elif track.strand == "-" and self.lend <= track.lend <= self.rend:
            cells[sc.clamp_col(track.lend)] = "<"

        label = "{}-{}".format(track.lend, track.rend)
        body = "".join(cells)
        if self.color and track.color:
            body = _paint(body, track.color)
        return "{:>{lw}}{} {} {:<{nw}} {}".format(
            label,
            track.strand,
            body,
            track.name,
            track.annot,
            lw=self.label_width,
            nw=self.name_width,
        ).rstrip()

    def _ruler(self, sc):
        ticks = [" "] * self.width
        labels = [" "] * (self.width + 24)
        if sc.mode == "compressed":
            pts = []
            for a, b, c0, c1 in sc.segments:
                if c1 - c0 >= 2:
                    pts += [a, b]
            pts = sorted(set(pts))
            if self.max_ticks and len(pts) > self.max_ticks:  # thin, keeping the ends
                step = len(pts) / self.max_ticks
                keep = {0, len(pts) - 1}
                keep |= {int(i * step) for i in range(self.max_ticks)}
                pts = [p for i, p in enumerate(pts) if i in keep]
        else:
            n = max(2, min(self.max_ticks or 8, self.width // 14))
            pts = [
                self.lend + round(i * (self.rend - self.lend) / (n - 1))
                for i in range(n)
            ]

        last_end = -2  # -2, not -1, so a label centered on column 0 still prints
        for p in pts:
            c = sc.col(p)
            if c is None:
                continue
            ticks[c] = "|"
            s = "{:,}".format(p)
            start = max(0, c - len(s) // 2)
            if start > last_end + 1 and start + len(s) < len(labels):
                labels[start : start + len(s)] = list(s)
                last_end = start + len(s)
        pad = " " * (self.label_width + 2)
        return [pad + "".join(labels).rstrip(), pad + "".join(ticks).rstrip()]

    def render(self, title=None):
        if not self.tracks:
            raise ValueError("no tracks to render")
        sc = self._scale()
        out = []
        head = title or "{}:{:,}-{:,}".format(self.contig_acc, self.lend, self.rend)
        mode = "intron-compressed" if sc.mode == "compressed" else "true scale"
        bp_per_col = (self.rend - self.lend + 1) / self.width
        out.append(
            "{}   [{}, {} cols{}]".format(
                head,
                mode,
                self.width,
                (
                    ", ~{:.0f} bp/col".format(bp_per_col)
                    if sc.mode == "proportional"
                    else ""
                ),
            )
        )
        for n in self.notes:
            out.append("  " + n)
        out += self._ruler(sc)

        rule_at = dict(self.rules)
        for i, t in enumerate(self.tracks):
            if i in rule_at:
                out.append("")
                if rule_at[i]:
                    out.append("  " + rule_at[i])
            out.append(self._row(t, sc))

        if self.marks:
            cells = [" "] * self.width
            for pos, ch, _ in self.marks:
                c = sc.col(pos)
                if c is not None:
                    cells[c] = ch
            out.append(" " * (self.label_width + 2) + "".join(cells).rstrip())
            for pos, ch, lbl in self.marks:
                # the coordinate is the legend; a label is extra text, not a repeat
                out.append(
                    "{:>{lw}}  {}{}".format(
                        pos, ch, " " + lbl if lbl else "", lw=self.label_width
                    )
                )
        return "\n".join(out)

    def __str__(self):
        return self.render()


#########################
# structural comparison
#########################


@dataclass
class StructureComparison:
    """Junction-level relationship between a query structure and a reference one."""

    relationship: str
    shared_introns: list
    query_only_introns: list
    ref_only_introns: list
    query_num_exons: int
    ref_num_exons: int
    overlap_bp: int

    def summary(self):
        if self.query_num_exons == 1 or self.ref_num_exons == 1:
            return "{} exons={}/{} overlap={}bp".format(
                self.relationship, self.query_num_exons, self.ref_num_exons, self.overlap_bp
            )
        return "{} introns shared={} query_only={} ref_only={}".format(
            self.relationship,
            len(self.shared_introns),
            len(self.query_only_introns),
            len(self.ref_only_introns),
        )


def _is_contiguous_subchain(sub, full):
    n = len(sub)
    if n == 0 or n > len(full):
        return False
    return any(full[i : i + n] == sub for i in range(len(full) - n + 1))


def _overlap_bp(exons_a, exons_b):
    total = 0
    for a0, a1 in exons_a:
        for b0, b1 in exons_b:
            lo, hi = max(a0, b0), min(a1, b1)
            if hi >= lo:
                total += hi - lo + 1
    return total


def compare_structures(query_exons, ref_exons):
    """Compare two exon structures at junction resolution.

    relationship is one of:
      identical_splice_pattern     - same intron chain
      contained_subchain           - query introns are a contiguous run of the reference's
      contains_reference_subchain  - the reference's introns are a contiguous run of the query's
      known_junctions_recombined   - every query intron is in the reference, but not as a run
      novel_junctions              - at least one query intron is absent from the reference
      no_shared_junctions          - spliced on both sides, nothing in common
      monoexonic_both              - neither side is spliced
      monoexonic_query             - query unspliced, reference spliced
      monoexonic_reference         - query spliced, reference unspliced
    """
    q_exons = sorted(tuple(e) for e in query_exons)
    r_exons = sorted(tuple(e) for e in ref_exons)
    q = intron_chain(q_exons)
    r = intron_chain(r_exons)
    r_set = set(r)
    q_set = set(q)
    shared = [i for i in q if i in r_set]
    q_only = [i for i in q if i not in r_set]
    r_only = [i for i in r if i not in q_set]
    overlap = _overlap_bp(q_exons, r_exons)

    if not q and not r:
        rel = "monoexonic_both"
    elif not q:
        rel = "monoexonic_query"
    elif not r:
        rel = "monoexonic_reference"
    elif q == r:
        rel = "identical_splice_pattern"
    elif _is_contiguous_subchain(q, r):
        rel = "contained_subchain"
    elif _is_contiguous_subchain(r, q):
        rel = "contains_reference_subchain"
    elif not shared:
        # ordered before the novel-junction test on purpose: a chain sharing
        # nothing with the reference is a different locus-level statement than a
        # chain that differs at one junction, and every disjoint chain trivially
        # has query-only introns, which would otherwise swallow this case whole
        rel = "no_shared_junctions"
    elif q_only:
        rel = "novel_junctions"
    else:
        rel = "known_junctions_recombined"

    return StructureComparison(
        relationship=rel,
        shared_introns=shared,
        query_only_introns=q_only,
        ref_only_introns=r_only,
        query_num_exons=len(q_exons),
        ref_num_exons=len(r_exons),
        overlap_bp=overlap,
    )


##########
# readers
##########


def transcripts_from_gtf(
    gtf_filename,
    contig_acc=None,
    lend=None,
    rend=None,
    strand=None,
    gene=None,
    transcript_ids=None,
    unstranded_as=None,
):
    """Parse a GTF into LRAA Transcript objects, optionally filtered.

    A region here selects every transcript OVERLAPPING [lend, rend]; the drawing
    clips to the window, so a model running off the edge is still visible and still
    recognisable.  That is deliberately not the containment rule
    parse_GTF_to_Transcripts enforces for lend_restrict/rend_restrict: containment
    exists so a chunked run never invents a truncated isoform, and applying it to a
    viewer would silently hide exactly the long models you are looking for.

    `gene` matches gene_id or gene_name, exactly or as a substring.
    `transcript_ids` matches transcript_id exactly or as a substring.
    `unstranded_as` supplies a strand for records whose GTF strand is not '+'/'-',
    which is what the LRAA multipath debug dumps (__mpgns.*.gtf) emit.
    """
    contig_to_transcripts = GTF_contig_to_transcripts.parse_GTF_to_Transcripts(
        gtf_filename,
        chr_restrict=contig_acc,
        strand_restrict=strand,
        strand_default=unstranded_as,
    )

    transcripts = []
    for contig, tlist in contig_to_transcripts.items():
        transcripts.extend(tlist)

    if lend is not None and rend is not None:
        transcripts = [
            t
            for t in transcripts
            if t.get_coords()[0] <= rend and t.get_coords()[1] >= lend
        ]

    if gene is not None:
        transcripts = [t for t in transcripts if _matches(gene, t.get_gene_id(), t.get_gene_name())]

    if transcript_ids:
        wanted = list(transcript_ids)
        transcripts = [
            t for t in transcripts if any(_matches(w, t.get_transcript_id()) for w in wanted)
        ]

    transcripts.sort(key=lambda t: (t.get_contig_acc(), t.get_coords()[0], t.get_coords()[1]))
    return transcripts


def _matches(query, *values):
    for v in values:
        if v is None:
            continue
        if query == v or query in v:
            return True
    return False


def exons_from_simple_path(splice_graph, simple_path):
    """LRAA simple path -> (exon blocks, spacer gaps, boundary site marks).

    Returns (exons, spacers, site_marks):
      exons      - merged exon blocks, adjacent splice-graph exon nodes fused
      spacers    - [(lend, rend), ...] unresolved SPACER gaps, drawn with '?'
      site_marks - [(pos, char), ...] 'T' for a TSS node, 'A' for a PolyA node
    """
    coord_blocks = []  # (kind, lend, rend); kind in {'exon', 'intron', 'spacer'}
    site_marks = []
    for node_id in simple_path:
        if node_id == SPACER:
            coord_blocks.append(["spacer", None, None])
            continue
        node_obj = splice_graph.get_node_obj_via_id(node_id)
        node_lend, node_rend = node_obj.get_coords()
        if node_id.startswith("TSS:"):
            site_marks.append((node_lend, "T"))
        elif node_id.startswith("POLYA:"):
            site_marks.append((node_lend, "A"))
        elif node_id.startswith("I:"):
            coord_blocks.append(["intron", node_lend, node_rend])
        else:
            coord_blocks.append(["exon", node_lend, node_rend])

    # a SPACER spans from the end of the preceding block to the start of the next
    for i, block in enumerate(coord_blocks):
        if block[0] != "spacer":
            continue
        prev_block = coord_blocks[i - 1] if i > 0 else None
        next_block = coord_blocks[i + 1] if i + 1 < len(coord_blocks) else None
        if prev_block is None or next_block is None:
            continue  # terminal spacer carries no coordinates; nothing to draw
        block[1] = prev_block[2] + 1
        block[2] = next_block[1] - 1

    exons = []
    for kind, block_lend, block_rend in coord_blocks:
        if kind != "exon":
            continue
        if exons and block_lend <= exons[-1][1] + 1:
            exons[-1][1] = max(exons[-1][1], block_rend)
        else:
            exons.append([block_lend, block_rend])

    spacers = [
        (b[1], b[2]) for b in coord_blocks if b[0] == "spacer" and b[1] is not None and b[2] >= b[1]
    ]
    return [tuple(e) for e in exons], spacers, site_marks


def track_from_simple_path(
    splice_graph, simple_path, name, strand=".", glyph="=", annot="", color=None
):
    """Build a Track from an LRAA simple path; SPACER gaps are drawn as '?'."""
    exons, spacers, site_marks = exons_from_simple_path(splice_graph, simple_path)
    if not exons:
        raise ValueError("simple path {} yields no exon blocks".format(name))
    junc_chars = {}
    for intron in intron_chain(exons):
        for spacer_lend, spacer_rend in spacers:
            if spacer_lend <= intron[1] and intron[0] <= spacer_rend:
                junc_chars[intron] = SPACER_CHAR
    return Track(
        name,
        exons,
        strand,
        glyph,
        annot=annot,
        color=color,
        junc_chars=junc_chars,
        site_marks=tuple(site_marks),
    )


############################
# high-level convenience API
############################


def illustrate(
    entries,
    reference=None,
    contig_acc=None,
    region=None,
    width=100,
    mode="compressed",
    color=False,
    max_ticks=8,
    pad_frac=0.02,
    title=None,
    highlight_char="*",
):
    """Render a stack of structures, annotated against a reference structure.

    `entries` is a list of LRAA Transcript objects, Track objects, or
    (name, exons, strand) / (name, exons, strand, annot) tuples.
    `reference` is the entry (or its name) to compare everything against; when
    given, each other row is annotated with compare_structures().summary() and
    its non-reference introns are drawn with `highlight_char`.

    Entries carrying a contig (Transcript objects) must all share one: a single
    coordinate axis cannot honestly hold two loci, and overlaying them would draw
    structures that align on the page while being unrelated in the genome.

    Returns the rendered text.
    """
    # the reference counts too: when it is a Transcript not already among the
    # entries it gets inserted as a row below, so leaving it out of this check
    # reopens exactly the false overlay the check exists to prevent
    contigs = {
        e.get_contig_acc()
        for e in list(entries) + ([reference] if reference is not None else [])
        if isinstance(e, Transcript)
    }
    if len(contigs) > 1:
        raise ValueError(
            "cannot draw transcripts from more than one contig on one axis: {}".format(
                ", ".join(sorted(contigs))
            )
        )
    if contig_acc is None and len(contigs) == 1:
        contig_acc = next(iter(contigs))

    tracks = [_as_track(e) for e in entries]
    if not tracks:
        raise ValueError("nothing to illustrate")

    ref_track = None
    if reference is not None:
        if isinstance(reference, str):
            matches = [t for t in tracks if t.name == reference]
            if not matches:
                matches = [t for t in tracks if reference in t.name]
            if not matches:
                raise ValueError("reference {} not among the entries".format(reference))
            ref_track = matches[0]
        else:
            ref_track = _as_track(reference)
            if ref_track.name not in [t.name for t in tracks]:
                tracks.insert(0, ref_track)

    if ref_track is not None:
        for track in tracks:
            if track is ref_track:
                if not track.annot:
                    track.annot = "(reference)"
                continue
            cmp_result = compare_structures(track.exons, ref_track.exons)
            track.annot = (track.annot + " " if track.annot else "") + cmp_result.summary()
            for intron in cmp_result.query_only_introns:
                track.junc_chars[intron] = highlight_char

    if region is not None:
        view_lend, view_rend = region
    else:
        view_lend = min(t.lend for t in tracks)
        view_rend = max(t.rend for t in tracks)
        pad = max(1, int((view_rend - view_lend + 1) * pad_frac))
        view_lend, view_rend = view_lend - pad, view_rend + pad

    if contig_acc is None:
        contig_acc = "."
    view = View(
        contig_acc,
        view_lend,
        view_rend,
        width=width,
        mode=mode,
        color=color,
        max_ticks=max_ticks,
    )
    for track in tracks:
        view.add(track)
    return view.render(title)


def _as_track(entry):
    if isinstance(entry, Track):
        return entry
    if isinstance(entry, Transcript):
        return Track(
            entry.get_transcript_id(), exons_of(entry), entry.get_strand()
        )
    name, exons, strand = entry[0], entry[1], entry[2]
    annot = entry[3] if len(entry) > 3 else ""
    return Track(name, list(exons), strand, annot=annot)
