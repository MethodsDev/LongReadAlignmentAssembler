#!/usr/bin/env python3
# encoding: utf-8

"""Collapse an LRAA GTF by splice pattern, and report the TSS/PolyA sites it describes.

Two products, from one parse of the input GTF:

  1. a collapsed GTF, where isoforms sharing BOTH a gene_id and an intron chain become
     one model spanning their outermost termini;
  2. a TSS BED and a PolyA BED, one row per distinct boundary site, carrying the support
     and (for PolyA) the polyadenylation-signal annotation of that site.

The BEDs are derived from the INPUT transcripts, not the collapsed ones. A merged model
represents several terminal variants, so no single boundary count describes it; reading
the uncollapsed side keeps every site's evidence intact and means the collapse never has
to invent one. The same evidence is attached to merged models as index-aligned per-site
lists, so the collapsed GTF stays self-contained.

Support semantics, stated once here because the column name cannot carry it: the value is
the GTF's TSS_read_count / PolyA_read_count. That is the splice graph's clustered sum of
per-read XW normalization weights, ALREADY TRUNCATED TO AN INTEGER by
Transcript.set_TSS_read_count, and for boundaries derived from an annotation rather than
from reads it is a synthetic seed (config['min_alignments_define_TSS_site']) rather than
anything observed. It equals a literal read count only for a read-derived site in a BAM
carrying no XW tag. Hence 'reported_boundary_support' and not 'read_count'.
"""

import os
import re
from collections import defaultdict

from Transcript import Transcript, GTF_contig_to_transcripts
import Util_funcs

# `,` separates sites within a list attribute, so a value that itself holds several
# alternatives must use something else or the lists stop being index-aligned. Only
# reachable for imported or externally merged GTFs, where one coordinate can arrive
# carrying two annotations.
SITE_SEP = ","
WITHIN_SITE_SEP = "|"

# Every per-site list is the same length; a site with no value for one field gets this
# rather than being skipped, which would shift every later index.
MISSING = "."

TSS_BED_COLUMNS = [
    "chrom",
    "start",
    "end",
    "name",
    "score",
    "strand",
    "reported_boundary_support",
    "num_transcripts",
    "transcript_ids",
]

POLYA_BED_COLUMNS = TSS_BED_COLUMNS + ["pas", "pas_offset", "internal_priming"]

_SUPPORT_SEMANTICS_COMMENT = (
    "# reported_boundary_support: the GTF's TSS_read_count/PolyA_read_count -- the "
    "splice graph's clustered sum of per-read XW normalization weights, truncated to "
    "an integer, and a synthetic seed rather than observed reads for "
    "annotation-derived boundaries. Not a literal read count unless the BAM carried "
    "no XW tag."
)


class GeneConflictError(Exception):
    """One intron chain is carried by more than one gene_id.

    Raised instead of collapsing, because the collapsed transcript_id is minted from the
    splice pattern: two such gene_ids would emit duplicate transcript_ids, silently, as
    the GTF format has no uniqueness constraint.
    """

    def __init__(self, conflicts, report_path):
        self.conflicts = conflicts
        self.report_path = report_path
        super().__init__(
            "{} intron pattern(s) are carried by more than one gene_id; see {}. "
            "Isoforms with identical splice patterns are one gene by definition, and "
            "the collapsed transcript_id is minted from the splice pattern, so "
            "collapsing these would emit DUPLICATE transcript_ids -- silently, as the "
            "gtf has no uniqueness check. Refusing rather than picking a gene_id, "
            "because that choice belongs to gene clustering and would desynchronize "
            "this gtf from the expression matrices, which take gene_id from the quant "
            "files. This normally indicates stale or externally produced input, or an "
            "upstream clustering invariant that did not hold -- note that "
            "merge_LRAA_GTFs.py logs a warning and proceeds with unreclustered "
            "transcripts if reclustering raises, which bypasses the "
            "guarantee.".format(len(conflicts), report_path)
        )


class SiteSupportConflictError(Exception):
    """Two transcripts report different support for the same boundary site.

    Every transcript ending at a site copies that support from the same splice-graph
    node, so a disagreement is a bug upstream rather than a case to reconcile here.
    Summing would inflate the site by however many isoforms happen to share it.
    """


class BoundarySite:
    """One distinct (contig, position, strand, type) boundary and its evidence."""

    __slots__ = (
        "contig",
        "position",
        "strand",
        "site_type",
        "transcript_ids",
        "_support",
        "_pas",
        "_pas_offset",
        "_internal_priming",
    )

    def __init__(self, contig, position, strand, site_type):
        self.contig = contig
        self.position = position
        self.strand = strand
        self.site_type = site_type
        self.transcript_ids = []
        self._support = None
        self._pas = set()
        self._pas_offset = set()
        self._internal_priming = set()

    def add(self, transcript_obj, support, pas, pas_offset, internal_priming):
        transcript_id = transcript_obj.get_transcript_id()
        self.transcript_ids.append(transcript_id)

        if support is not None:
            if self._support is not None and self._support != support:
                raise SiteSupportConflictError(
                    "{} site {}:{}({}) is reported with support {} and {} by "
                    "different transcripts ({}). Every transcript at a site copies the "
                    "support of one splice-graph node, so these cannot both be "
                    "right.".format(
                        self.site_type,
                        self.contig,
                        self.position,
                        self.strand,
                        self._support,
                        support,
                        ",".join(sorted(self.transcript_ids)),
                    )
                )
            self._support = support

        # PAS annotation describes the 3' end, so it is meaningful for PolyA sites only.
        # A TSS shared by transcripts with different 3' ends would otherwise look like a
        # site whose motif disagrees with itself.
        if self.site_type == "PolyA":
            if pas is not None:
                self._pas.add(str(pas))
            if pas_offset is not None:
                self._pas_offset.add(str(pas_offset))
            if internal_priming is not None:
                self._internal_priming.add(str(internal_priming))

    @property
    def support(self):
        return MISSING if self._support is None else str(self._support)

    @property
    def pas(self):
        return self._joined(self._pas)

    @property
    def pas_offset(self):
        return self._joined(self._pas_offset)

    @property
    def internal_priming(self):
        return self._joined(self._internal_priming)

    @staticmethod
    def _joined(values):
        if not values:
            return MISSING
        return WITHIN_SITE_SEP.join(sorted(values))

    def bed_row(self):
        name = "{}:{}:{}:{}".format(
            self.site_type, self.contig, self.position, self.strand
        )
        # BED is half open and 0 based; the site is a single genomic base.
        fields = [
            self.contig,
            str(self.position - 1),
            str(self.position),
            name,
            self._bed_score(),
            self.strand,
            self.support,
            str(len(self.transcript_ids)),
            ",".join(sorted(self.transcript_ids)),
        ]
        if self.site_type == "PolyA":
            fields += [self.pas, self.pas_offset, self.internal_priming]
        return "\t".join(fields)

    def _bed_score(self):
        """Column 5 exists for genome browsers, which clamp to 0-1000.

        The unclamped value is in reported_boundary_support; this column is display.
        """
        if self._support is None:
            return "0"
        try:
            return str(min(1000, int(round(float(self._support)))))
        except (TypeError, ValueError):
            return "0"


def boundary_sites(contig_to_transcripts):
    """Distinct TSS and PolyA sites of the supplied transcripts, keyed by site.

    Returns {(contig, position, strand, site_type): BoundarySite}. A transcript
    contributes a site only where it claims one: TSS "True" / PolyA "True".
    """
    sites = {}

    for contig in sorted(contig_to_transcripts.keys()):
        for transcript_obj in contig_to_transcripts[contig]:
            strand = transcript_obj.get_strand()
            lend, rend = transcript_obj.get_coords()

            if transcript_obj.has_TSS():
                position = lend if strand == "+" else rend
                _site(sites, contig, position, strand, "TSS").add(
                    transcript_obj,
                    transcript_obj.get_TSS_read_count(),
                    None,
                    None,
                    None,
                )

            if transcript_obj.has_PolyA():
                position = rend if strand == "+" else lend
                _site(sites, contig, position, strand, "PolyA").add(
                    transcript_obj,
                    transcript_obj.get_PolyA_read_count(),
                    transcript_obj.get_polyA_signal(),
                    transcript_obj.get_polyA_signal_offset(),
                    transcript_obj.get_likely_internal_primed(),
                )

    return sites


def _site(sites, contig, position, strand, site_type):
    key = (contig, position, strand, site_type)
    if key not in sites:
        sites[key] = BoundarySite(contig, position, strand, site_type)
    return sites[key]


def write_site_bed(path, sites, site_type, provenance_comments=()):
    """One row per site of the requested type, ordered by genomic position."""
    columns = TSS_BED_COLUMNS if site_type == "TSS" else POLYA_BED_COLUMNS
    selected = [site for site in sites.values() if site.site_type == site_type]
    selected.sort(key=lambda s: (s.contig, s.position, s.strand))

    with open(path, "wt") as ofh:
        for comment in provenance_comments:
            ofh.write(comment.rstrip("\n") + "\n")
        ofh.write(_SUPPORT_SEMANTICS_COMMENT + "\n")
        ofh.write("#" + "\t".join(columns) + "\n")
        for site in selected:
            ofh.write(site.bed_row() + "\n")

    return len(selected)


def collapse_gtf(
    input_gtf,
    output_gtf,
    merge_report=None,
    gene_conflicts_report=None,
    TSS_bed=None,
    PolyA_bed=None,
    nonfatal_gene_conflicts=False,
):
    """Collapse `input_gtf` by splice pattern and write the site BEDs beside it.

    Returns {output name: path} for what was actually written.

    On a cross-gene splice-pattern conflict the collapse cannot proceed. The BEDs can:
    they describe the input, not the collapse. So `nonfatal_gene_conflicts` decides
    between two whole contracts, and callers must pick one rather than blend them:

      False (the standalone CLI): write only the conflicts report and raise, leaving no
        partial output for a later reader to mistake for a collapse that worked.
      True (LRAA and the WDL task): write the conflicts report and both BEDs, skip the
        collapsed GTF and merge report, and return. A failed task publishes nothing, so
        raising here would destroy a completed run's primary outputs along with the very
        report needed to diagnose the conflict.
    """
    merge_report = merge_report or (output_gtf + ".isoform_merge_report.tsv")
    gene_conflicts_report = gene_conflicts_report or (output_gtf + ".gene_conflicts.tsv")
    TSS_bed = TSS_bed or (output_gtf + ".TSS.bed")
    PolyA_bed = PolyA_bed or (output_gtf + ".PolyA.bed")

    contig_to_input_transcripts = GTF_contig_to_transcripts.parse_GTF_to_Transcripts(
        input_gtf
    )

    provenance = _provenance_comments(input_gtf)
    sites = boundary_sites(contig_to_input_transcripts)

    # Checked before any output is opened, so the fatal policy really is all or nothing.
    #
    # Not reconciled here, deliberately. Unifying the gene_ids would be a gene
    # REASSIGNMENT, and this is not where gene identity is decided: the expression
    # matrices key their features on the gene_id carried by the quant files
    # (build_LRAA_expr_matrices.py:132-135), so a gene_id invented here would appear in
    # the collapsed gtf and nowhere else. Identical splice patterns are made one gene
    # upstream instead, by the clustering that assigns them
    # (GeneCommunityCluster._contract_identical_chains).
    conflicts = find_cross_gene_splice_patterns(contig_to_input_transcripts)
    write_gene_conflicts_report(gene_conflicts_report, conflicts)

    written = {"gene_conflicts_report": gene_conflicts_report}

    if conflicts and not nonfatal_gene_conflicts:
        raise GeneConflictError(conflicts, gene_conflicts_report)

    write_site_bed(TSS_bed, sites, "TSS", provenance)
    write_site_bed(PolyA_bed, sites, "PolyA", provenance)
    written["TSS_bed"] = TSS_bed
    written["PolyA_bed"] = PolyA_bed

    if conflicts:
        return written

    _write_collapsed_gtf(
        output_gtf,
        merge_report,
        contig_to_input_transcripts,
        sites,
        provenance,
    )
    written["collapsed_gtf"] = output_gtf
    written["merge_report"] = merge_report

    # Cheap, and it pins the property the conflict refusal exists for: the collapsed id
    # is derived from the splice pattern, so a duplicate here means two gene_ids shared
    # a pattern and the check missed it.
    assert_unique_transcript_ids(output_gtf)

    return written


def _write_collapsed_gtf(
    output_gtf, merge_report, contig_to_input_transcripts, sites, provenance
):
    with open(output_gtf, "wt") as ofh, open(merge_report, "wt") as merge_ofh:
        for comment in provenance:
            ofh.write(comment.rstrip("\n") + "\n")

        merge_ofh.write(
            "\t".join(
                [
                    "collapsed_transcript_id",
                    "gene_id",
                    "num_isoforms_merged",
                    "merged_transcript_ids",
                ]
            )
            + "\n"
        )

        for contig, transcript_obj_list in contig_to_input_transcripts.items():

            transcripts_to_output = list()

            # Grouped by (gene_id, splice_pattern_code): isoforms merge only when they
            # share the same gene_id AND the same intron chain. Cross-gene patterns were
            # refused above, so every pattern here belongs to exactly one gene_id.
            gene_splice_to_transcripts = defaultdict(list)

            for transcript_obj in transcript_obj_list:
                if transcript_obj.has_introns():
                    gene_splice_to_transcripts[
                        (
                            transcript_obj.get_gene_id(),
                            splice_pattern_code(transcript_obj),
                        )
                    ].append(transcript_obj)
                else:
                    transcripts_to_output.append(transcript_obj)

            for (
                gene_id,
                splice_pattern,
            ), transcripts_same_group_list in gene_splice_to_transcripts.items():

                if "^" in gene_id:
                    new_transcript_id = "^".join(
                        [gene_id.split("^")[0], splice_pattern]
                    )
                else:
                    new_transcript_id = splice_pattern

                member_ids = sorted(
                    [t.get_transcript_id() for t in transcripts_same_group_list]
                )

                if len(transcripts_same_group_list) == 1:
                    # Left as it is, scalar boundary attributes and all: a single-member
                    # group has exactly one TSS and one PolyA, so there is nothing to
                    # attribute and the per-site lists would only restate it.
                    transcript_obj = transcripts_same_group_list[0]
                    transcript_obj.set_transcript_id(new_transcript_id)
                    transcripts_to_output.append(transcript_obj)
                else:
                    merged_isoform = merge_isoforms(transcripts_same_group_list, sites)
                    merged_isoform.set_transcript_id(new_transcript_id)
                    transcripts_to_output.append(merged_isoform)

                merge_ofh.write(
                    "\t".join(
                        [
                            new_transcript_id,
                            gene_id,
                            str(len(member_ids)),
                            ",".join(member_ids),
                        ]
                    )
                    + "\n"
                )

            transcripts_to_output = sorted(
                transcripts_to_output, key=lambda x: x._exon_segments[0][0]
            )

            for transcript_obj in transcripts_to_output:
                ofh.write(transcript_obj.to_GTF_format(include_TPM=False) + "\n")


def _provenance_comments(input_gtf):
    """The input's own leading comments, plus a line naming this step.

    Reserializing transcripts drops the header the input carried, which is how a
    collapsed GTF ends up with no record of the LRAA version that produced its models.
    """
    comments = []
    with open(input_gtf, "rt") as fh:
        for line in fh:
            if not line.startswith("#"):
                break
            comments.append(line.rstrip("\n"))
    comments.append(
        "# splice-pattern collapsed from {}".format(os.path.basename(input_gtf))
    )
    return comments


def splice_pattern_code(transcript_obj):
    """The hash the collapsed transcript_id is minted from.

    Its input includes contig and strand (Transcript.get_introns_string), so the code
    identifies a splice pattern genome-wide rather than within a contig.
    """
    return Util_funcs.get_hash_code(transcript_obj.get_introns_string())


def find_cross_gene_splice_patterns(contig_to_input_transcripts):
    """[(contig, splice_pattern_code, {gene_id: [transcript_id, ...]})] for patterns
    carried by more than one gene_id. Empty list when the input is well formed."""
    found = []
    for contig, transcript_obj_list in contig_to_input_transcripts.items():
        pattern_to_genes = defaultdict(lambda: defaultdict(list))
        for transcript_obj in transcript_obj_list:
            if not transcript_obj.has_introns():
                continue
            pattern_to_genes[splice_pattern_code(transcript_obj)][
                transcript_obj.get_gene_id()
            ].append(transcript_obj.get_transcript_id())
        for code in sorted(pattern_to_genes):
            genes = pattern_to_genes[code]
            if len(genes) > 1:
                found.append((contig, code, genes))
    return found


def write_gene_conflicts_report(path, conflicts):
    """Written on every run, conflicts or not: it is a declared output, and a missing
    file is indistinguishable from a run that never checked."""
    with open(path, "wt") as fh:
        fh.write(
            "\t".join(
                [
                    "contig",
                    "splice_pattern_code",
                    "num_gene_ids",
                    "gene_id",
                    "num_isoforms",
                    "transcript_ids",
                ]
            )
            + "\n"
        )
        for contig, code, genes in conflicts:
            for gene_id in sorted(genes):
                transcript_ids = sorted(genes[gene_id])
                fh.write(
                    "\t".join(
                        [
                            contig,
                            code,
                            str(len(genes)),
                            gene_id,
                            str(len(transcript_ids)),
                            ",".join(transcript_ids),
                        ]
                    )
                    + "\n"
                )


def assert_unique_transcript_ids(gtf_file):
    seen = set()
    with open(gtf_file, "rt") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            cols = line.rstrip("\n").split("\t")
            if len(cols) < 9 or cols[2] != "transcript":
                continue
            match = re.search(r'transcript_id\s+"([^"]*)"', cols[8])
            if match is None:
                continue
            tid = match.group(1)
            if tid in seen:
                raise AssertionError(
                    "collapsed gtf {} contains duplicate transcript_id {}. This is a "
                    "bug in the collapse, not bad input: cross-gene splice patterns "
                    "are refused before writing.".format(gtf_file, tid)
                )
            seen.add(tid)


def merge_isoforms(transcript_obj_list, sites=None):
    """One model spanning the outermost termini of isoforms sharing a splice pattern.

    When `sites` is supplied, the merged model also carries per-site lists: the
    coordinates it absorbed, and the support and PAS annotation of each. Every list is
    built from one pass over the same sorted coordinates, so index i of every list
    describes the same site -- without that, a model with two PolyA sites at different
    support says which coordinates exist but not which count belongs to which.
    """

    first_transcript_obj = transcript_obj_list[0]
    template_exon_coords = first_transcript_obj.get_exon_segments()
    min_lend = template_exon_coords[0][0]
    max_rend = template_exon_coords[-1][1]
    contig_acc = first_transcript_obj.get_contig_acc()
    contig_strand = first_transcript_obj.get_strand()
    gene_id = first_transcript_obj.get_gene_id()
    has_TSS = first_transcript_obj.has_TSS()
    has_PolyA = first_transcript_obj.has_PolyA()

    merged_isoform_ids = [first_transcript_obj.get_transcript_id()]

    candidate_TSS_sites = set()
    candidate_PolyA_sites = set()

    if has_TSS:
        candidate_TSS_sites.add(min_lend if contig_strand == "+" else max_rend)

    if has_PolyA:
        candidate_PolyA_sites.add(max_rend if contig_strand == "+" else min_lend)

    for transcript_obj in transcript_obj_list[1:]:
        exon_coords = transcript_obj.get_exon_segments()
        lend = exon_coords[0][0]
        rend = exon_coords[-1][1]

        if lend < min_lend:
            min_lend = lend

            if contig_strand == "+":
                has_TSS = transcript_obj.has_TSS()
            else:
                has_PolyA = transcript_obj.has_PolyA()

        if rend > max_rend:
            max_rend = rend
            if contig_strand == "+":
                has_PolyA = transcript_obj.has_PolyA()
            else:
                has_TSS = transcript_obj.has_TSS()

        if transcript_obj.has_TSS():
            candidate_TSS_sites.add(lend if contig_strand == "+" else rend)

        if transcript_obj.has_PolyA():
            candidate_PolyA_sites.add(rend if contig_strand == "+" else lend)

        merged_isoform_ids.append(transcript_obj.get_transcript_id())

    template_exon_coords[0][0] = min_lend
    template_exon_coords[-1][1] = max_rend

    merged_transcript = Transcript(contig_acc, template_exon_coords, contig_strand)
    merged_transcript.set_gene_id(gene_id)

    merged_transcript._imported_has_TSS = has_TSS
    merged_transcript._imported_has_POLYA = has_PolyA

    merged_transcript.add_meta(
        "merged_isoforms_shared_splice_pattern", ",".join(sorted(merged_isoform_ids))
    )

    _add_site_meta(
        merged_transcript,
        contig_acc,
        contig_strand,
        "TSS",
        sorted(candidate_TSS_sites),
        sites,
    )
    _add_site_meta(
        merged_transcript,
        contig_acc,
        contig_strand,
        "PolyA",
        sorted(candidate_PolyA_sites),
        sites,
    )

    return merged_transcript


def _add_site_meta(
    merged_transcript, contig_acc, contig_strand, site_type, positions, sites
):
    """Attach `<type>_sites` and its parallel evidence lists, in one shared order.

    Positions are genomic ascending, which is what the coordinate list has always been;
    on the minus strand that means index 0 is the 3'-most TSS. The evidence lists are
    generated from this same sequence rather than sorted independently, so they cannot
    drift out of alignment.
    """
    if not positions:
        return

    merged_transcript.add_meta(
        "{}_sites".format(site_type), SITE_SEP.join(str(p) for p in positions)
    )

    if sites is None:
        return

    resolved = [
        sites.get((contig_acc, position, contig_strand, site_type))
        for position in positions
    ]

    merged_transcript.add_meta(
        "{}_site_support".format(site_type),
        SITE_SEP.join(MISSING if s is None else s.support for s in resolved),
    )

    if site_type != "PolyA":
        return

    for attribute, accessor in (
        ("PolyA_site_PAS", lambda s: s.pas),
        ("PolyA_site_PAS_offset", lambda s: s.pas_offset),
        ("PolyA_site_internal_priming", lambda s: s.internal_priming),
    ):
        merged_transcript.add_meta(
            attribute,
            SITE_SEP.join(MISSING if s is None else accessor(s) for s in resolved),
        )
