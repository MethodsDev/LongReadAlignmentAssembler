#!/usr/bin/env python3

import math
import sys, os, re
import Transcript
import MultiPath
import MultiPathCounter
import Simple_path_utils as SPU
from Quantify import Quantify
from collections import defaultdict
import LRAA_Globals
from LRAA_Globals import SPACER, DEBUG
import logging
from math import log
import intervaltree as itree

logger = logging.getLogger(__name__)


def filter_transcripts_by_min_length(transcripts, min_transcript_length):
    """Retain only transcripts meeting minimum cDNA length.

    A model containing an expressed reference transcript is exempt: the annotation
    asserts a transcript there, and its component was admitted for assembly on that
    basis, so dropping it here would discard its reads without any output row."""

    if min_transcript_length is None or min_transcript_length <= 0:
        return transcripts

    transcripts_retained = []
    for transcript in transcripts:
        if (
            transcript.contains_reference_model()
            and LRAA_Globals.config["ref_trans_filter_mode"] == "retain_expressed"
            and transcript.get_TPM() > 0
        ):
            transcripts_retained.append(transcript)
            continue

        if transcript.get_cdna_len() >= min_transcript_length:
            transcripts_retained.append(transcript)

    return transcripts_retained


def _cell_barcode_from_read_name(read_name):
    """Cell barcode encoded in a read name, or None for bulk input.

    Util_funcs.get_read_name_include_sc_encoding() emits
    "<cell_barcode>^<umi>^<read_name>" only when both the barcode and UMI tags are
    present on the alignment, and the bare read name otherwise. The separator is
    therefore what distinguishes single-cell input from bulk, and any filter keyed
    on it self-disables on bulk data.
    """

    separator_pos = read_name.find("^")
    if separator_pos <= 0:
        return None
    return read_name[:separator_pos]


def _input_has_cell_barcodes(transcripts, probe_limit=10000):
    """Whether these transcripts carry cell-encoded read names.

    Single-cell input answers on the first read examined; bulk input has to be
    probed, so the scan is bounded. Determined up front rather than during
    filtering, so the verdict cannot depend on the order transcripts are visited.
    """

    num_examined = 0
    for transcript in transcripts:
        for mp in transcript.get_multipaths_evidence_assigned():
            for read_name in mp.get_read_names():
                if _cell_barcode_from_read_name(read_name) is not None:
                    return True
                num_examined += 1
                if num_examined >= probe_limit:
                    return False
    return False


def filter_novel_monoexonic_isoforms_by_min_cells(
    transcripts, min_supporting_cells, cell_roster=None
):
    """Require a novel monoexonic model to be seen in a minimum number of cells.

    Single-cell only, and monoexonic only. A monoexonic model has no intron chain
    corroborating it, so the number of distinct cells it was observed in is the one
    independent axis of evidence available: reads all drawn from one cell describe a
    single amplification event, not a transcript the population expresses.

    A model containing a reference model is exempt, on the same basis as the
    reference reprieve in filter_monoexonic_isoforms_by_TPM_threshold(): the
    annotation already asserts a transcript there, so prevalence is not the evidence
    being asked for. Only novel monoexonic models must clear the bar.

    An absolute cell count rather than a fraction of the cells present. A fraction
    makes the bar scale with however many cells a cluster happens to contain, so the
    same model passes in a small cluster and fails in a large one; against 14 PBMC
    clusters spanning 122 to 1,506 cells, recovery of reference-matching monoexons
    under a fixed bar was 98% at 3 cells, 92% at 5 and 73% at 10, while a
    5%-of-roster rule recovered well under half of them.

    cell_roster, when given, is the set of barcodes accepted as real cells. Barcodes
    outside it are empty droplets or ambient signal; counting them would let ambient
    RNA, which is spread thinly across many droplets, masquerade as prevalence.
    """

    if min_supporting_cells is None or min_supporting_cells <= 0:
        return transcripts

    retain_reference = (
        LRAA_Globals.config["ref_trans_filter_mode"] == "retain_expressed"
    )

    if not _input_has_cell_barcodes(transcripts):
        # Bulk input encodes no barcode in read names; there is no cell axis to judge on.
        return transcripts

    transcripts_retained = list()
    num_filtered = 0

    for transcript in transcripts:
        if not transcript.is_monoexonic():
            transcripts_retained.append(transcript)
            continue

        if retain_reference and transcript.contains_reference_model():
            transcripts_retained.append(transcript)
            continue

        cells = set()
        for mp in transcript.get_multipaths_evidence_assigned():
            for read_name in mp.get_read_names():
                cell_barcode = _cell_barcode_from_read_name(read_name)
                if cell_barcode is None:
                    continue
                if cell_roster and cell_barcode not in cell_roster:
                    continue
                cells.add(cell_barcode)

        if len(cells) >= min_supporting_cells:
            transcripts_retained.append(transcript)
            continue

        num_filtered += 1
        logger.debug(
            "FILTERING novel monoexonic transcript %s: supported by %d cells, "
            "below min_monoexonic_supporting_cells=%d",
            transcript.get_transcript_id(),
            len(cells),
            min_supporting_cells,
        )

    if num_filtered:
        logger.info(
            "-filtered %d novel monoexonic transcripts supported by fewer than %d cells",
            num_filtered,
            min_supporting_cells,
        )

    return transcripts_retained


def _peak_overlapping_read_fraction(spans):
    """Largest fraction of these read spans that mutually overlap at one base.

    Equivalent to peak read depth divided by read count. A sweep is used rather than
    pairwise comparison so this stays linear in read count. Returns None when there is
    nothing to judge.
    """

    if not spans:
        return None

    events = []
    for lend, rend in spans:
        events.append((lend, 1))
        events.append((rend + 1, -1))
    events.sort()

    depth = 0
    peak = 0
    for _, delta in events:
        depth += delta
        if depth > peak:
            peak = depth

    return peak / len(spans)


def _collect_supporting_read_spans(transcript, read_name_to_span):
    """Genomic spans of the real reads assigned to this transcript."""

    read_names = set()
    for mp in transcript.get_multipaths_evidence_assigned():
        read_names.update(mp.get_read_names())

    spans = []
    for read_name in read_names:
        if read_name.startswith("reftranscript:") or read_name.startswith(
            "fake_for_merge"
        ):
            continue
        span = read_name_to_span.get(read_name)
        if span is not None:
            spans.append(span)

    return spans


def filter_monoexonic_isoforms_by_adjusted_TPM_ratio(
    transcripts, read_name_to_span, min_ratio
):
    """Drop monoexonic models whose read-count TPM and coverage-depth TPM disagree.

    The reported TPM counts reads, so it treats a read as evidence for the whole
    model no matter how little of it the read covers. The adjusted TPM instead
    divides the read bases landing on the model by the model's length, giving mean
    coverage depth in the same units:

        TPM     = read_count                    / num_total_reads * 1e6
        adjTPM  = (sum aligned bases / L_model) / num_total_reads * 1e6

    When the reads agree -- each spanning the model -- the two are similar. When they
    tile a long span, each contributing a fraction, adjTPM falls far below TPM. Their
    ratio is the mean fraction of the model an individual read covers, which is
    independent of library size and of read count, so it isolates agreement rather
    than abundance. Note that an absolute threshold on adjTPM does NOT work here: a
    20kb model tiled by 188 short reads still averages ~5x depth and would pass, while
    genuinely coherent 3-read models would fail.

    Spans are alignment extents, so a spliced read overlapping the model is credited
    its full overlap including any intron. That inflates the ratio, making the filter
    conservative, which is the safe direction.
    """

    if min_ratio is None or min_ratio <= 0:
        return transcripts

    transcripts_retained = list()
    num_filtered = 0

    for transcript in transcripts:
        if not transcript.is_monoexonic():
            transcripts_retained.append(transcript)
            continue

        spans = _collect_supporting_read_spans(transcript, read_name_to_span)
        model_lend, model_rend = transcript.get_coords()
        model_length = model_rend - model_lend + 1

        if not spans or model_length <= 0:
            transcripts_retained.append(transcript)
            continue

        covered_bases = 0
        for span_lend, span_rend in spans:
            overlap = min(model_rend, span_rend) - max(model_lend, span_lend) + 1
            if overlap > 0:
                covered_bases += overlap

        # ratio == (covered_bases / model_length) / read_count
        ratio = (covered_bases / model_length) / len(spans)

        if ratio >= min_ratio:
            transcripts_retained.append(transcript)
            continue

        num_filtered += 1
        logger.debug(
            "FILTERING monoexonic transcript %s: adjusted/plain TPM ratio %.3f over "
            "%d reads and %d bp (min_monoexonic_adjusted_TPM_ratio=%.2f)",
            transcript.get_transcript_id(),
            ratio,
            len(spans),
            model_length,
            min_ratio,
        )

    if num_filtered:
        logger.info(
            "-filtered %d monoexonic transcripts whose coverage-depth TPM diverges "
            "from their read-count TPM",
            num_filtered,
        )

    return transcripts_retained


def filter_monoexonic_isoforms_by_read_span_coherence(
    transcripts, read_name_to_span, min_peak_frac
):
    """Drop monoexonic models whose supporting reads do not stack on a common interval.

    Multi-exonic models are untouched: their intron chain is the corroborating
    evidence. A monoexonic model has none, so we require its reads to agree on where
    the transcript is. Reads that tile a span without overlapping describe contiguous
    coverage -- unspliced pre-mRNA, intronic signal, a gene body -- rather than one
    molecule, and collapsing that into a single model invents a transcript nothing
    observed end to end.
    """

    if min_peak_frac is None or min_peak_frac <= 0:
        return transcripts

    transcripts_retained = list()
    num_filtered = 0

    for transcript in transcripts:
        if not transcript.is_monoexonic():
            transcripts_retained.append(transcript)
            continue

        spans = _collect_supporting_read_spans(transcript, read_name_to_span)

        peak_frac = _peak_overlapping_read_fraction(spans)

        # No recorded spans means no basis to reject; leave the call to the other
        # filters rather than dropping a model over missing bookkeeping.
        if peak_frac is None or peak_frac >= min_peak_frac:
            transcripts_retained.append(transcript)
            continue

        num_filtered += 1
        logger.debug(
            "FILTERING monoexonic transcript %s: only %.2f of its %d supporting reads "
            "overlap at any single base (min_monoexonic_read_span_peak_frac=%.2f)",
            transcript.get_transcript_id(),
            peak_frac,
            len(spans),
            min_peak_frac,
        )

    if num_filtered:
        logger.info(
            "-filtered %d monoexonic transcripts lacking mutually overlapping read support",
            num_filtered,
        )

    return transcripts_retained


def filter_monoexonic_isoforms_by_TPM_threshold(transcripts, min_TPM):

    transcripts_retained = list()
    hifi_mode = LRAA_Globals.config.get("HiFi", False)

    if hifi_mode:
        logger.info("HiFi mode: single-exon transcripts must have TSS or PolyA annotation to be retained")

    for transcript in transcripts:
        tpm = transcript.get_TPM()

        # reftrans logic:
        if (
            transcript.contains_reference_model()
            and LRAA_Globals.config["ref_trans_filter_mode"] == "retain_expressed"
            and tpm > 0
        ):
            transcripts_retained.append(transcript)
            continue

        # regular filter logic
        if transcript.is_monoexonic():
            # In HiFi mode, require TSS or PolyA annotation for single-exon transcripts
            if hifi_mode and not (transcript.has_TSS() or transcript.has_PolyA()):
                # filter out: single-exon without boundary annotation in HiFi mode
                pass
            elif tpm < min_TPM:
                # filter out: below TPM threshold
                pass
            else:
                # keep: meets all criteria
                transcripts_retained.append(transcript)
        else:
            # keep: multi-exonic
            transcripts_retained.append(transcript)

    return transcripts_retained


def filter_isoforms_by_min_isoform_fraction(
    transcripts, min_isoform_fraction, run_EM, max_EM_iterations
):

    min_frac_gene_unique_reads = LRAA_Globals.config["min_frac_gene_unique_reads"]

    # Build contig/strand prefix from transcripts if available
    try:
        _pref_t = transcripts[0] if transcripts else None
        _ca = _pref_t.get_contig_acc() if _pref_t else None
        _cs = _pref_t.get_strand() if _pref_t else None
        _prefix = f"[{_ca}{_cs}] " if _ca and _cs else ""
    except Exception:
        _prefix = ""

    logger.info(
        f"{_prefix}Filtering transcripts according to min isoform fraction: {min_isoform_fraction}"
    )

    transcript_id_to_transcript_obj = dict(
        [(x.get_transcript_id(), x) for x in transcripts]
    )

    def get_isoform_unique_assigned_read_count(transcript_id, frac_read_assignments):
        num_unique_reads = 0
        for mp in frac_read_assignments[transcript_id]:
            if (
                frac_read_assignments[transcript_id][mp]
                >= LRAA_Globals.config["unique_read_filter_min_frac"]
            ):
                num_unique_reads += mp.get_read_count()

        return num_unique_reads

    gene_id_to_transcripts = _get_gene_id_to_transcripts(transcripts)

    all_transcripts_retained = list()

    # examine each gene separately:

    for gene_id, isoforms_of_gene in gene_id_to_transcripts.items():

        isoforms_were_filtered = True  # init for loop

        q = Quantify(run_EM, max_EM_iterations, quant_mode="draft")

        filtering_round = 0

        frac_read_assignments = None

        while isoforms_were_filtered:

            filtering_round += 1
            num_total_isoforms = 0
            num_filtered_isoforms = 0
            transcripts_retained = list()
            isoforms_were_filtered = (
                False  # update to True if we do filter an isoform out.
            )

            # Build a [contig+strand] prefix for this gene's isoforms
            try:
                _t0 = next(iter(isoforms_of_gene)) if isoforms_of_gene else None
                _ca = _t0.get_contig_acc() if _t0 else None
                _cs = _t0.get_strand() if _t0 else None
                _prefix = f"[{_ca}{_cs}] " if _ca and _cs else None
            except Exception:
                _prefix = None

            frac_read_assignments = q._estimate_isoform_read_support(
                isoforms_of_gene, prefix_str=_prefix
            )

            gene_id_to_read_count = Quantify.get_gene_read_counts(
                frac_read_assignments, transcript_id_to_transcript_obj
            )

            gene_read_count = gene_id_to_read_count[gene_id]

            num_isoforms_of_gene_filtered = 0

            isoforms_of_gene = sorted(
                isoforms_of_gene, key=lambda x: x.get_isoform_fraction()
            )
            num_isoforms_of_gene = len(isoforms_of_gene)

            for transcript in isoforms_of_gene:

                num_total_isoforms += 1
                transcript_id = transcript.get_transcript_id()

                transcript_unique_read_count = get_isoform_unique_assigned_read_count(
                    transcript_id, frac_read_assignments
                )

                frac_gene_unique_reads = (
                    transcript_unique_read_count / gene_read_count
                    if gene_read_count > 0
                    else 0
                )

                logger.debug(
                    "Transcript_id: {} has unique read frac of gene total reads: {}".format(
                        transcript_id, frac_gene_unique_reads
                    )
                )

                ## first check to see if we should retain a reftrans
                if (
                    transcript.contains_reference_model()
                    and LRAA_Globals.config["ref_trans_filter_mode"]
                    == "retain_expressed"
                    and transcript.get_TPM() > 0
                ):
                    transcripts_retained.append(transcript)
                    continue

                ## if tons of isoforms, allow pruning of multiple in a single round
                if (
                    num_isoforms_of_gene - num_isoforms_of_gene_filtered
                    > LRAA_Globals.config[
                        "min_isoform_count_aggressive_filtering_iso_fraction"
                    ]
                    and frac_gene_unique_reads < min_frac_gene_unique_reads
                    and transcript.get_isoform_fraction() < min_isoform_fraction
                ):
                    isoforms_were_filtered = True
                    num_filtered_isoforms += 1
                    num_isoforms_of_gene_filtered += 1

                # standard isoform fraction based filtering
                elif not isoforms_were_filtered and (
                    frac_gene_unique_reads < min_frac_gene_unique_reads
                    or transcript.get_isoform_fraction() < min_isoform_fraction
                ):

                    isoforms_were_filtered = True
                    num_filtered_isoforms += 1
                    num_isoforms_of_gene_filtered += 1

                    logger.debug(
                        "Filtering out transcript_id {} as low fraction of unique reads: {}".format(
                            transcript_id, frac_gene_unique_reads
                        )
                    )

                elif not isoforms_were_filtered and (
                    transcript.has_novel_splice_pattern() is True
                    and transcript_unique_read_count
                    < LRAA_Globals.config["min_unique_reads_novel_isoform"]
                ):
                    isoforms_were_filtered = True
                    num_filtered_isoforms += 1
                    num_isoforms_of_gene_filtered += 1

                    logger.debug(
                        "Filtering out transcript_id {} as novel isoform with too few unique reads: {}".format(
                            transcript_id, transcript_unique_read_count
                        )
                    )

                else:
                    transcripts_retained.append(transcript)

            logger.debug(
                "gene {} based isoform filtering round {} involved filtering of {} isoforms / {} total isoforms".format(
                    gene_id,
                    filtering_round,
                    num_filtered_isoforms,
                    num_total_isoforms,
                )
            )

            # reset list of transcripts
            isoforms_of_gene = transcripts_retained

        # done with individual gene-based isoform filtering
        all_transcripts_retained.extend(isoforms_of_gene)

    return all_transcripts_retained


def _get_gene_id_to_transcripts(transcripts):

    gene_id_to_transcripts = defaultdict(set)
    for transcript in transcripts:
        gene_id = transcript.get_gene_id()
        gene_id_to_transcripts[gene_id].add(transcript)

    return gene_id_to_transcripts


def prune_likely_degradation_products(transcripts, splice_graph, frac_read_assignments):

    try:
        _ca = splice_graph.get_contig_acc()
        _cs = splice_graph.get_contig_strand()
        _prefix = f"[{_ca}{_cs}] " if _ca and _cs else ""
    except Exception:
        _prefix = ""

    logger.info(f"{_prefix}Pruning likely degradation products")

    sg = splice_graph

    transcript_id_to_transcript_obj = dict(
        [(x.get_transcript_id(), x) for x in transcripts]
    )
    gene_read_counts = Quantify.get_gene_read_counts(
        frac_read_assignments, transcript_id_to_transcript_obj
    )

    transcripts_ret = list()  # transcripts not pruned and returned.

    # first organize by gene_id
    gene_id_to_transcripts = _get_gene_id_to_transcripts(transcripts)

    for gene_id, transcript_set in gene_id_to_transcripts.items():

        gene_read_count = gene_read_counts[gene_id]

        if gene_read_count == 0:
            # should explore why this is. ref-trans-only?
            # gbye!
            continue

        if len(transcript_set) == 1:
            transcripts_ret.extend(list(transcript_set))
            continue

        # compare isoforms for compatibility ignoring TSS and PolyA sites
        # if an isoform is fully contained by another and has substantially less expression, prune it.

        transcript_list = list(transcript_set)
        contig_strand = transcript_list[0].get_strand()

        """  #TODO: figure out sorting order to speed this up.
        transcript_list = sorted(transcript_list, key=lambda x: (x.get_coords()[0], x.get_coords()[1], x.get_left_boundary_sort_weight(), x.get_right_boundary_sort_weight()))

        if contig_strand == '-':
            transcript_list = list(sorted(transcript_list, key=lambda x: (-1 * x.get_coords()[1],
                                                                          -1 * x.get_coords()[0],
                                                                          -1 * x.get_right_boundary_sort_weight(),
                                                                          -1 * x.get_left_boundary_sort_weight()
                                                                       ) ) ) 
        """

        # Longest first, so a longer model is the candidate subsumer of a shorter one.
        # Equal cdna lengths are a real tie: terminal variants of one intron chain can
        # each qualify as the other's degradation product, and the structural tests
        # below run before the expression test, so whichever is visited first wins
        # outright. Transcript hashes by identity, so that used to be address order.
        # Break the tie on assigned reads, then on structure, so the better-supported
        # model survives and absorbs the weaker one. Read counts are a stable snapshot:
        # absorbing transfers multipath evidence, not assigned counts.
        transcript_list = sorted(
            transcript_list,
            key=lambda x: (
                -x.get_cdna_len(),
                -x.get_read_counts_assigned(),
                Transcript.Transcript.structural_sort_key(x),
            ),
        )

        transcript_prune_as_degradation = set()
        for i in range(len(transcript_list)):
            transcript_i = transcript_list[i]
            transcript_i_id = transcript_i.get_transcript_id()

            gene_i_id = transcript_i.get_gene_id()
            assert gene_i_id == gene_id, "Error: gene_i_id {} != gene_id {} ".format(
                gene_i_id, gene_id
            )

            if transcript_i in transcript_prune_as_degradation:
                continue

            transcript_i_simple_path = transcript_i.get_simple_path()
            i_path_trimmed, i_TSS_id, i_polyA_id = SPU.trim_TSS_and_PolyA(
                transcript_i_simple_path, contig_strand
            )
            transcript_i_read_counts_assigned = transcript_i.get_read_counts_assigned()

            frac_gene_expression_i = transcript_i_read_counts_assigned / gene_read_count

            # for j in range(i+1, len(transcript_list)):

            for j in range(len(transcript_list)):

                if i == j:
                    continue

                transcript_j = transcript_list[j]

                if transcript_j in transcript_prune_as_degradation:
                    continue

                transcript_j_id = transcript_j.get_transcript_id()
                gene_j_id = transcript_j.get_gene_id()

                assert gene_j_id == gene_id

                transcript_j_simple_path = transcript_j.get_simple_path()
                j_path_trimmed, j_TSS_id, j_polyA_id = SPU.trim_TSS_and_PolyA(
                    transcript_j_simple_path, contig_strand
                )
                transcript_j_read_counts_assigned = (
                    transcript_j.get_read_counts_assigned()
                )

                if len(j_path_trimmed) > len(i_path_trimmed):
                    # no way can i subsume j
                    continue

                frac_gene_expression_j = (
                    transcript_j_read_counts_assigned / gene_read_count
                )

                logger.debug(
                    "Exploring path_i: {} {} as subsuming path_j {} {}".format(
                        transcript_i_id,
                        transcript_i_simple_path,
                        transcript_j_id,
                        transcript_j_simple_path,
                    )
                )

                combined_frac_gene_expression = (
                    frac_gene_expression_j + frac_gene_expression_i
                )
                # Relative support is undefined when neither model carries assigned
                # reads, so such a pair is left out of expression-based pruning.
                frac_gene_expression_j_of_i_and_j = (
                    frac_gene_expression_j / combined_frac_gene_expression
                    if combined_frac_gene_expression > 0
                    else None
                )

                ##
                ## Evaluate containments
                ##
                ##
                ## variables of interest: (in LRAA_Globals.config[]
                #  TSS:
                #    - "max_frac_alt_TSS_from_degradation"
                #    - "min_frac_gene_alignments_define_TSS_site"
                #  PolyA:
                #    - "min_frac_alignments_define_polyA_site"

                paths_indicate_containment = SPU.path_A_contains_path_B(
                    i_path_trimmed, j_path_trimmed
                ) or SPU.are_overlapping_and_compatible_NO_gaps_in_overlap(
                    i_path_trimmed, j_path_trimmed
                )

                if paths_indicate_containment:

                    logger.debug(
                        "splice compatible & contained transcript_j_id {} has frac gene expression: {}".format(
                            transcript_j_id, frac_gene_expression_j
                        )
                    )

                    subsume_J = False  # init

                    # no TSS or PolyA on j - prune.
                    if j_TSS_id is None and j_polyA_id is None:
                        logger.debug(
                            "compatible/contained {} being pruned as lacking TSS or polyA annots.".format(
                                transcript_j_id
                            )
                        )
                        subsume_J = True

                    # no TSS on j and shares same PolyA as i - prune as redundant 5' truncation
                    elif (
                        j_TSS_id is None
                        and j_polyA_id is not None
                        and j_polyA_id == i_polyA_id
                    ):
                        logger.debug(
                            "compatible/contained {} being pruned as lacking TSS and sharing PolyA with {}".format(
                                transcript_j_id, transcript_i_id
                            )
                        )
                        subsume_J = True

                    # no PolyA on j and shares same TSS as i - prune as redundant 3' truncation
                    elif (
                        j_polyA_id is None
                        and j_TSS_id is not None
                        and j_TSS_id == i_TSS_id
                    ):
                        logger.debug(
                            "compatible/contained {} being pruned as lacking PolyA and sharing TSS with {}".format(
                                transcript_j_id, transcript_i_id
                            )
                        )
                        subsume_J = True

                    elif (
                        frac_gene_expression_j_of_i_and_j is not None
                        and frac_gene_expression_j_of_i_and_j
                        < LRAA_Globals.config["max_rel_frac_expr_alt_compat_contained"]
                    ):
                        # if relative fraction of support for both is below threshold, then prune.
                        logger.debug(
                            "compatible/contained {} being pruned as insufficiently relatively expressed.".format(
                                transcript_j_id
                            )
                        )
                        subsume_J = True

                    if subsume_J:
                        logger.debug(
                            "Pruning {} as likely degradation product of {}".format(
                                transcript_j, transcript_i
                            )
                        )
                        transcript_prune_as_degradation.add(transcript_j)
                        transcript_i.absorb_other_transcript_multipaths(transcript_j)

        # retain the ones not pruned
        for transcript in transcript_set:
            if transcript not in transcript_prune_as_degradation:
                transcripts_ret.append(transcript)
            else:
                if LRAA_Globals.DEBUG:
                    logger.debug(
                        "FILTERING transcript {} as a likely degradation product".format(
                            transcript
                        )
                    )

    return transcripts_ret


def filter_internally_primed_transcripts(
    transcripts,
    contig_seq_str,
    contig_strand,
    known_transcripts,
    restrict_filter_to_monoexonic,
    spare_monoexonic_with_known_3prime=False,
):

    # Build prefix using contig_strand and contig accession from first transcript if available
    try:
        _ca = transcripts[0].get_contig_acc() if transcripts else None
        _prefix = f"[{_ca}{contig_strand}] " if _ca and contig_strand else ""
    except Exception:
        _prefix = ""

    logger.info(f"{_prefix}filtering internally primed. Restricting to monoexonic: {restrict_filter_to_monoexonic}")

    known_polyA_dist_ok_window = LRAA_Globals.config["max_dist_between_alt_polyA_sites"]
    known_polyA_dist_ok_window_half = int(known_polyA_dist_ok_window / 2)

    # build a list of known/acceptable 3' ends that get a free pass
    known_ok_3prime_ends = set()
    if known_transcripts is not None:
        for known_transcript in known_transcripts:
            if known_transcript.get_strand() == contig_strand:
                transcript_lend, transcript_rend = known_transcript.get_coords()
                known_3prime_end = transcript_rend if contig_strand == "+" else transcript_lend
                known_ok_3prime_ends.add(known_3prime_end)

    known_ok_3prime_ends_itree = itree.IntervalTree()
    for ok_3prime_end in known_ok_3prime_ends:
        known_ok_3prime_ends_itree[
            max(1, ok_3prime_end - known_polyA_dist_ok_window_half) : (
                ok_3prime_end + known_polyA_dist_ok_window_half
            )
        ] = ok_3prime_end

    retained_transcripts = list()

    for transcript in transcripts:

        # evaluate whether transcript looks internally primed.
        transcript_lend, transcript_rend = transcript.get_coords()
        strand = transcript.get_orient()

        looks_internally_primed = _looks_internally_primed(
            transcript_lend, transcript_rend, strand, contig_seq_str
        )

        try:
            logger.info(
                f"{_prefix}Transcript {transcript} looks internally primed? {looks_internally_primed}"
            )
        except Exception:
            logger.info(
                "Transcript {} looks internally primed? {}".format(
                    transcript, looks_internally_primed
                )
            )

        # persist evaluation result (may later be overridden in special monoexonic restriction case)
        transcript.set_likely_internal_primed(looks_internally_primed)

        if restrict_filter_to_monoexonic and not transcript.is_monoexonic():
            # In this mode we retain multi-exonic transcripts but preserve the
            # evaluated annotation state for downstream GTF export.
            transcript.set_likely_internal_primed(looks_internally_primed)
            retained_transcripts.append(transcript)
            continue

        filter_flag = False

        # A 3' end that agrees with a reference annotation argues the polyA signal is
        # genuine rather than internal priming. Multi-exonic models have always been
        # spared on that basis. Monoexonic models are spared only when
        # spare_monoexonic_internal_priming_with_known_3prime is enabled, since a
        # monoexonic model has no intron chain corroborating it and is the more likely
        # internal-priming artifact.
        if looks_internally_primed:

            end_check = transcript_rend if contig_strand == "+" else transcript_lend
            agrees_with_known_3prime = (
                len(known_ok_3prime_ends_itree[end_check : end_check + 1]) > 0
            )

            if transcript.is_monoexonic() and not spare_monoexonic_with_known_3prime:
                logger.debug(
                    "FILTERING monoexonic transcript {} as likely internally primed".format(
                        transcript
                    )
                )
                filter_flag = True

            elif not agrees_with_known_3prime:
                logger.debug(
                    "FILTERING {}exonic transcript {} as likely internally primed".format(
                        "mono" if transcript.is_monoexonic() else "multi", transcript
                    )
                )
                filter_flag = True

            else:
                logger.info(
                    "Ignoring internal priming info for {} as found consistent with known 3' end".format(
                        transcript
                    )
                )

        if not filter_flag:
            # keep
            retained_transcripts.append(transcript)

    return retained_transcripts


def _looks_internally_primed(
    transcript_lend, transcript_rend, strand, contig_seq_str, check_length=20
):

    if strand not in {"+", "-"}:
        raise ValueError("Strand must be '+' or '-'")

    target_base = "A" if strand == "+" else "T"
    target_polyA_motif = target_base * 6

    contig_length = len(contig_seq_str)

    if strand == "+":
        start = transcript_rend + 1
        end = start + check_length - 1
    else:
        end = transcript_lend - 1
        start = end - check_length + 1

    # ensure coordinates within contig bounds
    start = max(1, start)
    end = min(end, contig_length)

    extracted_long_sequence = contig_seq_str[start - 1 : end].upper()
    extracted_short_sequence = (
        extracted_long_sequence[-8:] if strand == "-" else extracted_long_sequence[0:8]
    )

    has_flanking_polyA = (
        extracted_long_sequence.count(target_base) >= 12
        or target_polyA_motif in extracted_short_sequence
    )

    return has_flanking_polyA


def evaluate_splice_compatible_alt_isoforms(transcripts):

    transcript_id_to_splice_compatible_containments = defaultdict(set)
    transcript_id_to_splice_compatible_contained_by = defaultdict(set)

    # Prefix using transcript context if available
    try:
        _t = transcripts[0] if transcripts else None
        _ca = _t.get_contig_acc() if _t else None
        _cs = _t.get_strand() if _t else None
        _prefix = f"[{_ca}{_cs}] " if _ca and _cs else ""
    except Exception:
        _prefix = ""

    logger.info(f"{_prefix}-evaluationg splice compatible alt isoforms:")

    if len(transcripts) < 2:
        return (
            transcript_id_to_splice_compatible_containments,
            transcript_id_to_splice_compatible_contained_by,
        )

    if type(transcripts) == set:
        transcripts = list(transcripts)

    transcripts.sort(key=lambda x: x.get_TPM(), reverse=True)

    for i in range(len(transcripts) - 1):

        transcript_i = transcripts[i]
        transcript_i_sp = transcript_i.get_simple_path()
        transcript_i_introns = SPU.get_simple_path_introns(transcript_i_sp)
        transcript_i_id = transcript_i.get_transcript_id()

        if len(transcript_i_introns) == 0:
            continue

        for j in range(i + 1, len(transcripts)):
            transcript_j = transcripts[j]
            transcript_j_sp = transcript_j.get_simple_path()
            transcript_j_introns = SPU.get_simple_path_introns(transcript_j_sp)
            transcript_j_id = transcript_j.get_transcript_id()

            if len(transcript_j_introns) == 0:
                continue

            if len(transcript_j_introns - transcript_i_introns) == 0:

                transcript_id_to_splice_compatible_containments[transcript_i_id].add(
                    transcript_j_id
                )
                transcript_id_to_splice_compatible_contained_by[transcript_j_id].add(
                    transcript_i_id
                )

                try:
                    logger.info(
                        f"{_prefix}splice-compatible: {transcript_j_id} (TPM={transcript_j.get_TPM():.3f}, "
                        f"iso_frac={transcript_j.get_isoform_fraction():.3f}) contained by "
                        f"{transcript_i_id} (TPM={transcript_i.get_TPM():.3f}, "
                        f"iso_frac={transcript_i.get_isoform_fraction():.3f})"
                    )
                except Exception:
                    logger.debug(
                        "Splice compatible isoforms: {} expr: {} splice compatible with {} expr: {}".format(
                            transcript_i,
                            transcript_i.get_TPM(),
                            transcript_j,
                            transcript_j.get_TPM(),
                        )
                    )

            if len(transcript_i_introns - transcript_j_introns) == 0:

                transcript_id_to_splice_compatible_containments[transcript_j_id].add(
                    transcript_i_id
                )
                transcript_id_to_splice_compatible_contained_by[transcript_i_id].add(
                    transcript_j_id
                )

                try:
                    logger.info(
                        f"{_prefix}splice-compatible: {transcript_i_id} (TPM={transcript_i.get_TPM():.3f}, "
                        f"iso_frac={transcript_i.get_isoform_fraction():.3f}) contained by "
                        f"{transcript_j_id} (TPM={transcript_j.get_TPM():.3f}, "
                        f"iso_frac={transcript_j.get_isoform_fraction():.3f})"
                    )
                except Exception:
                    logger.debug(
                        "Splice compatible isoforms: {} expr: {} splice compatible with {} expr: {}".format(
                            transcript_j,
                            transcript_j.get_TPM(),
                            transcript_i,
                            transcript_i.get_TPM(),
                        )
                    )

    return (
        transcript_id_to_splice_compatible_containments,
        transcript_id_to_splice_compatible_contained_by,
    )


def filter_novel_isoforms_by_min_read_support(
    transcripts, min_reads_novel_isoform: int
):

    retained_transcripts = list()

    # Build prefix from first transcript if available
    try:
        _t = transcripts[0] if transcripts else None
        _ca = _t.get_contig_acc() if _t else None
        _cs = _t.get_strand() if _t else None
        _prefix = f"[{_ca}{_cs}] " if _ca and _cs else ""
    except Exception:
        _prefix = ""

    for transcript in transcripts:
        if transcript.has_novel_splice_pattern() is True:
            try:
                logger.info(
                    f"{_prefix}transcript {transcript} is a novel isoform with {transcript.get_read_counts_assigned()} read support"
                )
            except Exception:
                logger.info("transcript {} is a novel isoform with {} read support".format(transcript, transcript.get_read_counts_assigned()))
            if transcript.get_read_counts_assigned() >= min_reads_novel_isoform:
                retained_transcripts.append(transcript)
            else:
                # novel transcript gets pruned as insufficient evidence.
                if LRAA_Globals.DEBUG:
                    logger.debug(
                        "FILTERING {} as insufficient evidence: read_support {} vs. min required {}".format(
                            transcript,
                            transcript.get_read_counts_assigned(),
                            min_reads_novel_isoform,
                        )
                    )
                pass
        else:
            # known transcript, retaining.
            try:
                logger.info(f"{_prefix}transcript {transcript} is a known isoform with {transcript.get_read_counts_assigned()} read support")
            except Exception:
                logger.info("transcript {} is a known isoform with {} read support".format(transcript, transcript.get_read_counts_assigned()))
            retained_transcripts.append(transcript)

    return retained_transcripts
