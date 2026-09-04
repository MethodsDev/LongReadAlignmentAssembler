#!/usr/bin/env python3

import math
import sys, os, re
import Transcript
import MultiPath
import MultiPathCounter
import Simple_path_utils as SPU
import Util_funcs
from Quantify import Quantify
from collections import defaultdict
import LRAA_Globals
from LRAA_Globals import SPACER, DEBUG
import logging
from math import log
import intervaltree as itree

logger = logging.getLogger(__name__)


# Set LRAA_FILTER_DECISION_LOG to a path prefix to record one row per
# isoform-filtering decision. Written straight to a file rather than through the
# logger, because isoform filtering runs inside worker processes whose logger
# output never reaches the driver's stdout -- piping the driver silently collects
# nothing, which is how an earlier attempt produced an empty file after 23
# minutes of compute.
#
# One file per process rather than a shared append. O_APPEND on a regular file
# gives no cross-process guarantee of whole-line writes, and buffered text I/O
# can split a line at an arbitrary byte; concatenate the parts afterwards.
_FILTER_DECISION_LOG = os.environ.get("LRAA_FILTER_DECISION_LOG") or None

_FILTER_DECISION_FIELDS = (
    "verdict transcript_id unique_reads gene_reads frac_unique "
    "isoform_fraction FSM_reads n_introns carriers_standing novel_splice "
    "near_unique effectively_unique "
    "contig strand filtering_round chain"
).split()

_filter_decision_ofh = None
_filter_decision_failed = False


def _record_filter_decision(*values):
    global _filter_decision_ofh, _filter_decision_failed
    if _filter_decision_failed:
        return
    try:
        if _filter_decision_ofh is None:
            _filter_decision_ofh = open(
                "{}.{}.tsv".format(_FILTER_DECISION_LOG, os.getpid()), "w"
            )
            _filter_decision_ofh.write("\t".join(_FILTER_DECISION_FIELDS) + "\n")
        _filter_decision_ofh.write("\t".join(str(v) for v in values) + "\n")
        _filter_decision_ofh.flush()
    except Exception as e:
        # A diagnostic must not take down a run, but it must not fail silently
        # either -- an empty diagnostic that looks like "no decisions" is worse
        # than none at all.
        _filter_decision_failed = True
        logger.warning("filter-decision recording disabled after error: %s", e)

_PREFILTER_FIELDS = (
    "verdict transcript_id contig strand reads_assigned novel_splice"
).split()
_prefilter_ofh = None
_prefilter_failed = False


def _record_prefilter_decision(*values):
    """Log for min_reads_novel_isoform, which is NOT the same gate as the in-loop
    min_unique_reads_novel_isoform: it runs once before any iterative filtering and
    tests get_read_counts_assigned(), a FRACTIONAL total shared with competing
    isoforms, not a literal count of reads unique to this model. Kept in its own
    file so the two floors are never summed or confused.
    """
    global _prefilter_ofh, _prefilter_failed
    if _prefilter_failed or _FILTER_DECISION_LOG is None:
        return
    try:
        if _prefilter_ofh is None:
            _prefilter_ofh = open(
                "{}.prefilter.{}.tsv".format(_FILTER_DECISION_LOG, os.getpid()), "w"
            )
            _prefilter_ofh.write("\t".join(_PREFILTER_FIELDS) + "\n")
        _prefilter_ofh.write("\t".join(str(v) for v in values) + "\n")
        _prefilter_ofh.flush()
    except Exception as e:
        _prefilter_failed = True
        logger.warning("prefilter recording disabled after error: %s", e)


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


def filter_multiexonic_isoforms_by_TPM_threshold(transcripts, min_TPM):
    """Drop multi-exonic isoforms that quantification found no expression for.

    Supplied models reach path selection carrying a synthetic template read, so a
    reference structure is selectable even when a dominant overlapping isoform has
    already claimed every real read it could have had. That is deliberate: it defers
    the judgement to quantification, which can see read assignments the greedy
    selection order cannot. This is where the judgement is made.

    Retains a transcript when its TPM is strictly greater than min_TPM, so the
    default of 0 means "keep it if EM gave it any expression at all", matching what
    ref_trans_filter_mode=retain_expressed asks everywhere else. A structure no read
    supports quantifies to zero and is dropped here rather than being reported on the
    strength of its own annotation.

    Applies to every multi-exonic model, not only reference-containing ones: a novel
    model with no expression is no better evidenced than a supplied one. Monoexonic
    models are untouched; they have their own thresholds. A negative min_TPM disables
    the filter entirely.
    """

    if min_TPM is None or min_TPM < 0:
        return transcripts

    transcripts_retained = []
    for transcript in transcripts:
        if transcript.is_monoexonic():
            transcripts_retained.append(transcript)
            continue
        if transcript.get_TPM() > min_TPM:
            transcripts_retained.append(transcript)

    return transcripts_retained


def get_splice_pattern(transcript):
    """The isoform's intron chain, as splice-graph nodes.

    Introns are nodes of the gene's own graph, so identity of node sequence is
    identity of intron chain here. Reads reach this code as multipaths through
    the same graph, which makes the comparison exact without resolving either
    side back to coordinates. Pretty_alignment.get_splice_hashcode carries the
    coordinate form for reporting, where stability across graphs matters.
    """
    return tuple(
        node for node in transcript.get_simple_path() if node.startswith("I:")
    )


def get_isoform_FSM_read_count(transcript, frac_read_assignments):
    """Reads in this isoform's support whose splice pattern is exactly its own.

    A read is a full splice match when its intron chain is identical to the
    isoform's -- not a prefix, not a subset, not merely compatible. That is
    direct evidence for the whole structure, and neither test above can see
    it: isoform fraction measures how much of a gene EM apportioned here, and
    apportionment tracks abundance rather than evidence, while the unique-read
    fraction is measured against gene depth and so demands more of a minor
    isoform the better expressed its gene is.

    Counted within this isoform's own support, so a read that fits only some
    other terminal variant of the same intron chain cannot vouch for this one.
    Unioning across variants would also make the count depend on which of them
    had already been filtered when it was taken.

    The fraction EM assigned is deliberately ignored. It says which isoform
    currently explains a read best and is recomputed every filtering round, so
    gating on it would let a threshold protect a model in one round and drop
    it the next. Whether a read carries this intron chain does not change;
    whether EM prefers a neighbour does.
    """
    pattern = get_splice_pattern(transcript)
    if not pattern:
        return 0  # a monoexonic model has no splice pattern to match

    reads_by_path = dict()
    for mp in frac_read_assignments.get(transcript.get_transcript_id(), {}):
        if (
            tuple(node for node in mp.get_simple_path() if node.startswith("I:"))
            != pattern
        ):
            continue
        reads_by_path[mp.get_simple_path_tuple()] = mp.get_read_count()

    return sum(reads_by_path.values())


def has_insufficient_support(
    num_FSM_reads,
    frac_gene_unique_reads,
    is_spliced,
    min_FSM_reads_gate,
    min_frac_gene_unique_reads,
):
    """Whether a model is too weakly supported to keep.

    By default this is the isoform's uniquely-assigned reads as a fraction of the
    gene's total, which is relative to gene depth: a minor isoform of a deeply
    sequenced gene must clear a bar that rises with its neighbours. Setting
    min_FSM_reads_gate substitutes an absolute count of reads whose splice pattern
    is exactly this model's -- direct evidence for the structure rather than a
    share of the locus.

    Only for spliced models. A monoexonic model has no splice pattern, so its
    full-splice-match count is zero by definition and substituting it would delete
    every one of them outright; those keep the fraction.
    """
    if min_FSM_reads_gate > 0 and is_spliced:
        return num_FSM_reads < min_FSM_reads_gate

    return frac_gene_unique_reads < min_frac_gene_unique_reads


def filter_isoforms_by_min_isoform_fraction(
    transcripts, min_isoform_fraction, run_EM, max_EM_iterations
):

    min_frac_gene_unique_reads = LRAA_Globals.config["min_frac_gene_unique_reads"]
    min_FSM_reads_retain_isoform = LRAA_Globals.config["min_FSM_reads_retain_isoform"]
    min_FSM_reads_gate = LRAA_Globals.config["min_FSM_reads_gate"]
    min_unique_reads_novel_isoform = LRAA_Globals.config[
        "min_unique_reads_novel_isoform"
    ]
    _recording = _FILTER_DECISION_LOG is not None
    # The compatible-map inversion is needed on every path now: the novel-isoform
    # floor and frac_gene_unique_reads both read the unique counts, so there is no
    # configuration in which it can be skipped.

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

    def build_multipath_isoform_counts(frac_read_assignments):
        """multipath -> number of isoforms of this gene it can be assigned to.

        EM.run_EM writes an mp_assignments entry for every structurally compatible
        transcript, including ones whose E-step fraction is zero under the current
        theta, so membership in the map is the statement "this read could be assigned
        here". A read is UNIQUE to an isoform when exactly one map contains it.
        """
        n_compat = defaultdict(int)
        for mp_to_frac in frac_read_assignments.values():
            for mp in mp_to_frac:
                n_compat[mp] += 1
        return n_compat

    def get_isoform_unique_support(transcript_id, frac_read_assignments, n_compat):
        """Support from reads that can ONLY be assigned to this isoform.

        "Uniquely assigned" means exclusively assignable -- one compatible isoform --
        and that single definition now feeds every consumer: the literal count used by
        the novel-isoform floor, and the weighted share used by
        frac_gene_unique_reads. Previously both came from a fraction threshold
        (unique_read_filter_min_frac, 0.9995), which asks whether EM happened to
        concentrate a read's mass and so admits reads compatible with a competitor.

        The two returned magnitudes still differ on purpose. The WEIGHTED figure
        estimates how much of the original library this isoform uniquely explains and
        belongs in any ratio: frac_gene_unique_reads divides it by a gene total that
        Quantify.get_gene_read_counts also builds from get_read_weight(), so numerator
        and denominator agree. The LITERAL figure counts observations, because "two
        unique reads" is a confidence statement about having seen a structure twice.

        near_unique / effectively_unique are returned for the decision log only, so the
        divergence between this definition and the old threshold stays observable.
        """
        num_unique_reads = 0
        weighted_unique_support = 0.0
        near_unique = 0
        effectively_unique = 0
        gate = LRAA_Globals.config["unique_read_filter_min_frac"]
        for mp, frac in frac_read_assignments[transcript_id].items():
            n = mp.get_read_count()
            if n_compat[mp] == 1:
                num_unique_reads += n
                weighted_unique_support += mp.get_read_weight()
            if frac >= gate:
                near_unique += n
                if n_compat[mp] > 1:
                    # clears the old threshold yet is NOT exclusively assignable
                    effectively_unique += n
        return (
            num_unique_reads,
            weighted_unique_support,
            near_unique,
            effectively_unique,
        )

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

            n_compat_mp = build_multipath_isoform_counts(frac_read_assignments)

            gene_id_to_read_count = Quantify.get_gene_read_counts(
                frac_read_assignments, transcript_id_to_transcript_obj
            )

            gene_read_count = gene_id_to_read_count[gene_id]

            num_isoforms_of_gene_filtered = 0

            # Bare get_isoform_fraction() is a float with no tie-break, and this
            # is a weakest-first greedy removal: two isoforms on the same
            # fraction must not be able to swap places.
            isoforms_of_gene = sorted(
                isoforms_of_gene,
                key=lambda x: (
                    x.get_isoform_fraction(),
                    Transcript.Transcript.structural_sort_key(x),
                ),
            )
            num_isoforms_of_gene = len(isoforms_of_gene)

            # Isoforms carrying each splice pattern, and how many are still
            # standing. Filtering runs weakest-first, so carriers fall away in
            # ascending order of isoform fraction and the last one left is the
            # strongest. The tally has to be live rather than taken at round
            # start: the aggressive branch below removes several isoforms in a
            # single pass, which can take every carrier of a pattern at once.
            pattern_to_carriers = defaultdict(list)
            for isoform in isoforms_of_gene:
                pattern = get_splice_pattern(isoform)
                if pattern:
                    pattern_to_carriers[pattern].append(isoform)
            pattern_carriers_standing = {
                pattern: len(carriers)
                for pattern, carriers in pattern_to_carriers.items()
            }

            for transcript in isoforms_of_gene:

                num_total_isoforms += 1
                transcript_id = transcript.get_transcript_id()

                (
                    transcript_unique_read_count,
                    transcript_unique_weighted_support,
                    near_unique_reads,
                    effectively_uniq_reads,
                ) = get_isoform_unique_support(
                    transcript_id, frac_read_assignments, n_compat_mp
                )

                # Weighted over weighted: gene_read_count is itself built from
                # get_read_weight() (Quantify.get_gene_read_counts), so the literal count
                # must not be the numerator here or the ratio is deflated by the
                # acceptance rate. The literal count is what the absolute gate below uses.
                frac_gene_unique_reads = (
                    transcript_unique_weighted_support / gene_read_count
                    if gene_read_count > 0
                    else 0
                )

                logger.debug(
                    "Transcript_id: {} has unique read frac of gene total reads: {}".format(
                        transcript_id, frac_gene_unique_reads
                    )
                )

                transcript_pattern = get_splice_pattern(transcript)

                # Reads whose splice pattern is exactly this model's. Only needed
                # when the exemption is enabled or when recording decisions.
                recording = _FILTER_DECISION_LOG is not None
                num_FSM_reads = (
                    get_isoform_FSM_read_count(transcript, frac_read_assignments)
                    if (
                        min_FSM_reads_retain_isoform > 0
                        or min_FSM_reads_gate > 0
                        or recording
                    )
                    else 0
                )

                # Which quantity decides that a model is too weakly supported.
                #
                # By default it is the isoform's uniquely-assigned reads as a
                # fraction of the gene's total, which is relative to gene depth:
                # a minor isoform of a deeply sequenced gene must clear a bar
                # that rises with its neighbours. Setting min_FSM_reads_gate
                # substitutes an absolute count of reads whose splice pattern is
                # exactly this model's -- direct evidence for the structure
                # rather than a share of the locus.
                #
                # Measured on chr20 ref-guided, gating at 2 full-splice-match
                # reads instead: precision 0.329 -> 0.363, 268 fewer false chains
                # for 10 fewer true ones. Tightening the existing fraction to
                # 0.02 reaches the same precision but costs 58 more true chains,
                # so the improvement comes from the quantity and not the cut.
                # Ranked by permutation importance over held-out gene folds the
                # unique-read fraction places tenth while this count places
                # first; the filter currently uses the former and not the latter.
                # Only for spliced models: a monoexonic model has no splice
                # pattern, so its full-splice-match count is zero by definition
                # and substituting it would delete every one of them outright.
                # Those keep the fraction.
                insufficient_support = has_insufficient_support(
                    num_FSM_reads,
                    frac_gene_unique_reads,
                    bool(transcript_pattern),
                    min_FSM_reads_gate,
                    min_frac_gene_unique_reads,
                )

                # A splice pattern that reads attest to end to end keeps one model:
                # the last still standing, which is the strongest because filtering
                # runs weakest-first.
                retain_on_FSM_evidence = (
                    min_FSM_reads_retain_isoform > 0
                    and transcript_pattern
                    and pattern_carriers_standing[transcript_pattern] == 1
                    and num_FSM_reads >= min_FSM_reads_retain_isoform
                )

                if (
                    transcript.contains_reference_model()
                    and LRAA_Globals.config["ref_trans_filter_mode"]
                    == "retain_expressed"
                    and transcript.get_TPM() > 0
                ):
                    transcripts_retained.append(transcript)
                    verdict = "retained_reference_model"

                elif retain_on_FSM_evidence:
                    transcripts_retained.append(transcript)
                    verdict = "retained_FSM_evidence"

                ## if tons of isoforms, allow pruning of multiple in a single round
                elif (
                    num_isoforms_of_gene - num_isoforms_of_gene_filtered
                    > LRAA_Globals.config[
                        "min_isoform_count_aggressive_filtering_iso_fraction"
                    ]
                    and insufficient_support
                    and transcript.get_isoform_fraction() < min_isoform_fraction
                ):
                    verdict = "filtered_aggressive"

                # standard isoform fraction based filtering
                elif not isoforms_were_filtered and (
                    insufficient_support
                    or transcript.get_isoform_fraction() < min_isoform_fraction
                ):
                    verdict = "filtered_low_unique_or_isoform_fraction"

                elif not isoforms_were_filtered and (
                    transcript.has_novel_splice_pattern() is True
                    and transcript_unique_read_count
                    < min_unique_reads_novel_isoform
                ):
                    verdict = "filtered_novel_too_few_unique_reads"

                else:
                    transcripts_retained.append(transcript)
                    verdict = "retained"

                if verdict.startswith("filtered"):
                    isoforms_were_filtered = True
                    num_filtered_isoforms += 1
                    num_isoforms_of_gene_filtered += 1
                    if transcript_pattern:
                        pattern_carriers_standing[transcript_pattern] -= 1

                    # kept verbatim: existing tooling and analyses grep these
                    if verdict == "filtered_low_unique_or_isoform_fraction":
                        logger.debug(
                            "Filtering out transcript_id {} as low fraction of unique reads: {}".format(
                                transcript_id, frac_gene_unique_reads
                            )
                        )
                    elif verdict == "filtered_novel_too_few_unique_reads":
                        logger.debug(
                            "Filtering out transcript_id {} as novel isoform with too few unique reads: {}".format(
                                transcript_id, transcript_unique_read_count
                            )
                        )

                if recording:
                    _record_filter_decision(
                        verdict,
                        transcript_id,
                        transcript_unique_read_count,
                        gene_read_count,
                        frac_gene_unique_reads,
                        transcript.get_isoform_fraction(),
                        num_FSM_reads,
                        len(transcript_pattern) if transcript_pattern else 0,
                        (
                            pattern_carriers_standing.get(transcript_pattern, 0)
                            if transcript_pattern
                            else 0
                        ),
                        transcript.has_novel_splice_pattern(),
                        near_unique_reads,
                        effectively_uniq_reads,
                        # Coordinate chain plus contig/strand and round: the node-tuple
                        # form (get_splice_pattern) is only identity-stable within one
                        # gene's graph, so it must not key cross-run joins. The round is
                        # required because the same transcript_id is reconsidered as EM
                        # support changes, and a bare id would silently collapse rounds.
                        transcript.get_contig_acc(),
                        transcript.get_strand(),
                        filtering_round,
                        ",".join(
                            "{}-{}".format(lend, rend)
                            for lend, rend in transcript.get_introns()
                        )
                        if transcript.get_introns()
                        else "monoexonic",
                    )

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

    # gene_id -> list of Transcript, in the caller's order.  Was a
    # defaultdict(set); Transcript has no __hash__, so that set iterated in
    # memory-address order and every consumer of this mapping -- the
    # weakest-first isoform-fraction filter and the degradation-product pruner,
    # both of which decide which transcript survives -- inherited it.
    gene_id_to_transcripts = defaultdict(dict)
    for transcript in transcripts:
        gene_id = transcript.get_gene_id()
        gene_id_to_transcripts[gene_id][transcript.get_transcript_id()] = transcript

    return {
        gene_id: list(by_id.values())
        for gene_id, by_id in gene_id_to_transcripts.items()
    }


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


def annotate_polyA_signal(transcripts, contig_seq_str, contig_strand):
    """Record the canonical PAS upstream of each model's own 3' terminus.

    Pure annotation: nothing is filtered, reordered or removed, and the pass is not
    gated on any config key, so the attribute is present on every emitted model.

    Deliberately separate from internal-priming filtering. That asks whether the
    genome DOWNSTREAM of the terminus looks like an oligo-dT template; this asks
    whether the signal a genuine cleavage site requires is present UPSTREAM. The two
    are independent evidence about the same site and a terminus can carry both, so
    they are reported as separate attributes rather than folded into one verdict.
    """

    try:
        _ca = transcripts[0].get_contig_acc() if transcripts else None
        _prefix = f"[{_ca}{contig_strand}] " if _ca and contig_strand else ""
    except Exception:
        _prefix = ""

    motifs, window = Util_funcs.resolve_polyA_signal_settings()

    n_with_signal = 0
    for transcript in transcripts:
        transcript_lend, transcript_rend = transcript.get_coords()
        strand = transcript.get_orient()
        three_prime_pos = transcript_rend if strand == "+" else transcript_lend

        hexamer, offset = Util_funcs.find_polyA_signal(
            contig_seq_str, three_prime_pos, strand, window=window, hexamers=motifs
        )
        n_with_signal += hexamer is not None
        transcript.set_polyA_signal(hexamer, offset)

    logger.info(
        f"{_prefix}polyA signal {list(motifs)} in window {list(window)} upstream of the "
        f"3' end: {n_with_signal} of {len(transcripts)} transcripts"
    )

    return transcripts


def filter_internally_primed_transcripts(
    transcripts,
    contig_seq_str,
    contig_strand,
    known_transcripts,
    restrict_filter_to_monoexonic,
    spare_monoexonic_with_known_3prime=False,
    delete=True,
):
    """Annotate every transcript's 3' terminus for A-rich genomic context, and delete
    the models that policy says should not survive it.

    Two stages now apply the same rule, deliberately and for different reasons.
    `Splice_graph._incorporate_PolyA_objects` rejects an A-rich READ-DERIVED candidate
    before it can become a graph vertex, so the constraint is never created. This stage
    then judges the model's own emitted 3' terminus, which reconstruction may have
    placed somewhere no rejected candidate ever sat.

    Both are wanted. On chr20 the gate alone leaves 123 monoexonic models terminating
    in A-rich context, and they are an artifact population rather than real single-exon
    transcripts: 2.4% end within 25 bp of a measured PolyASite cleavage site against
    52.4% for unflagged monoexonic models, and 1 of 123 ends inside an annotated exon.

    The annotation is set on every transcript regardless of `delete`, because
    `Transcript.py` reads an `InternalPriming` attribute back off an input GTF in order
    to re-export it, so the attribute has to round-trip even when nothing is being
    removed. Consumers should read it as a statement about the emitted model's
    terminus, not about a rejected PolyA candidate.

    The reprieve in `spare_monoexonic_with_known_3prime` keys on proximity to a
    `known_transcripts` 3' end -- the reference annotation handed to this function --
    and not on any measured cleavage atlas. Note that its window is asymmetric by one
    base -- intervals go in as [K - half, K + half) but are probed at [E, E + 1), so a
    known end spares a model ending at E over [E - half + 1, E + half]. That reads as an
    off-by-one rather than an intended tolerance. It predates this change and is
    unreachable at shipped defaults, so it is left alone here rather than fixed in
    passing; the tests deliberately do not assert that edge.
    """

    try:
        _ca = transcripts[0].get_contig_acc() if transcripts else None
        _prefix = f"[{_ca}{contig_strand}] " if _ca and contig_strand else ""
    except Exception:
        _prefix = ""

    logger.info(
        f"{_prefix}internal priming: delete={delete}, "
        f"restricting deletion to monoexonic: {restrict_filter_to_monoexonic}"
    )

    known_polyA_dist_ok_window = LRAA_Globals.config["max_dist_between_alt_polyA_sites"]
    known_polyA_dist_ok_window_half = int(known_polyA_dist_ok_window / 2)

    known_ok_3prime_ends = set()
    if known_transcripts is not None:
        for known_transcript in known_transcripts:
            if known_transcript.get_strand() == contig_strand:
                transcript_lend, transcript_rend = known_transcript.get_coords()
                known_3prime_end = (
                    transcript_rend if contig_strand == "+" else transcript_lend
                )
                known_ok_3prime_ends.add(known_3prime_end)

    known_ok_3prime_ends_itree = itree.IntervalTree()
    for ok_3prime_end in known_ok_3prime_ends:
        known_ok_3prime_ends_itree[
            max(1, ok_3prime_end - known_polyA_dist_ok_window_half) : (
                ok_3prime_end + known_polyA_dist_ok_window_half
            )
        ] = ok_3prime_end

    retained_transcripts = list()
    n_flagged = 0
    n_deleted = 0

    for transcript in transcripts:

        transcript_lend, transcript_rend = transcript.get_coords()
        strand = transcript.get_orient()
        three_prime_pos = transcript_rend if strand == "+" else transcript_lend

        looks_internally_primed = Util_funcs.looks_internally_primed(
            contig_seq_str, three_prime_pos, strand
        )
        n_flagged += bool(looks_internally_primed)

        # Set on every transcript, deleted or not, so the attribute round-trips.
        transcript.set_likely_internal_primed(looks_internally_primed)

        if not delete:
            retained_transcripts.append(transcript)
            continue

        if restrict_filter_to_monoexonic and not transcript.is_monoexonic():
            retained_transcripts.append(transcript)
            continue

        filter_flag = False

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
            retained_transcripts.append(transcript)
        else:
            n_deleted += 1

    logger.info(
        f"{_prefix}internal priming: {n_flagged} of {len(transcripts)} transcripts "
        f"flagged, {n_deleted} deleted, {len(retained_transcripts)} retained"
    )

    return retained_transcripts


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
            _record_prefilter_decision(
                (
                    "prefilter_novel_retained"
                    if transcript.get_read_counts_assigned()
                    >= min_reads_novel_isoform
                    else "prefilter_novel_removed"
                ),
                transcript.get_transcript_id(),
                _ca,
                _cs,
                transcript.get_read_counts_assigned(),
                True,
            )
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
