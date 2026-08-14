#!/usr/bin/env python3

import sys, os, re
import Transcript
import MultiPath
import MultiPathCounter
import Simple_path_utils as SPU
from collections import defaultdict
import LRAA_Globals
from LRAA_Globals import SPACER, DEBUG
import logging
from math import log
from GenomeFeature import Exon
import Util_funcs
import EM
import time

logger = logging.getLogger(__name__)


class Quantify:

    def __init__(self, run_EM, max_EM_iterations, quant_mode="final"):

        self._run_EM = run_EM  # boolean
        self._max_EM_iterations = max_EM_iterations

        self._path_node_id_to_gene_ids = defaultdict(set)

        # transcript_id -> Transcript, per gene.  A dict rather than a set: no
        # class in pylib defines __hash__, so a set of Transcript objects hashes
        # by id() and iterates in memory-address order, which differs between
        # processes.  Keying by transcript id keeps set semantics (one entry per
        # transcript) while iterating in a stable insertion order.
        self._gene_id_to_transcript_objs = defaultdict(dict)

        self._read_name_to_multipath = dict()

        self._mp_to_transcripts = dict()
        self._unassigned_mp_count_pairs = list()

        # Component identity, for consumers that re-derive a per-read isoform
        # split.  Rebuilt per quantify() call rather than updated; see there.
        self._transcript_id_to_component_id = dict()
        self._component_id_to_gene_ids = dict()

        # Whether the two maps above describe a completed quantification.  A flag
        # rather than an emptiness check: a successful quantify() always yields at
        # least one component, so empty does mean invalid today, but that follows
        # from an assert 200 lines away and would stop being true the moment
        # someone relaxed it.  Stating it directly costs one line.
        self._component_identity_valid = False

        self._quant_mode = quant_mode

        return

    def quantify(self, splice_graph, transcripts, mp_counter):

        assert type(transcripts) == list
        assert type(transcripts[0]) == Transcript.Transcript
        assert type(mp_counter) == MultiPathCounter.MultiPathCounter

        # Every container below describes ONE quantification and is rebuilt by the
        # assignment steps in this call.  quantify() runs up to three times on one
        # object over a shrinking transcript set -- de novo calls it after
        # degradation pruning and again after isoform-fraction filtering -- so an
        # entry that survives describes transcripts this call does not contain.
        #
        # _mp_to_transcripts is the one that changes an answer.  Component
        # construction unions the genes of every retained entry, and the guard
        # there only skips genes absent from this call; a stale multipath bridging
        # two genes that both still have isoforms fuses them into one component
        # this call's reads no longer bridge.  It can only over-fuse, so theta is
        # normalized across genes that should not share a denominator.
        #
        # The gene indexes fuse nothing but keep filtered-out transcripts alive as
        # read-assignment candidates.  _unassigned_mp_count_pairs already reset
        # itself per call, which is the semantics the rest of these now match.
        self._path_node_id_to_gene_ids = defaultdict(set)
        self._gene_id_to_transcript_objs = defaultdict(dict)
        self._read_name_to_multipath = dict()
        self._mp_to_transcripts = dict()
        self._unassigned_mp_count_pairs = list()

        # Read after the call by consumers recomputing a per-read isoform split.
        # Cleared here rather than where they are filled so a call that raises
        # earlier answers empty: empty is a visible failure, stale is not.
        self._transcript_id_to_component_id = dict()
        self._component_id_to_gene_ids = dict()
        # Invalidated BEFORE any work that can fail, never after.  A raise between
        # here and the successful exit below must leave the accessors refusing, not
        # answering with the previous call's components.
        self._component_identity_valid = False

        contig_acc = splice_graph.get_contig_acc()
        contig_strand = splice_graph.get_contig_strand()

        # init transcript quant info
        gene_to_transcripts = defaultdict(list)
        for transcript in transcripts:
            transcript.init_quant_info()
            gene_id = transcript.get_gene_id()
            gene_to_transcripts[gene_id].append(transcript)

        try:
            logger.info(
                f"[{contig_acc}{contig_strand}] have {len(gene_to_transcripts)} genes and {len(transcripts)} isoforms to quantify."
            )
        except Exception:
            logger.info(
                "have {} genes and {} isoforms to quantify.".format(
                    len(gene_to_transcripts), len(transcripts)
                )
            )

        # assign path nodes to gene
        # also assign gene_id to transcript objs
        self._assign_path_nodes_to_gene(transcripts)

        self._assign_reads_to_transcripts(splice_graph, mp_counter)

        # The EM unit is a read-sharing component of genes, not a single gene: a
        # read compatible with transcripts in two genes has to be apportioned by
        # one converged EM rather than handed to each gene by two independent
        # ones.  Reporting below stays strictly gene-scoped.
        gene_components = self._build_read_sharing_gene_components(gene_to_transcripts)

        num_multigene_components = sum(1 for x in gene_components if len(x) > 1)
        if num_multigene_components > 0:
            logger.info(
                "[%s%s] %d of %d quant component(s) hold more than one gene; "
                "each such component is quantified as one joint EM.",
                contig_acc,
                contig_strand,
                num_multigene_components,
                len(gene_components),
            )

        # Component identity, exposed because theta is normalized over a whole
        # component and not over a gene.  A consumer that renormalizes per gene
        # gives a read compatible with transcripts of two genes in one component
        # a split summing to 2, which is arithmetically plausible and silent.
        # Both maps were cleared at entry above.
        for component_gene_ids in gene_components:
            # Leftmost gene of the component, made canonical by the sort in
            # _build_read_sharing_gene_components.  Components partition the
            # genes, so this is unique; treat it as an opaque handle.
            component_id = component_gene_ids[0]
            self._component_id_to_gene_ids[component_id] = list(component_gene_ids)
            for gene_id in component_gene_ids:
                for transcript in gene_to_transcripts[gene_id]:
                    self._transcript_id_to_component_id[
                        transcript.get_transcript_id()
                    ] = component_id

        transcript_to_fractional_read_assignment = dict()

        for component_gene_ids in gene_components:

            transcripts_list = list()
            for gene_id in component_gene_ids:
                transcripts_list.extend(gene_to_transcripts[gene_id])

            trans_coords = list()
            for transcript in transcripts_list:
                trans_coords.extend(transcript.get_coords())

            trans_coords = sorted(trans_coords)
            component_lend = trans_coords[0]
            component_rend = trans_coords[-1]

            # Single-gene components keep the original log line verbatim so an
            # ordinary run reads exactly as before; a joint EM is named as such,
            # with its member genes and combined span, because that is precisely
            # what someone chasing an unexpected number needs to see.
            if len(component_gene_ids) > 1:
                component_descr = "read-sharing component of {} genes {{{}}}".format(
                    len(component_gene_ids), ",".join(component_gene_ids)
                )
            else:
                component_descr = component_gene_ids[0]

            try:
                logger.info(
                    "[%s%s] quant estimates for isoforms of %s %s%s:%d-%d",
                    contig_acc,
                    contig_strand,
                    component_descr,
                    contig_acc,
                    contig_strand,
                    component_lend,
                    component_rend,
                )
            except Exception:
                logger.info(
                    "quant estimates for isoforms of %s %s%s:%d-%d",
                    component_descr,
                    contig_acc,
                    contig_strand,
                    component_lend,
                    component_rend,
                )

            # build a contig/strand prefix for downstream logs
            try:
                prefix_str = f"[{contig_acc}{contig_strand}] " if contig_acc and contig_strand else None
            except Exception:
                prefix_str = None

            # One joint EM over every transcript of every gene in the component.
            # _estimate_isoform_read_support() computes isoform fractions per
            # gene_id internally, so widening the EM does not widen them.
            component_transcript_to_fractional_read_assignment = (
                self._estimate_isoform_read_support(transcripts_list, prefix_str=prefix_str)
            )
            # copy over to the full data structure
            for transcript_id in component_transcript_to_fractional_read_assignment:
                transcript_to_fractional_read_assignment[transcript_id] = (
                    component_transcript_to_fractional_read_assignment[transcript_id]
                )

        # see documentation for _estimate_isoform_read_support() below

        # Only here: every path that reaches this point built the maps from this
        # call's components.
        self._component_identity_valid = True

        return transcript_to_fractional_read_assignment

    def _assign_path_nodes_to_gene(self, transcripts):

        for transcript in transcripts:

            simplepath = transcript._simplepath

            if simplepath is None:
                logger.warn(
                    "simplepath is not avaialble for transcript: {}".format(transcript)
                )
                continue

            # assert simplepath is not None, "Error, simplepath not set for transcript obj: {}".format(transcript)

            transcript_id = transcript.get_transcript_id()
            gene_id = transcript.get_gene_id()
            self._gene_id_to_transcript_objs[gene_id][transcript_id] = transcript

            for node_id in simplepath:
                if node_id != SPACER:
                    self._path_node_id_to_gene_ids[node_id].add(gene_id)

        return

    def _assign_reads_to_transcripts(
        self,
        splice_graph,
        mp_counter,
        fraction_read_align_overlap=None,
    ):
        if fraction_read_align_overlap is None:
            fraction_read_align_overlap = LRAA_Globals.config[
                "fraction_read_align_overlap"
            ]

        assert (
            fraction_read_align_overlap >= 0 and fraction_read_align_overlap <= 1.0
        ), "Error, fraction_read_align_overlap must be between 0 and 1.0"
        self._unassigned_mp_count_pairs = list()

        try:
            ca = splice_graph.get_contig_acc()
            cs = splice_graph.get_contig_strand()
        except Exception:
            ca, cs = None, None

        try:
            if ca is not None and cs is not None:
                logger.info(f"[{ca}{cs}] # Assigning reads to transcripts")
            else:
                logger.info("# Assigning reads to transcripts")
        except Exception:
            logger.info("# Assigning reads to transcripts")

        local_debug = False

        if local_debug is True:
            LRAA_orig_setting = LRAA_Globals.DEBUG
            logging_orig_setting = logging.DEBUG if LRAA_Globals.DEBUG else logging.INFO
            LRAA_Globals.DEBUG = True
            logging.getLogger().setLevel(logging.DEBUG)

        # assign to gene based on majority voting of nodes.
        # TODO:// might want or need this to involve length and/or feature type weighted shared node voting

        mp_count_pairs = mp_counter.get_all_MultiPathCountPairs()

        num_mp_count_pairs = len(mp_count_pairs)
        try:
            if ca is not None and cs is not None:
                logger.info(f"[{ca}{cs}] - have {num_mp_count_pairs} mp_count_pairs")
            else:
                logger.info("- have {} mp_count_pairs".format(num_mp_count_pairs))
        except Exception:
            logger.info("- have {} mp_count_pairs".format(num_mp_count_pairs))

        # progress monitoring configuration
        show_progress = LRAA_Globals.config.get("show_progress_quant_assign", True)
        prefer_tqdm = LRAA_Globals.config.get("use_tqdm_progress", True)
        progress_every_n = LRAA_Globals.config.get("progress_update_every_n", 1000)
        progress_interval = LRAA_Globals.config.get(
            "progress_update_interval_sec", 5.0
        )
        start_time = time.time()
        last_progress_time = start_time

        # initialize tqdm if available and preferred
        tqdm = None
        pbar = None
        if show_progress and prefer_tqdm and num_mp_count_pairs > 0:
            try:
                import importlib

                _mod = importlib.import_module("tqdm")
                tqdm = getattr(_mod, "tqdm", None)
                if tqdm is not None:
                    pbar = tqdm(
                        total=num_mp_count_pairs,
                        desc="quant-assign",
                        unit="mp",
                        leave=False,
                        dynamic_ncols=True,
                    )
            except Exception:
                tqdm = None
                pbar = None

        # Read paths that match no gene's splice-graph nodes are not assignable and
        # are not offered to rescue, so their reads leave quantification without an
        # output row. Count them so the loss is visible.
        num_paths_unanchored = 0
        num_read_counts_unanchored = 0

        num_paths_total = 0
        num_read_counts_total = 0

        num_paths_anchored_to_gene = 0
        num_read_counts_anchored_to_gene = 0

        num_paths_assigned = 0
        num_read_counts_assigned = 0

        mp_seen = set()

        num_mp_count_pairs_processed = 0

        def _maybe_report_progress():
            nonlocal last_progress_time
            # if tqdm is active, let it handle progress display
            if pbar is not None:
                return
            if not show_progress or num_mp_count_pairs == 0:
                return
            now = time.time()
            should_by_count = (
                progress_every_n is not None
                and progress_every_n > 0
                and num_mp_count_pairs_processed % progress_every_n == 0
            )
            should_by_time = (
                progress_interval is not None
                and progress_interval > 0
                and now - last_progress_time >= progress_interval
            )
            if should_by_count or should_by_time:
                elapsed = now - start_time
                frac = num_mp_count_pairs_processed / num_mp_count_pairs
                rate = num_mp_count_pairs_processed / elapsed if elapsed > 0 else 0.0
                remaining = max(num_mp_count_pairs - num_mp_count_pairs_processed, 0)
                eta_sec = remaining / rate if rate > 0 else float("inf")
                eta_txt = (
                    "{:02.0f}m{:02.0f}s".format(eta_sec // 60, eta_sec % 60)
                    if eta_sec != float("inf")
                    else "inf"
                )
                msg = (
                    f"\r[quant-assign] {num_mp_count_pairs_processed}/{num_mp_count_pairs} "
                    f"({frac*100:5.1f}%) | {rate:6.1f}/s | ETA {eta_txt}    "
                )
                # write to stderr to avoid interfering with stdout data files
                try:
                    sys.stderr.write(msg)
                    sys.stderr.flush()
                except Exception:
                    # fall back to logger in case stderr is unavailable
                    try:
                        if ca is not None and cs is not None:
                            logger.info(f"[{ca}{cs}] {msg.strip()}")
                        else:
                            logger.info(msg.strip())
                    except Exception:
                        pass
                last_progress_time = now

        for mp_count_pair in mp_count_pairs:

            num_mp_count_pairs_processed += 1
            if pbar is not None:
                try:
                    pbar.update(1)
                except Exception:
                    pass
            else:
                _maybe_report_progress()

            mp, count = mp_count_pair.get_multipath_and_count()

            mp_id = mp.get_id()

            if mp_id in mp_seen:
                raise RuntimeError("multipath already evaluated - error. " + str(mp))
            mp_seen.add(mp_id)

            num_paths_total += 1
            num_read_counts_total += count

            sp = mp.get_simple_path()

            top_genes = self._get_all_genes_with_node_matches_to_simplepath(sp)

            if top_genes is None:
                num_paths_unanchored += 1
                num_read_counts_unanchored += count
                logger.debug("mp_count_pair unanchored: " + str(mp_count_pair))

                continue

            logger.debug(
                "mp_count_pair {} anchored to genes: {}".format(
                    mp_count_pair, top_genes
                )
            )
            num_paths_anchored_to_gene += 1
            num_read_counts_anchored_to_gene += count

            ## assign reads to transcripts
            # Ordered dedupe across the (possibly overlapping) top genes; the
            # resulting sequence is what read-to-transcript compatibility is
            # built from, so it must not depend on address order.
            gene_isoforms_by_id = dict()
            for top_gene in top_genes:
                gene_isoforms_by_id.update(self._gene_id_to_transcript_objs[top_gene])
            gene_isoforms = list(gene_isoforms_by_id.values())

            logger.debug(
                "mp_count_pair {} assigned to gene {} with isoforms to test for read assignment:\n\t{}".format(
                    mp,
                    top_gene,
                    "\n\t".join(
                        [
                            "{}\t{}".format(x.get_transcript_id(), x._simplepath)
                            for x in gene_isoforms
                        ]
                    ),
                )
            )

            # most stringent test - exact match including PolyA and TSS where present.
            transcripts_assigned = self._assign_path_to_transcript(
                splice_graph,
                mp,
                gene_isoforms,
                test_type="exact",
                fraction_read_align_overlap=fraction_read_align_overlap,
                trim_TSS_polyA=False,
                anchor_PolyA_TSS=True,
            )

            if transcripts_assigned is None:
                # keep TSS,PolyA allow inexact but compatible and read alignment coverage check.
                transcripts_assigned = self._assign_path_to_transcript(
                    splice_graph,
                    mp,
                    gene_isoforms,
                    test_type="COMPATIBLE_CONTAINED",
                    fraction_read_align_overlap=fraction_read_align_overlap,
                    trim_TSS_polyA=False,
                    anchor_PolyA_TSS=True,
                )

            if transcripts_assigned is None:
                # keep TSS,PolyA allow inexact but compatible and read alignment coverage check.
                transcripts_assigned = self._assign_path_to_transcript(
                    splice_graph,
                    mp,
                    gene_isoforms,
                    test_type="INTRONS_CONTAINED",
                    fraction_read_align_overlap=fraction_read_align_overlap,
                    trim_TSS_polyA=False,
                    anchor_PolyA_TSS=True,
                )

            if transcripts_assigned is None:
                # keep TSS,PolyA allow inexact but compatible and read alignment coverage check.
                transcripts_assigned = self._assign_path_to_transcript(
                    splice_graph,
                    mp,
                    gene_isoforms,
                    test_type="other",
                    fraction_read_align_overlap=fraction_read_align_overlap,
                    trim_TSS_polyA=False,
                    anchor_PolyA_TSS=True,
                )
            if transcripts_assigned is None:
                # keep TSS and PolyA features but disable required anchoring of TSS/PolyA, allow inexact but compatible
                transcripts_assigned = self._assign_path_to_transcript(
                    splice_graph,
                    mp,
                    gene_isoforms,
                    test_type="other",
                    fraction_read_align_overlap=fraction_read_align_overlap,
                    trim_TSS_polyA=False,
                    anchor_PolyA_TSS=False,
                )

            ##############################
            ## With TSS and PolyA trimming
            ##############################

            if transcripts_assigned is None:
                # TSS and polyA trimmed, and exact splice path matching required
                transcripts_assigned = self._assign_path_to_transcript(
                    splice_graph,
                    mp,
                    gene_isoforms,
                    test_type="exact",
                    fraction_read_align_overlap=fraction_read_align_overlap,
                    trim_TSS_polyA=True,
                    anchor_PolyA_TSS=False,
                )

            if transcripts_assigned is None:
                # keep TSS,PolyA allow inexact but compatible and read alignment coverage check.
                transcripts_assigned = self._assign_path_to_transcript(
                    splice_graph,
                    mp,
                    gene_isoforms,
                    test_type="COMPATIBLE_CONTAINED",
                    fraction_read_align_overlap=fraction_read_align_overlap,
                    trim_TSS_polyA=True,
                    anchor_PolyA_TSS=False,
                )

            if transcripts_assigned is None:
                # keep TSS,PolyA allow inexact but compatible and read alignment coverage check.
                transcripts_assigned = self._assign_path_to_transcript(
                    splice_graph,
                    mp,
                    gene_isoforms,
                    test_type="INTRONS_CONTAINED",
                    fraction_read_align_overlap=fraction_read_align_overlap,
                    trim_TSS_polyA=True,
                    anchor_PolyA_TSS=False,
                )

            if transcripts_assigned is None:
                # compatibile allowed with read alignment coverage check
                transcripts_assigned = self._assign_path_to_transcript(
                    splice_graph,
                    mp,
                    gene_isoforms,
                    test_type="other",
                    fraction_read_align_overlap=fraction_read_align_overlap,
                    trim_TSS_polyA=True,
                    anchor_PolyA_TSS=False,
                )

            if (
                transcripts_assigned is None
                and LRAA_Globals.config["aggressively_assign_reads"]
            ):
                # last resort, do majority voting
                transcripts_assigned = (
                    self._assign_path_to_transcript_by_majority_voting(
                        splice_graph, mp, gene_isoforms
                    )
                )

            if transcripts_assigned is None:
                logger.debug(
                    "mp_count_pair {} maps to gene but no isoform(transcript)".format(
                        mp_count_pair
                    )
                )
                self._unassigned_mp_count_pairs.append(mp_count_pair)
            else:
                logger.debug(
                    "mp_count_pair {} maps to transcripts:\n{}".format(
                        mp_count_pair,
                        "\n".join(
                            [
                                "{}\t{}".format(x.get_transcript_id(), x._simplepath)
                                for x in transcripts_assigned
                            ]
                        ),
                    )
                )

                self._mp_to_transcripts[mp] = transcripts_assigned

                transcript_read_weights = (
                    self._assign_read_weights_based_on_read_end_agreement(
                        splice_graph, mp, transcripts_assigned
                    )
                )

                for transcript in transcripts_assigned:
                    transcript_id = transcript.get_transcript_id()
                    mp_read_weight = transcript_read_weights[transcript_id]
                    # read_names = mp.get_read_names()
                    logger.debug(
                        "Assigning {} mp: {} read weights as: {}".format(
                            transcript.get_transcript_id(),
                            mp.toShortDescr(),
                            mp_read_weight,
                        )
                    )

                    # assign mp and weight to transcript
                    transcript.add_multipaths_evidence_assigned(mp)
                    transcript.set_multipaths_evidence_weights({mp: mp_read_weight})

                num_paths_assigned += 1
                num_read_counts_assigned += count

        if num_paths_total == 0:
            num_paths_total = 1e-5  # make nonzero to avoid div-by-zero below

        if num_read_counts_total == 0:
            num_read_counts_total = 1e-5  # ditto above

        ## audit summary
        lines = [
            "num_paths_total: {}, num_read_counts_total: {}".format(
                num_paths_total, num_read_counts_total
            ),
            "\tnum_paths_anchored_to_gene: {} = {:.2f}%, num_read_counts_anchored_to_gene: {} = {:.2f}%\n".format(
                num_paths_anchored_to_gene,
                num_paths_anchored_to_gene / num_paths_total * 100,
                num_read_counts_anchored_to_gene,
                num_read_counts_anchored_to_gene / num_read_counts_total * 100,
            ),
            "\tnum_paths_assigned_to_trans: {} = {:.2f}%, num_read_counts_assigned_to_trans: {} = {:.2f}%\n".format(
                num_paths_assigned,
                num_paths_assigned / num_paths_total * 100,
                num_read_counts_assigned,
                num_read_counts_assigned / num_read_counts_total * 100,
            ),
            "\tnum_paths_unanchored: {} = {:.2f}%, num_read_counts_unanchored: {} = {:.2f}%\n".format(
                num_paths_unanchored,
                num_paths_unanchored / num_paths_total * 100,
                num_read_counts_unanchored,
                num_read_counts_unanchored / num_read_counts_total * 100,
            ),
        ]

        try:
            prefix = f"[{ca}{cs}] " if ca is not None and cs is not None else ""
        except Exception:
            prefix = ""

        audit_txt = "\n".join([prefix + x for x in lines])

        logger.debug(audit_txt)
        logger.info(audit_txt)

        if num_read_counts_unanchored > 0:
            logger.warning(
                "%s%d read path(s) carrying %s reads (%.1f%% of reads here) matched no gene "
                "in the splice graph; those reads are excluded from quantification because "
                "no isoform exists at those loci",
                prefix,
                num_paths_unanchored,
                num_read_counts_unanchored,
                num_read_counts_unanchored / num_read_counts_total * 100,
            )

        # finish the progress display cleanly
        if pbar is not None:
            try:
                pbar.close()
            except Exception:
                pass
        elif show_progress and num_mp_count_pairs > 0:
            try:
                sys.stderr.write("\n")
                sys.stderr.flush()
            except Exception:
                pass

        if local_debug is True:
            LRAA_Globals.DEBUG = LRAA_orig_setting
            logging.getLogger().setLevel(logging_orig_setting)

        return

    def _build_read_sharing_gene_components(self, gene_to_transcripts):
        """Group gene_ids into components of genes that share at least one read.

        The EM apportions each read across the transcripts it is compatible with,
        so genes whose transcripts compete for the same read have to be solved in
        one joint EM.  Quantifying them one gene at a time hands the whole read to
        each gene independently, which a deleted post-hoc pass used to try to
        repair after the fact.  Overlapping annotation is enough to produce this,
        and is what does produce it in practice: an ordinary read compatible with
        one transcript of a readthrough gene model and one transcript of its
        constituent gene, e.g. PEDS1-UBE2V1 against UBE2V1.  No unusual alignment
        is involved -- on a chr20 measurement of that population the median
        aligned span was 1,325 bp with a median longest intron of zero.

        Genes are joined when one multipath is compatible with transcripts in both
        of them, read off the compatibility _assign_reads_to_transcripts already
        produced.  There is no second pass over the BAM and no recomputed
        compatibility, and genes sharing no read stay singletons, which is the
        overwhelming majority.

        Returns a list of lists of gene_ids.  Components are ordered by leftmost
        transcript coordinate then gene_id, and genes within a component likewise,
        so neither the EM order nor the transcript order within an EM depends on
        dict or set iteration order.
        """

        gene_id_to_lend = dict()
        for gene_id, transcripts_of_gene in gene_to_transcripts.items():
            gene_id_to_lend[gene_id] = min(
                transcript.get_coords()[0] for transcript in transcripts_of_gene
            )

        # union-find over gene_ids
        gene_id_to_parent = dict((gene_id, gene_id) for gene_id in gene_to_transcripts)

        def _find(gene_id):
            root = gene_id
            while gene_id_to_parent[root] != root:
                root = gene_id_to_parent[root]
            while gene_id_to_parent[gene_id] != root:
                gene_id_to_parent[gene_id], gene_id = root, gene_id_to_parent[gene_id]
            return root

        def _union(gene_id_a, gene_id_b):
            root_a, root_b = _find(gene_id_a), _find(gene_id_b)
            if root_a != root_b:
                gene_id_to_parent[root_b] = root_a

        for transcripts_assigned in self._mp_to_transcripts.values():
            anchor_gene_id = None
            for transcript in transcripts_assigned:
                gene_id = transcript.get_gene_id()
                if gene_id not in gene_id_to_parent:
                    # transcript held over from an earlier quantify() call on this
                    # object; not part of this quantification.
                    continue
                if anchor_gene_id is None:
                    anchor_gene_id = gene_id
                else:
                    _union(anchor_gene_id, gene_id)

        root_to_gene_ids = defaultdict(list)
        for gene_id in gene_id_to_parent:
            root_to_gene_ids[_find(gene_id)].append(gene_id)

        gene_components = list()
        for component_gene_ids in root_to_gene_ids.values():
            component_gene_ids.sort(
                key=lambda gene_id: (gene_id_to_lend[gene_id], gene_id)
            )
            gene_components.append(component_gene_ids)

        # leftmost gene of each component is now its first, so this orders
        # components by leftmost coordinate then gene_id.
        gene_components.sort(
            key=lambda component_gene_ids: (
                gene_id_to_lend[component_gene_ids[0]],
                component_gene_ids[0],
            )
        )

        return gene_components

    def _get_gene_with_best_node_matches_to_simplepath(self, simplepath):

        gene_ranker = defaultdict(int)

        for node in simplepath:
            if node != SPACER:
                if node in self._path_node_id_to_gene_ids:
                    gene_set = self._path_node_id_to_gene_ids[node]
                    for gene_id in gene_set:
                        gene_ranker[gene_id] += 1

        if len(gene_ranker) == 0:
            return None
        else:
            # (-count, gene_id) is a total order.  Ranking on the bare count was
            # a stable sort over dict insertion order, and that insertion order
            # comes from iterating _path_node_id_to_gene_ids[node], a set of
            # gene-id STRINGS whose iteration order is randomised per process by
            # PYTHONHASHSEED.  Equally-matching genes could therefore swap the
            # top slot between runs, and this function picks genes_ranked[0].
            genes_ranked = sorted(
                gene_ranker.keys(), key=lambda x: (-gene_ranker[x], x)
            )
            return genes_ranked[0]

    def _get_all_genes_with_node_matches_to_simplepath(self, simplepath):

        gene_ranker = defaultdict(int)

        for node in simplepath:
            if node != SPACER:
                if node in self._path_node_id_to_gene_ids:
                    gene_set = self._path_node_id_to_gene_ids[node]
                    for gene_id in gene_set:
                        gene_ranker[gene_id] += 1

        if len(gene_ranker) == 0:
            return None
        else:
            genes_ranked = sorted(
                gene_ranker.keys(), key=lambda x: (-gene_ranker[x], x)
            )
            return genes_ranked

    def _assign_path_to_transcript(
        self,
        splice_graph,
        mp,
        transcripts,
        test_type: str,  # choices: ["exact", "FSM", "other"]
        fraction_read_align_overlap,
        trim_TSS_polyA: bool,
        anchor_PolyA_TSS: bool,
    ):

        assert type(mp) == MultiPath.MultiPath
        # A sequence, not a set: the caller's order is now meaningful and must
        # be preserved.  Sorted defensively so the compatibility list this
        # builds is a function of transcript structure alone, whatever order a
        # future caller supplies.
        transcripts = sorted(transcripts, key=Transcript.Transcript.structural_sort_key)
        assert transcripts, "Error, no transcripts supplied"
        assert type(transcripts[0]) == Transcript.Transcript
        assert (
            fraction_read_align_overlap >= 0 and fraction_read_align_overlap <= 1.0
        ), "Error, fraction_read_align_overlap must be between 0 and 1.0"

        test_type_choices = ["exact", "FSM", "COMPATIBLE_CONTAINED", "INTRONS_CONTAINED", "other"]

        assert (
            test_type in test_type_choices
        ), "Error, not recognizing test_type {}".format(test_type)

        contig_strand = splice_graph.get_contig_strand()

        read_sp = mp.get_simple_path()
        mp_id = mp.get_id()

        if trim_TSS_polyA:
            read_sp, read_TSS_id, read_polyA_id = SPU.trim_TSS_and_PolyA(
                read_sp, contig_strand
            )

        # store read name to mp for later debugging.
        for read_name in mp.get_read_names():
            self._read_name_to_multipath[read_name] = mp

        def is_PolyA_or_TSS(simple_node):
            if re.match("^(TSS|POLYA):", simple_node):
                return True
            else:
                return False

        # Precompute read ordered introns and read span (lend, rend) for containment checks
        read_introns_ordered = [node for node in read_sp if re.match("^I:", node)]
        _trimmed_read_no_spacers = SPU.trim_terminal_spacers(read_sp.copy())
        read_span = (
            splice_graph.get_node_obj_via_id(_trimmed_read_no_spacers[0]).get_coords()[0],
            splice_graph.get_node_obj_via_id(_trimmed_read_no_spacers[-1]).get_coords()[1],
        )

        transcripts_compatible_with_read = list()

        logger.debug("** Assessing transcript compatibility for: {}".format(mp))

        for i, transcript in enumerate(transcripts):
            transcript_sp = transcript._simplepath

            transcript_id = transcript.get_transcript_id()
            mp_descr = mp.toShortDescr()

            logger.debug(
                "* evaluating transcript {} compatibility with mp: {}".format(
                    transcript_id, mp_descr
                )
            )

            assert transcript_sp is not None

            logger.debug(
                "[{} trim_TSS_polyA={} test_type={} anchor_PolyA_TSS={}] -evaluating [{}/{}] transcript: {} {}".format(
                    mp_descr,
                    trim_TSS_polyA,
                    test_type,
                    anchor_PolyA_TSS,
                    i + 1,
                    len(transcripts),
                    transcript.get_transcript_id(),
                    transcript_sp,
                )
            )

            if trim_TSS_polyA:
                transcript_sp, transcript_TSS_id, transcript_polyA_id = (
                    SPU.trim_TSS_and_PolyA(transcript_sp, contig_strand)
                )

            else:
                if anchor_PolyA_TSS:

                    #######################################################################################
                    ## Testing first and last positions of read and transcript, position matching required.
                    #######################################################################################

                    fail_msg = None

                    ##########################
                    ## testing first positions
                    ##########################

                    # read first position is TSS or PolyA but transcript lacks it.
                    if is_PolyA_or_TSS(read_sp[0]) and not is_PolyA_or_TSS(
                        transcript_sp[0]
                    ):
                        fail_msg = "read TSS or polyA pos[0] of {} inconsistent w/ transcript {}".format(
                            mp_descr, transcript_id
                        )

                    # both read and isoform have TSS or PolyA at first position, but they don't match.
                    elif (
                        is_PolyA_or_TSS(transcript_sp[0])
                        and is_PolyA_or_TSS(read_sp[0])
                        and transcript_sp[0] != read_sp[0]
                    ):
                        fail_msg = "read TSS or polyA pos[0] of {} inconsistent w/ transcript {}".format(
                            mp_descr, transcript_id
                        )

                    ############################################
                    ## test last position of read and transcript
                    ############################################

                    # read last position is TSS or PolyA, but transcript is not.
                    elif is_PolyA_or_TSS(read_sp[-1]) and not is_PolyA_or_TSS(
                        transcript_sp[-1]
                    ):
                        fail_msg = "read TSS or polyA pos[-1] of {} inconsistent w/ transcript {}".format(
                            mp_descr, transcript_id
                        )

                    # both read and isoform last position is TSS or PolyA, but don't match up.
                    elif (
                        is_PolyA_or_TSS(read_sp[-1])
                        and is_PolyA_or_TSS(transcript_sp[-1])
                        and transcript_sp[-1] != read_sp[-1]
                    ):
                        fail_msg = "read TSS or polyA pos[-1] of {} incosistent w/ transcript {}".format(
                            mp_descr, transcript_id
                        )

                    if fail_msg is not None:
                        logger.debug(
                            "[{} trim_TSS_polyA={} test_type={} anchor_PolyA_TSS={}] -evaluating [{}/{}] transcript: {} {}, FAIL MSG: {}".format(
                                mp_descr,
                                trim_TSS_polyA,
                                test_type,
                                anchor_PolyA_TSS,
                                i + 1,
                                len(transcripts),
                                transcript.get_transcript_id(),
                                transcript_sp,
                                fail_msg,
                            )
                        )
                        continue

            if test_type == "exact":

                ##################################################
                ## Test for exact match of splice paths end-to-end
                ##################################################

                if transcript_sp == read_sp:
                    logger.debug(
                        "{} [trim_TSS_polyA={} test_type={} anchor_PolyA_TSS={}]  Read {} IDENTICAL with transcript {}".format(
                            mp_descr,
                            trim_TSS_polyA,
                            test_type,
                            anchor_PolyA_TSS,
                            read_sp,
                            transcript_sp,
                        )
                    )
                    transcripts_compatible_with_read.append(transcript)
                else:
                    logger.debug(
                        "{} [trim_TSS_polyA={} test_type={} anchor_PolyA_TSS={}]  Read {} NOT_identical with transcript {}".format(
                            mp_descr,
                            trim_TSS_polyA,
                            test_type,
                            anchor_PolyA_TSS,
                            read_sp,
                            transcript_sp,
                        )
                    )

            elif test_type == "FSM":

                if (
                    SPU.simple_paths_have_identical_intron_representation(
                        read_sp, transcript_sp
                    )
                    and SPU.fraction_read_overlap(splice_graph, read_sp, transcript_sp)
                    >= fraction_read_align_overlap
                ):

                    logger.debug(
                        "{} [trim_TSS_polyA={} test_type={} anchor_PolyA_TSS={}]  Read {} FSM with transcript {}".format(
                            mp_descr,
                            trim_TSS_polyA,
                            test_type,
                            anchor_PolyA_TSS,
                            read_sp,
                            transcript_sp,
                        )
                    )
                    # print("Read {} compatible with transcript {}".format(read_sp, transcript_sp))
                    transcripts_compatible_with_read.append(transcript)

            elif test_type == "COMPATIBLE_CONTAINED":

                # check for transcript path fully containing the splice pattern of the read
                if (SPU.path_A_contains_path_B(transcript_sp, read_sp)):
                    logger.debug(
                        "{} [trim_TSS_polyA={} test_type={} anchor_PolyA_TSS={}]  Read {} COMPATIBLE and CONTAINED by transcript {}".format(
                            mp_descr,
                            trim_TSS_polyA,
                            test_type,
                            anchor_PolyA_TSS,
                            read_sp,
                            transcript_sp,
                        )
                    )
                    # print("Read {} compatible with transcript {}".format(read_sp, transcript_sp))
                    transcripts_compatible_with_read.append(transcript)

            elif test_type == "INTRONS_CONTAINED":

                # ordered intron match within read span (no extra transcript introns within span) + meets fraction read overlap criteria
                if (
                    SPU.read_introns_match_transcript_introns_overlapping_read_span(
                        splice_graph, read_sp, transcript_sp
                    )
                    and SPU.fraction_read_overlap(
                        splice_graph, read_sp, transcript_sp
                    )
                    >= fraction_read_align_overlap
                ):

                    logger.debug(
                        "{} [trim_TSS_polyA={} test_type={} anchor_PolyA_TSS={}]  Read {} COMPATIBLE (ordered introns within read span + overlap) with transcript {}".format(
                            mp_descr,
                            trim_TSS_polyA,
                            test_type,
                            anchor_PolyA_TSS,
                            read_sp,
                            transcript_sp,
                        )
                    )
                    transcripts_compatible_with_read.append(transcript)

            elif test_type == "other":

                ######################################################################################
                # Test for compatibilty, no gaps within region of overlap, and sufficient read overlap
                ######################################################################################

                # compatible and meets fraction read overlap criteria
                if (
                    SPU.are_overlapping_and_compatible_NO_gaps_in_overlap(
                        transcript_sp, read_sp
                    )
                    and SPU.fraction_read_overlap(splice_graph, read_sp, transcript_sp)
                    >= fraction_read_align_overlap
                ):

                    logger.debug(
                        "{} [trim_TSS_polyA={} test_type={} anchor_PolyA_TSS={}]  Read {} COMPATIBLE with transcript {}".format(
                            mp_descr,
                            trim_TSS_polyA,
                            test_type,
                            anchor_PolyA_TSS,
                            read_sp,
                            transcript_sp,
                        )
                    )
                    # print("Read {} compatible with transcript {}".format(read_sp, transcript_sp))
                    transcripts_compatible_with_read.append(transcript)

                else:
                    logger.debug(
                        "{} [trim_TSS_polyA={} test_type={} anchor_PolyA_TSS={}]  Read {} NOT_compatible with transcript {}".format(
                            mp_descr,
                            trim_TSS_polyA,
                            test_type,
                            anchor_PolyA_TSS,
                            read_sp,
                            transcript_sp,
                        )
                    )
            else:
                RuntimeError("not recognizing test type: " + test_type)

        if len(transcripts_compatible_with_read) == 0:
            logger.debug(
                "{} NO TRANSCRIPTS FOUND COMPATIBLE WITH READ.".format(mp_descr)
            )
            return None
        else:
            logger.debug(
                "{} FOUND COMPATIBLE WITH\n{}".format(
                    mp_descr,
                    "\n".join([str(x) for x in transcripts_compatible_with_read]),
                )
            )
            return transcripts_compatible_with_read

    def _assign_path_to_transcript_by_majority_voting(
        self, splice_graph, mp, transcripts
    ):

        assert type(mp) == MultiPath.MultiPath
        transcripts = sorted(transcripts, key=Transcript.Transcript.structural_sort_key)
        assert transcripts, "Error, no transcripts supplied"
        assert type(transcripts[0]) == Transcript.Transcript

        contig_strand = splice_graph.get_contig_strand()

        read_sp = mp.get_simple_path()

        ## For majority voting, let's trim TSS and polyA so it doesn't contribute to the scoring.
        # read_sp, read_TSS_id, read_polyA_id = SPU.trim_TSS_and_PolyA(
        #    read_sp, contig_strand
        # )

        # store read name to mp for later debugging.
        for read_name in mp.get_read_names():
            self._read_name_to_multipath[read_name] = mp

        scored_transcripts = list()

        for transcript in transcripts:
            transcript_sp = transcript._simplepath

            shared_simple_nodes = [
                simple_node for simple_node in read_sp if simple_node in transcript_sp
            ]

            if len(shared_simple_nodes) > 0:
                scored_transcripts.append([len(shared_simple_nodes), transcript])

        if len(scored_transcripts) > 0:
            scored_transcripts = sorted(
                scored_transcripts, key=lambda x: x[0], reverse=True
            )
            logger.debug(
                "Majority Voting: Candidate order for {} is:\n{} ".format(
                    mp.toShortDescr(),
                    "\n".join(
                        [
                            "score:{}\t{}".format(str(x[0]), str(x[1]))
                            for x in scored_transcripts
                        ]
                    ),
                )
            )

            top_transcript_score_pair = scored_transcripts.pop(0)
            top_transcript_score, top_transcript = top_transcript_score_pair
            top_transcripts = [top_transcript]
            # capture ties
            for alt_top_transcript in scored_transcripts:
                if alt_top_transcript[0] == top_transcript_score:
                    top_transcripts.append(alt_top_transcript[1])

            logger.debug(
                "Majority Voting CHOOSING TOP CANDIDATE for {} as: {}".format(
                    mp.toShortDescr(), str(top_transcripts)
                )
            )
            return top_transcripts

        else:
            return None

    def _assign_read_weights_based_on_read_end_agreement(
        self, splice_graph, mp, transcripts_assigned
    ):

        mp_sp = mp.get_simple_path()

        # Determine which end corresponds to the 3' end given strand.
        # Simple paths are ordered by genomic coordinate; 3' is 'rend' on '+' and 'lend' on '-'.
        try:
            contig_strand = splice_graph.get_contig_strand()
        except Exception:
            contig_strand = "+"
        three_prime_end_key = "rend" if contig_strand == "+" else "lend"

        # Weight only by agreement of the read 3' end (rend) with the transcript 3' end.
        # The previous behavior used both 5' (lend) and 3' (rend) distances; now we use only 3' end.
        transcript_id_to_sum_end_dists = dict()
        sum_dists = 0
        for transcript in transcripts_assigned:
            transcript_id = transcript.get_transcript_id()
            transcript_sp = transcript.get_simple_path()
            dist_rend = self._get_simple_path_dist_to_termini(
                splice_graph, mp_sp, transcript_sp, three_prime_end_key
            )
            """
            logger.debug(
                "MP {} lend dist for {} = {}".format(
                    mp.get_id(), transcript_id, dist_lend
                )
            )
            logger.debug(
                "MP {} rend dist for {} = {}".format(
                    mp.get_id(), transcript_id, dist_rend
                )
            )
            """
            # Only consider the 3' end distance for weighting
            sum_dist = dist_rend
            transcript_id_to_sum_end_dists[transcript_id] = sum_dist
            sum_dists += sum_dist

        # determine relative weightings
        # Pseudocount sized to the alt-PolyA aggregation window: sub-window distance
        # differences should not drive any transcript's weight to zero.
        transcript_id_to_mp_weights = dict()
        sum_weights = 0.0
        num_transcripts_assigned = len(transcripts_assigned)
        end_dist_pseudocount = LRAA_Globals.config["max_dist_between_alt_polyA_sites"]
        adj_sum_dists = sum_dists + num_transcripts_assigned * end_dist_pseudocount
        for transcript in transcripts_assigned:
            transcript_id = transcript.get_transcript_id()
            dist = transcript_id_to_sum_end_dists[transcript_id]
            logger.debug(
                "transcript {} has sum_end_dist: {} and total_sum_dists_all_trans: {}".format(
                    transcript_id, dist, sum_dists
                )
            )
            weight = (
                1.0 - (dist + end_dist_pseudocount) / adj_sum_dists
                if adj_sum_dists > 0
                else 1.0
            )
            transcript_id_to_mp_weights[transcript_id] = weight
            sum_weights += weight

        # renormalize weights
        for transcript, weight in transcript_id_to_mp_weights.items():
            transcript_id_to_mp_weights[transcript] = (
                weight / sum_weights
                if sum_weights > 0
                else 1 / num_transcripts_assigned
            )

        return transcript_id_to_mp_weights

    def _get_simple_path_dist_to_termini(
        self, splice_graph, source_sp, target_sp, from_which_end
    ):

        assert from_which_end in (
            "lend",
            "rend",
        ), "from_which_end must be 'lend' or 'rend'"

        if from_which_end == "rend":
            # reverse it so we can walk from left instead of from right.
            target_sp = list(reversed(target_sp))

        matching_node_idx = None
        for idx, node_id in enumerate(target_sp):
            if node_id in source_sp:
                matching_node_idx = idx
                break

        assert (
            matching_node_idx is not None
        ), "Error, cannot find matching node in target path {} from source path {}".format(
            target_sp, source_sp
        )

        dist = 0
        for node_id in target_sp[0:matching_node_idx]:
            node = splice_graph.get_node_obj_via_id(node_id)
            if type(node) == Exon:
                dist += node.get_feature_length()

        return dist

    def _estimate_isoform_read_support(self, transcripts, prefix_str=None):
        """

        Given the reads assigned to the transcript (accessed with transcript.get_read_names() )
            Read counts are assigned to transcripts taking into account multiple read mappings
            Without EM (run_EM is False), read counts are equally divided among the isoforms they're assigned.
            With EM, they're assigned fractionally according to inferred isoform expression levels in EM cycles.

            Final read counts and isoform fraction values are stored in the transcript objects themselves, and accessed as:
                transcript.get_read_counts_assigned() and transcript.get_isoform_fraction()

        returns transcript_to_fractional_read_assignment
             with structure [transcript_id][read_name] = frac_read_assigned

        """

        try:
            if prefix_str:
                logger.info(
                    f"{prefix_str}-estimating isoform read support for {len(transcripts)} transcripts."
                )
            else:
                logger.info(
                    "-estimating isoform read support for {} transcripts.".format(
                        len(transcripts)
                    )
                )
        except Exception:
            logger.info(
                "-estimating isoform read support for {} transcripts.".format(
                    len(transcripts)
                )
            )

        start_time = time.time()

        transcript_to_expr_val = defaultdict(float)
        transcript_to_fractional_mp_assignment = defaultdict(dict)
        transcript_to_read_count = defaultdict(float)

        if self._run_EM:

            (
                transcript_to_expr_val,
                transcript_to_fractional_mp_assignment,
                transcript_to_read_count,
            ) = EM.run_EM(transcripts, self._max_EM_iterations, prefix_str=prefix_str)

        else:
            # simple equal fractional assignment of reads to compatible transcripts

            # first populate read_name to list of transcripts compatible with read
            all_mps = set()
            mp_to_transcripts = defaultdict(set)
            for transcript in transcripts:
                mps = transcript.get_multipaths_evidence_assigned()
                for mp in mps:
                    mp_to_transcripts[mp].add((transcript))
                    all_mps.add(mp)

            # get total read count
            total_mapped_reads = 0
            for mp in all_mps:
                total_mapped_reads += mp.get_read_count()

            for transcript in transcripts:
                transcript_read_count_total = 0
                transcript_id = transcript.get_transcript_id()

                mps = transcript.get_multipaths_evidence_assigned()
                for mp in mps:

                    frac_mp_assignment = 0
                    all_transcripts_with_mp = mp_to_transcripts[mp]

                    num_transcripts_with_mp = len(all_transcripts_with_mp)

                    # split read equally across all copatible reads.

                    frac_mp_assignment = (
                        1 / num_transcripts_with_mp
                        if num_transcripts_with_mp > 0
                        else 0
                    )

                    num_reads_in_mp = mp.get_read_count()

                    transcript_read_count_total += frac_mp_assignment * num_reads_in_mp
                    transcript_to_fractional_mp_assignment[transcript_id][
                        mp
                    ] = frac_mp_assignment

                transcript_to_read_count[transcript_id] = transcript_read_count_total
                transcript_to_expr_val[transcript_id] = (
                    transcript_read_count_total / total_mapped_reads
                    if total_mapped_reads > 0
                    else 0
                )
                logger.debug(
                    f"-assigning transcript {transcript_id} read count: {transcript_read_count_total} and expr val {transcript_read_count_total}/{total_mapped_reads} = {transcript_to_expr_val[transcript_id]}"
                )

        ## assign final read counts to each transcript object.
        # The transcripts handed here are a read-sharing component of genes, which
        # is the EM unit but is NOT the reporting unit.  isoform_fraction is a
        # transcript's share of its OWN gene, so the regrouping by gene_id below is
        # load-bearing rather than incidental bookkeeping: it is what keeps the
        # fraction gene-scoped while the EM that produced the counts was not.
        gene_to_transcripts = defaultdict(list)
        for transcript in transcripts:
            transcript_id = transcript.get_transcript_id()
            transcript_read_count = transcript_to_read_count[transcript_id]
            transcript.set_read_counts_assigned(transcript_read_count)
            gene_id = transcript.get_gene_id()
            gene_to_transcripts[gene_id].append(transcript)

        ## isoform isoform fraction
        for gene_id in gene_to_transcripts:
            transcripts_of_gene = gene_to_transcripts[gene_id]

            # evaluate isoform fraction
            sum_gene_reads = 0
            for transcript_of_gene in transcripts_of_gene:
                transcript_read_count = transcript_of_gene.get_read_counts_assigned()
                sum_gene_reads += transcript_read_count

            logger.debug(
                "gene_id {} has total reads: {}".format(gene_id, sum_gene_reads)
            )

            for transcript_of_gene in transcripts_of_gene:
                transcript_id = transcript_of_gene.get_transcript_id()
                transcript_read_count = transcript_of_gene.get_read_counts_assigned()
                isoform_frac = (
                    transcript_read_count / sum_gene_reads if sum_gene_reads > 0 else 0
                )
                logger.debug(
                    "\ttranscript_id {} has {} reads = {} isoform fraction of {}".format(
                        transcript_id, transcript_read_count, isoform_frac, gene_id
                    )
                )
                transcript_of_gene.set_isoform_fraction(isoform_frac)

        end_time = time.time()

        try:
            if prefix_str:
                logger.info(
                    f"{prefix_str}Time for quant of {len(transcripts)} transcripts = {(end_time - start_time) / 60:.2f} minutes"
                )
            else:
                logger.info(
                    "Time for quant of {} transcripts = {:.2f} minutes".format(
                        len(transcripts), (end_time - start_time) / 60
                    )
                )
        except Exception:
            logger.info(
                "Time for quant of {} transcripts = {:.2f} minutes".format(
                    len(transcripts), (end_time - start_time) / 60
                )
            )

        return transcript_to_fractional_mp_assignment

    def get_mp_to_transcripts(self):
        return self._mp_to_transcripts

    def get_unassigned_mp_count_pairs(self):
        return list(self._unassigned_mp_count_pairs)

    def _assert_component_identity_valid(self):
        if not self._component_identity_valid:
            raise RuntimeError(
                "component identity is not available: quantify() has not completed "
                "on this object. Returning an empty mapping instead would be read as "
                "'no components', and a consumer renormalizing per gene on that "
                "basis inflates any read compatible with two genes."
            )

    def get_transcript_id_to_component_id(self):
        """transcript_id -> opaque id of the read-sharing component holding it.

        The EM unit is a component of genes joined by shared reads, so theta is
        normalized over a component and not over a gene.  A consumer computing a
        per-read isoform split must renormalize over a component's transcripts;
        renormalizing over a gene's gives a read spanning two genes of one
        component a split summing to 2.

        Every quantified transcript appears, singletons included.  Rebuilt by each
        quantify() call, so it never carries an entry from a previous call.

        Raises before the first call and after a failed one, rather than returning
        an empty mapping a caller could mistake for a real answer.
        """
        self._assert_component_identity_valid()
        return dict(self._transcript_id_to_component_id)

    def get_component_id_to_gene_ids(self):
        """component_id -> the gene_ids it holds, ordered leftmost gene first."""
        self._assert_component_identity_valid()
        return dict(
            (component_id, list(gene_ids))
            for component_id, gene_ids in self._component_id_to_gene_ids.items()
        )

    def get_unassigned_read_names(self):
        read_names = set()
        for mp_count_pair in self._unassigned_mp_count_pairs:
            try:
                mp, _ = mp_count_pair.get_multipath_and_count()
                read_names.update(mp.get_read_names())
            except Exception:
                continue
        return read_names

    def dump_mp_to_transcripts_to_file(self, output_filename, contig_acc, strand):

        logger.debug(
            "Writing mp to read and transcript assignment files for {} {}".format(
                contig_acc, strand
            )
        )

        with open(output_filename, "at") as ofh_mp:
            with open(output_filename + ".abridged", "at") as ofh_mp_abridged:

                mp_to_transcripts = self.get_mp_to_transcripts()
                for mp, transcripts in mp_to_transcripts.items():
                    read_names = mp.get_read_names()
                    print(
                        "\t".join(
                            [
                                contig_acc,
                                strand,
                                mp.get_id() + ":" + str(mp.get_simple_path()),
                                str(len(transcripts)),
                                ";".join([x.get_transcript_id() for x in transcripts]),
                                str(len(read_names)),
                            ]
                        ),
                        file=ofh_mp_abridged,
                    )

                    print(
                        "\t".join(
                            [
                                contig_acc,
                                strand,
                                mp.get_id() + ":" + str(mp.get_simple_path()),
                                str(len(transcripts)),
                                ";".join([x.get_transcript_id() for x in transcripts]),
                                str(len(read_names)),
                                ";".join(list(read_names)),
                            ]
                        ),
                        file=ofh_mp,
                    )

        return

    def report_quant_results(
        self,
        transcripts,
        transcript_to_fractional_mp_assignment,
        ofh_quant_vals,
        ofh_read_tracking,
        splice_compatible_containments=None,
        splice_compatible_contained_by=None,
    ):

        ## generate final report.

        ## sort descendingly by read support
        transcripts = sorted(
            transcripts,
            key=lambda x: (x.get_read_counts_assigned(), x.get_transcript_id()),
            reverse=True,
        )

        # first, get sum of reads per gene
        gene_to_read_count = defaultdict(int)
        total_reported_read_count = 0
        for transcript in transcripts:
            gene_id = transcript.get_gene_id()
            counts = transcript.get_read_counts_assigned()
            gene_to_read_count[gene_id] += counts
            total_reported_read_count += counts

        for transcript in transcripts:
            transcript_id = transcript.get_transcript_id()
            num_exons = transcript.get_num_exon_segments()
            transcript_splice_hash_code = (
                Util_funcs.get_hash_code(transcript.get_introns_string())
                if num_exons > 1
                else transcript_id
            )
            gene_id = transcript.get_gene_id()
            output_gene_id = transcript.get_output_gene_id()
            counts = transcript.get_read_counts_assigned()
            isoform_frac = transcript.get_isoform_fraction()

            multipaths = transcript.get_multipaths_evidence_assigned()

            tpm = (
                counts / total_reported_read_count * 1e6
                if total_reported_read_count > 0
                else 0
            )
            rpm_total_reads = transcript.get_TPM()

            num_uniquely_assigned_reads = 0

            for mp in multipaths:
                frac_read_assigned = transcript_to_fractional_mp_assignment[
                    transcript_id
                ][mp]

                mp_id = mp.get_id()

                # get_read_names() returns a set of str, whose iteration order
                # is randomised per process by PYTHONHASHSEED.  These become one
                # output row each, so unsorted it reorders quant.tracking.
                for readname in sorted(mp.get_read_names()):

                    tracking_report_info = [
                        output_gene_id,
                        transcript_id,
                        transcript_splice_hash_code,
                        f"{num_exons}",
                        mp_id,
                        readname,
                        "{:.6f}".format(frac_read_assigned),
                    ]

                    read_weight = (
                        transcript.get_multipath_weight(mp)
                        if LRAA_Globals.config["weight_reads_by_3prime_agreement"]
                        else 1.0
                    )
                    tracking_report_info.append("{:.6f}".format(read_weight))

                    # Always emit read tracking rows for robustness; downstream annotator filters/consumes as needed
                    print("\t".join(tracking_report_info), file=ofh_read_tracking)

                    if (
                        frac_read_assigned
                        >= LRAA_Globals.config["unique_read_report_min_frac"]
                    ):
                        num_uniquely_assigned_reads += 1

            gene_read_count = gene_to_read_count[gene_id]
            unique_gene_read_fraction = (
                num_uniquely_assigned_reads / gene_read_count
                if gene_read_count > 0
                else 0
            )

            report_vals = [
                output_gene_id,
                transcript_id,
                f"{num_uniquely_assigned_reads}",
                f"{counts:.1f}",
                f"{isoform_frac:.3f}",
                f"{unique_gene_read_fraction:0.3f}",
                f"{tpm:.3f}",
                transcript.get_exons_string(),
                (transcript.get_introns_string() if num_exons > 1 else ""),
                transcript_splice_hash_code,
            ]

            if splice_compatible_containments is not None:
                # These are sets of transcript ids; emit them sorted so the report is
                # reproducible between runs rather than following set iteration order.
                def _format_id_set(id_sets, transcript_id):
                    if transcript_id not in id_sets:
                        return ""
                    return "{" + ",".join(sorted(id_sets[transcript_id])) + "}"

                report_vals.append(
                    _format_id_set(splice_compatible_containments, transcript_id)
                )
                report_vals.append(
                    _format_id_set(splice_compatible_contained_by, transcript_id)
                )

            report_vals.append(f"{rpm_total_reads:.3f}")

            report_txt = "\t".join(report_vals)

            logger.debug(report_txt)
            print(report_txt, file=ofh_quant_vals)

            """
            if (DEBUG):
                print("transcript_id\t{}\n{}".format(transcript_id, transcript._simplepath), file=ofh_read_tracking)
                for readname in readnames:
                    print("read:\t{}\n{}".format(readname,
                                                 self._read_name_to_multipath[readname].get_simple_path()),
                          file=ofh_read_tracking)
                print("\n", file=ofh_read_tracking)
            else:
                print("\t".join([gene_id, transcript_id, ",".join(readnames)]), file=ofh_read_tracking)
            """

        return

    @staticmethod
    def get_gene_read_counts(frac_read_assignments, transcript_id_to_transcript_obj):
        gene_id_to_read_count = defaultdict(int)
        for (
            transcript_id,
            transcript_read_frac_assignments,
        ) in frac_read_assignments.items():
            gene_id = transcript_id_to_transcript_obj[transcript_id].get_gene_id()
            for mp, frac_assigned in transcript_read_frac_assignments.items():
                gene_id_to_read_count[gene_id] += frac_assigned * mp.get_read_count()

        return gene_id_to_read_count

    def get_mp_read_count(self, mp):
        return mp.get_read_count()
