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
import lmdb
from GenomeFeature import Exon
import EM
import time

logger = logging.getLogger(__name__)


class Quantify:

    def __init__(self, run_EM, max_EM_iterations, quant_mode="final"):

        self._run_EM = run_EM  # boolean
        self._max_EM_iterations = max_EM_iterations

        self._path_node_id_to_gene_ids = defaultdict(set)

        self._gene_id_to_transcript_objs = defaultdict(set)

        self._read_name_to_multipath = dict()

        self._mp_to_transcripts = dict()

        self._quant_mode = quant_mode

        return

    def quantify(self, splice_graph, transcripts, mp_counter):

        assert type(transcripts) == list
        assert type(transcripts[0]) == Transcript.Transcript
        assert type(mp_counter) == MultiPathCounter.MultiPathCounter

        contig_acc = splice_graph.get_contig_acc()
        contig_strand = splice_graph.get_contig_strand()

        # init transcript quant info
        gene_to_transcripts = defaultdict(list)
        for transcript in transcripts:
            transcript.init_quant_info()
            gene_id = transcript.get_gene_id()
            gene_to_transcripts[gene_id].append(transcript)

        logger.info(
            "have {} genes and {} isoforms to quantify.".format(
                len(gene_to_transcripts), len(transcripts)
            )
        )

        # assign path nodes to gene
        # also assign gene_id to transcript objs
        self._assign_path_nodes_to_gene(transcripts)

        self._assign_reads_to_transcripts(splice_graph, mp_counter)

        transcript_to_fractional_read_assignment = dict()

        for gene_id, transcripts_list in gene_to_transcripts.items():

            trans_coords = list()
            for transcript in transcripts_list:
                trans_coords.extend(transcript.get_coords())

            trans_coords = sorted(trans_coords)
            gene_lend = trans_coords[0]
            gene_rend = trans_coords[-1]
            logger.info(
                f"quant estimates for isoforms of {gene_id} {contig_acc}{contig_strand}:{gene_lend}-{gene_rend}"
            )

            gene_transcript_to_fractional_read_assignment = (
                self._estimate_isoform_read_support(transcripts_list)
            )
            # copy over to the full data structure
            for transcript_id in gene_transcript_to_fractional_read_assignment:
                transcript_to_fractional_read_assignment[transcript_id] = (
                    gene_transcript_to_fractional_read_assignment[transcript_id]
                )

        # see documentation for _estimate_isoform_read_support() below

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
            self._gene_id_to_transcript_objs[gene_id].add(transcript)

            for node_id in simplepath:
                if node_id != SPACER:
                    self._path_node_id_to_gene_ids[node_id].add(gene_id)

        return

    def _assign_reads_to_transcripts(
        self,
        splice_graph,
        mp_counter,
        fraction_read_align_overlap=LRAA_Globals.config["fraction_read_align_overlap"],
    ):

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
        logger.info("- have {} mp_count_pairs".format(num_mp_count_pairs))

        gene_unanchored_mp_count_pairs = list()

        num_paths_total = 0
        num_read_counts_total = 0

        num_paths_anchored_to_gene = 0
        num_read_counts_anchored_to_gene = 0

        num_paths_assigned = 0
        num_read_counts_assigned = 0

        mp_seen = set()

        num_mp_count_pairs_processed = 0

        for mp_count_pair in mp_count_pairs:

            num_mp_count_pairs_processed += 1
            if num_mp_count_pairs_processed % 100 == 0:
                print(
                    "\r{}/{} = {:.3f} mp_count_pairs processed   ".format(
                        num_mp_count_pairs_processed,
                        num_mp_count_pairs,
                        num_mp_count_pairs_processed / num_mp_count_pairs,
                    ),
                    end="",
                )

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
                gene_unanchored_mp_count_pairs.append(mp_count_pair)
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
            gene_isoforms = set()
            for top_gene in top_genes:
                isoforms = self._gene_id_to_transcript_objs[top_gene]
                for isoform in isoforms:
                    gene_isoforms.add(isoform)

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
        audit_txt = "\n".join(
            [
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
            ]
        )

        logger.debug(audit_txt)
        logger.info(audit_txt)

        if local_debug is True:
            LRAA_Globals.DEBUG = LRAA_orig_setting
            logging.getLogger().setLevel(logging_orig_setting)

        return

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
            genes_ranked = sorted(
                gene_ranker.keys(), key=lambda x: gene_ranker[x], reverse=True
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
                gene_ranker.keys(), key=lambda x: gene_ranker[x], reverse=True
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
        assert (
            type(transcripts) == set
        ), "Error, type(transcripts) is {} not set ".format(type(transcripts))
        assert type(list(transcripts)[0]) == Transcript.Transcript
        assert (
            fraction_read_align_overlap >= 0 and fraction_read_align_overlap <= 1.0
        ), "Error, fraction_read_align_overlap must be between 0 and 1.0"

        test_type_choices = ["exact", "FSM", "other"]

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

            else:  # test_type == "other"

                ######################################################################################
                # Test for compatibilty, no gaps within region of overlap, and sufficient read overlap
                ######################################################################################

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
        assert (
            type(transcripts) == set
        ), "Error, type(transcripts) is {} not set ".format(type(transcripts))
        assert type(list(transcripts)[0]) == Transcript.Transcript

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

        transcript_id_to_sum_end_dists = dict()
        sum_dists = 0
        for transcript in transcripts_assigned:
            transcript_id = transcript.get_transcript_id()
            transcript_sp = transcript.get_simple_path()
            dist_lend = self._get_simple_path_dist_to_termini(
                splice_graph, mp_sp, transcript_sp, "lend"
            )
            dist_rend = self._get_simple_path_dist_to_termini(
                splice_graph, mp_sp, transcript_sp, "rend"
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
            sum_dist = dist_lend + dist_rend
            transcript_id_to_sum_end_dists[transcript_id] = sum_dist
            sum_dists += sum_dist

        # determine relative weightings
        transcript_id_to_mp_weights = dict()
        sum_weights = 0.0
        num_transcripts_assigned = len(transcripts_assigned)
        for transcript in transcripts_assigned:
            transcript_id = transcript.get_transcript_id()
            dist = transcript_id_to_sum_end_dists[transcript_id]
            logger.debug(
                "transcript {} has sum_end_dist: {} and total_sum_dists_all_trans: {}".format(
                    transcript_id, dist, sum_dists
                )
            )
            weight = 1.0 - (dist / sum_dists) if sum_dists > 0 else 1.0
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

    def _estimate_isoform_read_support(self, transcripts):
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
            ) = EM.run_EM(transcripts, self._max_EM_iterations)

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

        logger.info(
            "Time for quant of {} transcripts = {:.2f} minutes".format(
                len(transcripts), (end_time - start_time) / 60
            )
        )

        return transcript_to_fractional_mp_assignment

    def get_mp_to_transcripts(self):
        return self._mp_to_transcripts

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
        ofh_quant_read_tracking_lmdb=None,
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
        for transcript in transcripts:
            gene_id = transcript.get_gene_id()
            counts = transcript.get_read_counts_assigned()
            gene_to_read_count[gene_id] += counts

        for transcript in transcripts:
            transcript_id = transcript.get_transcript_id()
            gene_id = transcript.get_gene_id()
            counts = transcript.get_read_counts_assigned()
            isoform_frac = transcript.get_isoform_fraction()

            multipaths = transcript.get_multipaths_evidence_assigned()

            tpm = transcript.get_TPM()

            num_uniquely_assigned_reads = 0

            for mp in multipaths:
                frac_read_assigned = transcript_to_fractional_mp_assignment[
                    transcript_id
                ][mp]

                mp_id = mp.get_id()

                for readname in mp.get_read_names():

                    tracking_report_info = [
                        gene_id,
                        transcript_id,
                        mp_id,
                        readname,
                        "{:.3f}".format(frac_read_assigned),
                    ]

                    if LRAA_Globals.DEBUG:
                        tracking_report_info.append(
                            "{:.3f}".format(transcript.get_multipath_weight(mp))
                        )

                    if frac_read_assigned > 1e-3 or LRAA_Globals.DEBUG:
                        print("\t".join(tracking_report_info), file=ofh_read_tracking)

                    if frac_read_assigned == 1:
                        num_uniquely_assigned_reads += 1

            gene_read_count = gene_to_read_count[gene_id]
            unique_gene_read_fraction = (
                num_uniquely_assigned_reads / gene_read_count
                if gene_read_count > 0
                else 0
            )

            report_vals = [
                gene_id,
                transcript_id,
                f"{num_uniquely_assigned_reads}",
                f"{counts:.1f}",
                f"{isoform_frac:.3f}",
                f"{unique_gene_read_fraction:0.3f}",
                f"{tpm:.3f}",
                transcript.get_exons_string(),
                (
                    transcript.get_introns_string()
                    if transcript.get_num_exon_segments() > 1
                    else ""
                ),
            ]

            if splice_compatible_containments is not None:
                splice_compat_containment_vals = (
                    str(splice_compatible_containments[transcript_id])
                    if transcript_id in splice_compatible_containments
                    else ""
                )
                report_vals.append(splice_compat_containment_vals)

                splice_compat_contained_by_vals = (
                    str(splice_compatible_contained_by[transcript_id])
                    if transcript_id in splice_compatible_contained_by
                    else ""
                )
                report_vals.append(splice_compat_contained_by_vals)

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
