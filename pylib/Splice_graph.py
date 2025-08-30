#!/usr/bin/env python3
# encoding: utf-8

import sys, os, re
import subprocess
import logging
import string
import pysam
from collections import defaultdict
import networkx as nx
import intervaltree as itree
from GenomeFeature import *
from Pretty_alignment_manager import Pretty_alignment_manager
import LRAA_Globals
import statistics

logger = logging.getLogger(__name__)


class Splice_graph:

    # ---------------
    # class variables
    # ---------------

    _read_aln_gap_merge_int = 10  ## internal alignment gap merging
    _inter_exon_segment_merge_dist = 50  ## unbranched exon segments within this range get merged into single segments.
    _max_genomic_contig_length = 1e10

    # noise filtering params
    _min_alt_splice_freq = 0.01
    _min_alt_unspliced_freq = 0.20
    _max_intron_length_for_exon_segment_filtering = 10000
    _min_intron_support = 2
    _min_terminal_splice_exon_anchor_length = 15

    _remove_unspliced_introns = False  # True

    def __init__(self):

        # ------------------
        # instance variables
        # self._alignments_bam_filename = ""

        self._contig_acc = ""
        self._contig_strand = None
        self._contig_seq_str = ""
        self._contig_len = 0

        self._contig_base_cov = list()

        self._splice_graph = None  # becomes networkx digraph TODO://rename this var as confusing when using _splice_graph as this obj in other modules

        self._node_id_to_node = dict()
        self._itree_exon_segments = (
            None  # itree.IntervalTree() # stores Exon, TSS, and PolyAsite objects
        )

        self._intron_objs = dict()  # "lend:rend" => intron_obj
        self._itree_introns = None  # itree.IntervalTree()

        ## connected components (genes)
        self._components = list()  # ordered list of graph components
        self._node_id_to_component = dict()  # node_id -> component index

        self._TSS_objs = list()  # TSS objects
        self._PolyA_objs = list()  # PolyA site objects

        self._input_transcript_lend_boundaries = set()
        self._input_transcript_rend_boundaries = set()
        self._input_transcript_exon_coords_itree = itree.IntervalTree()
        self._input_transcripts_introns = dict()  # "lend:rend" => intron_obj

        return

    @classmethod
    def init_sg_params(cls, params):
        cls._read_aln_gap_merge_int = params[
            "read_aln_gap_merge_int"
        ]  ## internal alignment gap merging
        cls._inter_exon_segment_merge_dist = params[
            "inter_exon_segment_merge_dist"
        ]  ## unbranched exon segments within this range get merged into single segments.
        cls._max_genomic_contig_length = params["max_genomic_contig_length"]

        # noise filtering params
        cls._min_alt_splice_freq = params["min_alt_splice_freq"]
        cls._min_alt_unspliced_freq = params["min_alt_unspliced_freq"]
        cls._max_intron_length_for_exon_segment_filtering = params[
            "max_intron_length_for_exon_segment_filtering"
        ]
        cls._min_intron_support = params["min_intron_support"]
        cls._min_terminal_splice_exon_anchor_length = params[
            "min_terminal_splice_exon_anchor_length"
        ]

        cls._remove_unspliced_introns = params["remove_unspliced_introns"]

    def get_contig_acc(self):
        return self._contig_acc

    def get_contig_strand(self):
        return self._contig_strand

    def set_read_aln_gap_merge(self, read_aln_gap_merge_int):

        self._read_aln_gap_merge_int = read_aln_gap_merge_int

        return

    def get_intron_node_obj(self, intron_lend, intron_rend):

        intron_token = "{}:{}".format(intron_lend, intron_rend)

        if intron_token in self._intron_objs:
            return self._intron_objs[intron_token]

        return None

    def get_node_obj_via_id(self, node_id):

        if node_id not in self._node_id_to_node:
            raise RuntimeError(
                "Error, node_id: {} not recognized in splice graph".format(node_id)
            )

        return self._node_id_to_node[node_id]

    def get_successors(self, node):
        return list(self._splice_graph.successors(node))

    def get_predecessors(self, node):
        return list(self._splice_graph.predecessors(node))

    def get_overlapping_exon_segments(
        self, range_lend, range_rend, min_frac_feature_overlap=0.0
    ):

        overlapping_exon_segments = list()

        for overlapping_interval in self._itree_exon_segments[
            range_lend : range_rend + 1
        ]:

            node = overlapping_interval.data
            lend, rend = node.get_coords()
            feature_len = rend - lend + 1

            if lend <= range_rend and rend >= range_lend:
                # require substantial feature overlap

                if min_frac_feature_overlap > 0:
                    coords = sorted([lend, rend, range_lend, range_rend])
                    overlap_len = coords[2] - coords[1] + 1
                    frac_feature_len_overlap = overlap_len / feature_len

                    target_range_len = range_rend - range_lend + 1
                    frac_target_range_len_overlap = overlap_len / target_range_len

                    logger.debug(
                        "\t\tRange {}-{} overlaps {}, overlap_len = {}, frac_feature_overlap = {:.3f}, frac_range_overlap = {:.3f}".format(
                            range_lend,
                            range_rend,
                            node,
                            overlap_len,
                            frac_feature_len_overlap,
                            frac_target_range_len_overlap,
                        )
                    )

                    if (
                        frac_feature_len_overlap >= min_frac_feature_overlap
                        or frac_target_range_len_overlap >= min_frac_feature_overlap
                    ):  # LRAA_Globals.config['min_feature_frac_overlap']:
                        overlapping_exon_segments.append(node)

                else:
                    overlapping_exon_segments.append(node)

        return overlapping_exon_segments

    def get_overlapping_introns(self, lend, rend):

        overlapping_introns = list()

        for overlapping_intron_interval in self._itree_introns[lend : rend + 1]:
            overlapping_introns.append(overlapping_intron_interval.data)

        return overlapping_introns

    def build_splice_graph_for_contig(
        self,
        contig_acc,
        contig_strand,
        contig_seq_str,
        alignments_bam_file,
        region_lend,
        region_rend,
        input_transcripts=None,
        quant_mode=False,
        restrict_splice_type=None,
        SE_read_encapsulation_mask=None,
    ):

        logger.info(
            "creating splice graph for {} leveraging bam {}, strand {}".format(
                contig_acc, alignments_bam_file, contig_strand
            )
        )

        self._contig_acc = contig_acc
        self._contig_strand = contig_strand
        self._contig_seq_str = contig_seq_str
        self._contig_seq_len = len(contig_seq_str)
        self._region_lend = region_lend
        self._region_rend = region_rend
        self._restrict_splice_type = restrict_splice_type

        ## do the work:

        self._initialize_contig_coverage()  # populates self._contig_base_cov

        ##---------------------------------------------------
        ## intron extraction and exonic base coverage defined.
        # -requires min intron-supporting reads
        # -excludes introns w/ heavily unbalanced splice site support (via Splice_graph._min_alt_splice_freq setting)
        # - stores intron objs in self._intron_objs
        # - base coverage incremented under self._contig_base_cov
        # - defines and stores TSS and PolyA objects

        if alignments_bam_file is not None:
            self._populate_exon_coverage_and_extract_introns(
                alignments_bam_file, contig_acc, contig_strand, quant_mode, restrict_splice_type,
                SE_read_encapsulation_mask
            )

        # incorporate guide structures if provided
        if input_transcripts:
            self._integrate_input_transcript_structures(
                input_transcripts, contig_acc, contig_strand
            )

        ##--------------------------------------------------------------------------------
        # initializes self._splice_graph
        # -defines exons by segmenting genomic coverage based on intron splice coordinates
        # - constructs self._splice_graph as nx.DiGraph()
        self._build_draft_splice_graph()

        if self.is_empty():
            return None

        if LRAA_Globals.DEBUG:
            self.write_intron_exon_splice_graph_bed_files("__prefilter_r1", pad=0)
            self.describe_graph("__prefilter_r1.graph")

        ##----------------------------------------------
        ## Refine the splice graph
        ## removes exon segments, not introns, also removes unspliced introns if Splice_graph._remove_unspliced_introns flag is set.

        if not quant_mode:
            self._prune_lowly_expressed_intron_overlapping_exon_segments()

        self._merge_neighboring_proximal_unbranched_exon_segments()

        if Splice_graph._remove_unspliced_introns and not quant_mode:
            self._prune_unspliced_introns()

        self._prune_low_support_introns()

        self._prune_disconnected_introns()

        self._merge_neighboring_proximal_unbranched_exon_segments()

        # populates self._itree_exon_segments for overlap queries
        self._finalize_splice_graph()

        ## incorporate TSS and PolyA features
        if len(self._TSS_objs) > 0:
            self._incorporate_TSS()

        if len(self._PolyA_objs) > 0:
            self._incorporate_PolyAsites()

        self._finalize_splice_graph()  # do again after TSS and PolyA integration

        if LRAA_Globals.DEBUG:
            self.write_intron_exon_splice_graph_bed_files("__prefilter_r2", pad=0)
            self.describe_graph("__prefilter_r2.graph")

        connected_components = list(
            nx.connected_components(self._splice_graph.to_undirected())
        )

        # filter out TSS and PolyA features based on min isoform fraction
        if len(self._TSS_objs) > 0 or len(self._PolyA_objs) > 0:

            for connected_component in connected_components:
                if len(self._TSS_objs) > 0:
                    self._eliminate_low_support_TSS(connected_component)
                if len(self._PolyA_objs) > 0:
                    self._eliminate_low_support_PolyA(connected_component)

            # revise again
            self._merge_neighboring_proximal_unbranched_exon_segments()

            if LRAA_Globals.DEBUG:
                self.write_intron_exon_splice_graph_bed_files("__prefilter_r3", pad=0)
                self.describe_graph("__prefilter_r3.graph")

            # if not quant_mode:
            self._prune_exon_spurs_at_introns()

            self._finalize_splice_graph()  # do again after TSS and PolyA integration
            connected_components = list(
                nx.connected_components(self._splice_graph.to_undirected())
            )

        else:
            # if not quant_mode:
            self._prune_exon_spurs_at_introns()
            self._finalize_splice_graph()

        if LRAA_Globals.DEBUG:
            self.write_intron_exon_splice_graph_bed_files("__final_graph", pad=0)
            self.describe_graph("__final.graph")
            append_log_file("__TSS_info.bed", self._TSS_objs)
            append_log_file("__PolyAsite_info.bed", self._PolyA_objs)

        self._components = connected_components
        for i, component in enumerate(self._components):
            for node in component:
                id = node.get_id()
                logger.debug(f"assigning node {id} to component {i}")
                self._node_id_to_component[id] = i

        return self._splice_graph

    def _node_has_successors(self, node):

        if len(list(self._splice_graph.successors(node))) > 0:
            return True
        else:
            return False

    def _node_has_predecessors(self, node):

        if len(list(self._splice_graph.predecessors(node))) > 0:
            return True
        else:
            return False

    def _get_exon_and_intron_nodes(self):

        intron_objs = list()
        exon_segment_objs = list()

        for node in self._splice_graph:
            if type(node) == Intron:
                intron_objs.append(node)
            elif type(node) == Exon:
                exon_segment_objs.append(node)
            elif type(node) in (TSS, PolyAsite):
                pass

            else:
                raise RuntimeError(
                    "Error, not identifying node: {} as Exon or Intron type - instead {} ".format(
                        node, type(node)
                    )
                )

        exon_segment_objs = sorted(exon_segment_objs, key=lambda x: x._lend)
        intron_objs = sorted(intron_objs, key=lambda x: x._lend)

        return (exon_segment_objs, intron_objs)

    def _initialize_contig_coverage(self):

        ## Contig Depth Array Capture
        # get genome contig sequence
        contig_seq_str = self._contig_seq_str
        contig_len = self._contig_seq_len
        logging.info("initing coverage array of len: {}".format(contig_len))
        if contig_len > Splice_graph._max_genomic_contig_length:
            raise RuntimeError(
                "genomic contig length {} exceeds maximum allowed {}".format(
                    contig_len, Splice_graph._max_genomic_contig_length
                )
            )

        # init depth of coverage array
        self._contig_base_cov = [0 for i in range(0, contig_len + 1)]

        return

    def _populate_exon_coverage_and_extract_introns(
        self, bam_filename, contig_acc, contig_strand, quant_mode, restrict_splice_type, SE_read_encapsulation_mask
    ):

        ## Intron Capture

        intron_counter = defaultdict(int)

        intron_to_read_types = defaultdict(set)

        intron_splice_site_support = defaultdict(int)

        # get read alignments
        # - illumina and pacbio reads filtered based on tech-specific min per_id
        # - pretty alignments: store the pysam alignment record along with inferred transcript exons segments.
        pretty_alignment_manager = Pretty_alignment_manager(self)
        pretty_alignments = pretty_alignment_manager.retrieve_pretty_alignments(
            contig_acc,
            contig_strand,
            self._contig_seq_str,
            bam_filename,
            region_lend=self._region_lend,
            region_rend=self._region_rend,
            use_cache=False,
            try_correct_alignments=False,
            per_id_QC_raise_error=True,
            restrict_splice_type=restrict_splice_type,
            SE_read_encapsulation_mask=SE_read_encapsulation_mask,
        )

        assert contig_strand in ("+", "-")

        logger.info("-got {} pretty alignments.".format(len(pretty_alignments)))

        total_read_alignments_used = 0

        ## do not count boundaries with soft clipping for TSS sites.  TODO:// examine polyA on reads for polyA sites.
        TSS_position_counter = defaultdict(int)
        polyA_position_counter = defaultdict(int)

        if LRAA_Globals.DEBUG:
            TSS_reads_ofh = open("__TSS_read_support.tsv", "at")
            POLYA_reads_ofh = open("__POLYA_read_support.tsv", "at")

        if contig_strand == "+":
            pretty_alignments = sorted(
                pretty_alignments, key=lambda x: x.get_alignment_span()
            )
        else:
            pretty_alignments = sorted(
                pretty_alignments, key=lambda x: list(reversed(x.get_alignment_span()))
            )

        for pretty_alignment in pretty_alignments:

            read_name = pretty_alignment.get_read_name()
            align_lend, align_rend = pretty_alignment.get_alignment_span()

            TSS_pos, polyA_pos = (
                (align_lend, align_rend)
                if contig_strand == "+"
                else (align_rend, align_lend)
            )
            TSS_pos_soft_clipping = (
                pretty_alignment.left_soft_clipping
                if contig_strand == "+"
                else pretty_alignment.right_soft_clipping
            )
            polyA_pos_soft_clipping = (
                pretty_alignment.right_soft_clipping
                if contig_strand == "+"
                else pretty_alignment.left_soft_clipping
            )

            ######################
            ## Capture TSS support

            if TSS_pos_soft_clipping <= LRAA_Globals.config["max_soft_clip_at_TSS"]:
                TSS_position_counter[TSS_pos] += 1
                if LRAA_Globals.DEBUG:
                    print(
                        "\t".join(
                            [
                                contig_acc,
                                contig_strand,
                                str(TSS_pos),
                                read_name,
                                pretty_alignment.get_pretty_alignment_string(
                                    contig_acc
                                ),
                            ]
                        ),
                        file=TSS_reads_ofh,
                    )

            ########################
            ## Capture PolyA support

            if polyA_pos_soft_clipping <= LRAA_Globals.config["max_soft_clip_at_PolyA"]:
                polyA_position_counter[polyA_pos] += 1
                if LRAA_Globals.DEBUG:
                    print(
                        "\t".join(
                            [
                                contig_acc,
                                contig_strand,
                                str(polyA_pos),
                                read_name,
                                f"polyA_sc:{polyA_pos_soft_clipping}",
                                pretty_alignment.get_pretty_alignment_string(
                                    contig_acc
                                ),
                            ]
                        ),
                        file=POLYA_reads_ofh,
                    )

            alignment_segments = pretty_alignment.get_pretty_alignment_segments()
            # print("Pretty alignment segments: " + str(alignment_segments))
            logger.debug(
                "Pretty alignment for read {} : {}".format(
                    read_name, str(alignment_segments)
                )
            )

            read_type = pretty_alignment.get_read_type()

            if len(alignment_segments) > 1:
                # ensure proper consensus splice sites.
                introns_list = self._get_introns_matching_splicing_consensus(
                    alignment_segments
                )

                # print("introns list: " + str(introns_list))
                logger.debug(
                    "\tread {} : introns found: {}".format(read_name, str(introns_list))
                )
                for intron in introns_list:
                    intron_counter[intron] += 1
                    intron_to_read_types[intron].add(read_type)

                    intron_lend, intron_rend, splice_orient = intron
                    intron_splice_site_support[intron_lend] += 1
                    intron_splice_site_support[intron_rend] += 1

            total_read_alignments_used += 1
            # add to coverage
            for segment in alignment_segments:
                for i in range(segment[0], segment[1] + 1):
                    if i > self._contig_seq_len:
                        break

                    self._contig_base_cov[i] += 1

        logger.info(
            "-total read alignments used: {}".format(total_read_alignments_used)
        )

        # retain only those introns that meet the min threshold
        for intron_coords, count in intron_counter.items():

            read_types = intron_to_read_types[intron_coords]

            if count >= Splice_graph._min_intron_support:
                ## check splice support
                intron_lend, intron_rend, intron_orient = intron_coords
                splice_support_left = intron_splice_site_support[intron_lend]
                splice_support_right = intron_splice_site_support[intron_rend]

                min_support = min(splice_support_left, splice_support_right)
                max_support = max(splice_support_left, splice_support_right)

                frac_intron_support = min_support / max_support
                intron_coords_key = "{}:{}".format(intron_lend, intron_rend)
                # if frac_intron_support >= Splice_graph._min_alt_splice_freq:
                intron_obj = Intron(
                    self._contig_acc, intron_lend, intron_rend, intron_orient, count
                )
                intron_obj.add_read_types(list(read_types))

                logger.debug(
                    "Adding intron {} with frac_intron_support {}".format(
                        intron_coords_key, frac_intron_support
                    )
                )
                self._intron_objs[intron_coords_key] = intron_obj
                # else:
                #    logger.debug("Excluding intron {} with insufficient frac_intron_support {} (below threshold Splice_graph._min_alt_splice_freq {})".format(intron_coords_key, frac_intron_support, Splice_graph._min_alt_splice_freq))

        # Define TSS and PolyA sites
        #    unless it's derived from the input annotation in quant-only mode, in which case we'll limit ourselves to that.
        if not (quant_mode and (len(self._TSS_objs) > 0 or len(self._PolyA_objs) > 0)):

            if LRAA_Globals.config["infer_TSS"]:
                self._incorporate_TSS_objects(
                    contig_acc, contig_strand, TSS_position_counter
                )

            if LRAA_Globals.config["infer_PolyA"]:
                self._incorporate_PolyA_objects(
                    contig_acc, contig_strand, polyA_position_counter
                )

        return

    def _incorporate_TSS_objects(self, contig_acc, contig_strand, TSS_position_counter):

        if LRAA_Globals.DEBUG:
            write_pos_counter_info(
                "__prelim_TSS_raw_counts.tsv",
                TSS_position_counter,
                contig_acc,
                contig_strand,
            )

        TSS_grouped_positions = aggregate_sites_within_window(
            TSS_position_counter,
            LRAA_Globals.config["max_dist_between_alt_TSS_sites"],
            LRAA_Globals.config["min_alignments_define_TSS_site"],
        )

        TSS_grouped_positions = filter_non_peaky_positions(
            TSS_grouped_positions, TSS_position_counter, contig_acc, contig_strand
        )

        for TSS_peak in TSS_grouped_positions:
            position, count = TSS_peak
            self._TSS_objs.append(
                TSS(contig_acc, position, position, contig_strand, count)
            )

        if LRAA_Globals.DEBUG:
            append_log_file("__prefilter_TSS_info.bed", self._TSS_objs)

        return

    def _incorporate_PolyA_objects(
        self, contig_acc, contig_strand, polyA_position_counter
    ):

        if LRAA_Globals.DEBUG:
            write_pos_counter_info(
                "__prelim_polyA_raw_counts.tsv",
                polyA_position_counter,
                contig_acc,
                contig_strand,
            )

        PolyA_grouped_positions = aggregate_sites_within_window(
            polyA_position_counter,
            LRAA_Globals.config["max_dist_between_alt_polyA_sites"],
            LRAA_Globals.config["min_alignments_define_polyA_site"],
        )

        for polyA_site_grouping in PolyA_grouped_positions:
            position, count = polyA_site_grouping
            self._PolyA_objs.append(
                PolyAsite(contig_acc, position, position, contig_strand, count)
            )

        if LRAA_Globals.DEBUG:
            append_log_file("__prefilter_PolyAsite_info.bed", self._PolyA_objs)

        return

    def _integrate_input_transcript_structures(
        self, transcripts, contig_acc, contig_strand
    ):
        """
        Fold in the reference annotations:
        - add introns where they're missing
        - add base coverage where it's missing.
        - if has TSS or PolyA info in the meta, then incorporate with evidence according to TPM val.

        """

        logger.info("Integrating input transcript structures.")

        TSS_evidence_counter = defaultdict(int)
        PolyA_evidence_counter = defaultdict(int)

        ## Exon Coverage and Intron Capture

        for transcript in transcripts:
            orient = transcript.get_orient()
            transcript_id = transcript.get_transcript_id()

            if contig_strand is not None and orient != contig_strand:
                continue

            exon_segments = transcript.get_exon_segments()
            last_rend = None

            logger.debug("Integrating transcript: {}".format(transcript))

            trans_lend, trans_rend = transcript.get_coords()

            # Note: under LRAA MERGE mode, the mp_counter is separately updated for read counts to match min novel count setting.
            #       The polyA and TSS are reconfigured here according to their min settings under MERGE mode.

            if transcript.has_annotated_TSS():
                TSS_coord = trans_lend if orient == "+" else trans_rend
                if transcript.has_annotated_TPM():
                    TSS_evidence_counter[TSS_coord] += round(transcript.get_TPM())
                else:
                    if LRAA_Globals.LRAA_MODE == "MERGE":
                        TSS_evidence_counter[TSS_coord] += LRAA_Globals.config[
                            "min_alignments_define_TSS_site"
                        ]
                    else:
                        TSS_evidence_counter[TSS_coord] = 1

            if transcript.has_annotated_PolyA():
                polyA_coord = trans_rend if orient == "+" else trans_lend
                if transcript.has_annotated_TPM():
                    PolyA_evidence_counter[polyA_coord] += round(transcript.get_TPM())
                else:
                    if LRAA_Globals.LRAA_MODE == "MERGE":
                        PolyA_evidence_counter[polyA_coord] += LRAA_Globals.config[
                            "min_alignments_define_polyA_site"
                        ]
                    else:
                        PolyA_evidence_counter[polyA_coord] = 1

            if (
                LRAA_Globals.config["fracture_splice_graph_at_input_transcript_bounds"]
                is True
            ):
                self._input_transcript_lend_boundaries.add(trans_lend)
                self._input_transcript_rend_boundaries.add(trans_rend)

            # add coverage for exonic region
            for exon_segment in exon_segments:
                lend, rend = exon_segment
                self._input_transcript_exon_coords_itree[lend : rend + 1] = (
                    transcript_id
                )

                logger.debug("-ensuring coverage for exon: {}-{}".format(lend, rend))
                for i in range(lend, rend + 1):
                    if i > self._contig_seq_len:
                        break
                    if self._contig_base_cov[i] < 1:
                        self._contig_base_cov[i] = 1

                # add missing introns.:
                if last_rend is not None:
                    intron_lend = last_rend + 1
                    intron_rend = lend - 1
                    intron_coords_key = "{}:{}".format(intron_lend, intron_rend)

                    if intron_coords_key not in self._intron_objs:
                        intron_obj = Intron(
                            self._contig_acc, intron_lend, intron_rend, orient, 1
                        )
                        intron_obj.add_read_types(["ref_transcript"])
                        self._intron_objs[intron_coords_key] = intron_obj
                        logger.debug("adding intron {}".format(intron_coords_key))
                        self._input_transcripts_introns[intron_coords_key] = intron_obj
                    else:
                        logger.debug(
                            "intron {} already in splice graph".format(
                                intron_coords_key
                            )
                        )
                        # ensure read type represented.
                        intron_obj = self._intron_objs[intron_coords_key]
                        intron_obj.add_read_types(["ref_transcript"])
                        self._input_transcripts_introns[intron_coords_key] = intron_obj

                last_rend = rend

        if len(TSS_evidence_counter) > 0:
            self._incorporate_TSS_objects(
                contig_acc, contig_strand, TSS_evidence_counter
            )

        if len(PolyA_evidence_counter) > 0:
            self._incorporate_PolyA_objects(
                contig_acc, contig_strand, PolyA_evidence_counter
            )

        return

    def _incorporate_TSS(self):

        contig_acc = self._contig_acc
        contig_strand = self._contig_strand
        assert contig_strand in ("+", "-")
        sg = self._splice_graph

        for TSS_obj in self._TSS_objs:

            sg.add_node(TSS_obj)

            TSS_coord, _ = TSS_obj.get_coords()
            exon_intervals = self.get_overlapping_exon_segments(
                TSS_coord, TSS_coord, min_frac_feature_overlap=0.0
            )
            logger.debug("TSS {} overlaps {}".format(TSS_obj, exon_intervals))
            assert (
                len(exon_intervals) <= 1
            ), "Error, TSS {} overlaps multiple intervals: {}".format(
                TSS_obj, exon_intervals
            )

            if len(exon_intervals) == 0:
                logger.warning(
                    "TSS {} candidate has no exon overlap in the splice graph.".format(
                        TSS_obj
                    )
                )
                continue

            exon_interval = exon_intervals[0]
            # split exon interval at TSS
            exon_lend, exon_rend = exon_interval.get_coords()

            exon_coverage = exon_interval.get_read_support()

            exon_predecessors = sg.predecessors(exon_interval)
            exon_successors = sg.successors(exon_interval)

            if contig_strand == "+" and exon_lend == TSS_coord:
                # just add edge
                logger.debug(
                    "Prefixing TSS {} to exon {}".format(TSS_obj, exon_interval)
                )
                sg.add_edge(TSS_obj, exon_interval)

            elif contig_strand == "-" and exon_rend == TSS_coord:
                # just add edge
                logger.debug(
                    "Postfixing TSS {} to exon {}".format(TSS_obj, exon_interval)
                )
                sg.add_edge(exon_interval, TSS_obj)
            else:
                # must split exon:

                if contig_strand == "+":
                    break_lend, break_rend = TSS_coord - 1, TSS_coord
                else:
                    break_lend, break_rend = TSS_coord, TSS_coord + 1

                new_split_exon_left = Exon(
                    contig_acc,
                    exon_lend,
                    break_lend,
                    contig_strand,
                    self._get_mean_coverage(exon_lend, break_lend),
                )
                for exon_predecessor in exon_predecessors:
                    sg.add_edge(exon_predecessor, new_split_exon_left)

                new_split_exon_right = Exon(
                    contig_acc,
                    break_rend,
                    exon_rend,
                    contig_strand,
                    self._get_mean_coverage(break_rend, exon_rend),
                )
                for exon_successor in exon_successors:
                    sg.add_edge(new_split_exon_right, exon_successor)

                sg.add_edge(new_split_exon_left, new_split_exon_right)

                if contig_strand == "+":
                    sg.add_edge(TSS_obj, new_split_exon_right)
                else:
                    sg.add_edge(new_split_exon_left, TSS_obj)

                # remove interval that got split.
                sg.remove_node(exon_interval)
                self._itree_exon_segments.remove(
                    itree.Interval(exon_lend, exon_rend + 1, exon_interval)
                )
                self._itree_exon_segments[exon_lend : break_lend + 1] = (
                    new_split_exon_left
                )
                self._itree_exon_segments[break_rend : exon_rend + 1] = (
                    new_split_exon_right
                )

                logger.debug(
                    "TSS incorporation: Split {} into {} and {}".format(
                        exon_interval, new_split_exon_left, new_split_exon_right
                    )
                )

        return

    def _incorporate_PolyAsites(self):

        contig_acc = self._contig_acc
        contig_strand = self._contig_strand
        assert contig_strand in ("+", "-")
        sg = self._splice_graph

        for polyA_obj in self._PolyA_objs:

            sg.add_node(polyA_obj)

            polyA_coord, _ = polyA_obj.get_coords()
            exon_intervals = self.get_overlapping_exon_segments(
                polyA_coord, polyA_coord, min_frac_feature_overlap=0.0
            )
            logger.debug("PolyA {} overlaps {}".format(polyA_obj, exon_intervals))
            assert (
                len(exon_intervals) <= 1
            ), "Error, PolyA {} overlaps multiple exon intervals: {}".format(
                polyA_obj, exon_intervals
            )

            if len(exon_intervals) == 0:
                logger.warning(
                    "Warning, polyA site {} no longer overlaps an exon candidate in the graph".format(
                        polyA_obj
                    )
                )
                continue

            exon_interval = exon_intervals[0]
            # split exon interval at TSS
            exon_lend, exon_rend = exon_interval.get_coords()

            exon_coverage = exon_interval.get_read_support()

            exon_predecessors = sg.predecessors(exon_interval)
            exon_successors = sg.successors(exon_interval)

            if contig_strand == "+" and exon_rend == polyA_coord:
                # just add edge
                logger.debug(
                    "Postfixing polyA {} to exon {}".format(polyA_obj, exon_interval)
                )
                sg.add_edge(exon_interval, polyA_obj)

            elif contig_strand == "-" and exon_lend == polyA_coord:
                # just add edge
                logger.debug(
                    "Prefixing polyA {} to exon {}".format(polyA_obj, exon_interval)
                )
                sg.add_edge(polyA_obj, exon_interval)
            else:
                # must split exon:

                if contig_strand == "+":
                    break_lend, break_rend = polyA_coord, polyA_coord + 1
                else:
                    break_lend, break_rend = polyA_coord - 1, polyA_coord

                new_split_exon_left = Exon(
                    contig_acc,
                    exon_lend,
                    break_lend,
                    contig_strand,
                    self._get_mean_coverage(exon_lend, break_lend),
                )
                for exon_predecessor in exon_predecessors:
                    sg.add_edge(exon_predecessor, new_split_exon_left)

                new_split_exon_right = Exon(
                    contig_acc,
                    break_rend,
                    exon_rend,
                    contig_strand,
                    self._get_mean_coverage(break_rend, exon_rend),
                )
                for exon_successor in exon_successors:
                    sg.add_edge(new_split_exon_right, exon_successor)

                sg.add_edge(new_split_exon_left, new_split_exon_right)

                if contig_strand == "+":
                    sg.add_edge(new_split_exon_left, polyA_obj)
                else:
                    sg.add_edge(polyA_obj, new_split_exon_right)

                # remove interval that got split.
                sg.remove_node(exon_interval)
                self._itree_exon_segments.remove(
                    itree.Interval(exon_lend, exon_rend + 1, exon_interval)
                )
                self._itree_exon_segments[exon_lend : break_lend + 1] = (
                    new_split_exon_left
                )
                self._itree_exon_segments[break_rend : exon_rend + 1] = (
                    new_split_exon_right
                )

                logger.debug(
                    "PolyA incorporation: Split {} into {} and {}".format(
                        exon_interval, new_split_exon_left, new_split_exon_right
                    )
                )

    def _get_introns_matching_splicing_consensus(self, alignment_segments):

        genome_seq = self._contig_seq_str

        contig_acc = self.get_contig_acc()
        contig_strand = self.get_contig_strand()

        top_strand_agreement_count = 0
        bottom_strand_agreement_count = 0

        introns_list = list()

        num_introns_with_splice_signals = 0
        num_introns_conflicting_splice_signals = 0

        for i in range(len(alignment_segments) - 1):
            seg_left_rend = alignment_segments[i][1]  # exon coord not inclusive
            seg_right_lend = alignment_segments[i + 1][0]  # exon coord inclusive

            intron_lend = seg_left_rend + 1
            intron_rend = seg_right_lend - 1

            splice_type = Intron.check_canonical_splicing(
                intron_lend, intron_rend, genome_seq
            )

            if splice_type is not None:
                num_introns_with_splice_signals += 1
                if splice_type == contig_strand:
                    introns_list.append((intron_lend, intron_rend, splice_type))

                else:
                    num_introns_conflicting_splice_signals += 1
                    logger.debug(
                        "Splice type for intron {}:{}-{}[{}] is not consistent with read alignment orientation: {}".format(
                            contig_acc,
                            intron_lend,
                            intron_rend,
                            splice_type,
                            contig_strand,
                        )
                    )

        # logger.debug("Of the {} intron candidates, {} = {:.1f} % have conflicting orientation to contig.".format(
        #    num_introns_with_splice_signals,
        #    num_introns_conflicting_splice_signals,
        #    num_introns_conflicting_splice_signals/(num_introns_with_splice_signals+0.001)*100))

        return introns_list

    def _build_draft_splice_graph(self):

        ## do some intron pruning
        self._prune_likely_false_introns()

        ## segment genome coverage
        exon_segments = self._segment_exon_by_coverage_n_splicing()

        draft_splice_graph = nx.DiGraph()

        ## add intron nodes.
        lend_to_intron = defaultdict(list)
        rend_to_intron = defaultdict(list)

        for intron in self._intron_objs.values():
            intron_lend, intron_rend = intron.get_coords()
            count = intron.get_read_support()

            lend_to_intron[intron_lend].append(intron)
            rend_to_intron[intron_rend].append(intron)

            draft_splice_graph.add_node(intron)

        ## add exon nodes and connect to introns.
        ## also connect adjacent exon segments.

        prev_exon_obj = None
        for exon in exon_segments:
            exon_lend, exon_rend = exon

            exon_mean_cov = self._get_mean_coverage(exon_lend, exon_rend)

            exon_obj = Exon(
                self._contig_acc,
                exon_lend,
                exon_rend,
                self._contig_strand,
                exon_mean_cov,
            )

            draft_splice_graph.add_node(exon_obj)

            # connect adjacent exon segments.
            if prev_exon_obj is not None:
                if prev_exon_obj._rend + 1 == exon_obj._lend:
                    # adjacent segments.
                    draft_splice_graph.add_edge(prev_exon_obj, exon_obj)

            prev_exon_obj = exon_obj

            ## connect to introns where appropriate
            candidate_splice_left = exon_lend - 1
            if candidate_splice_left in rend_to_intron:
                introns = rend_to_intron[candidate_splice_left]
                for intron_obj in introns:
                    draft_splice_graph.add_edge(intron_obj, exon_obj)

            candidate_splice_right = exon_rend + 1
            if candidate_splice_right in lend_to_intron:
                introns = lend_to_intron[candidate_splice_right]
                for intron_obj in lend_to_intron[candidate_splice_right]:
                    draft_splice_graph.add_edge(exon_obj, intron_obj)

        self._splice_graph = draft_splice_graph

        return

    def _prune_likely_false_introns(self):

        self._prune_introns_too_long()
        self._prune_spurious_introns_shared_boundary("left")
        self._prune_spurious_introns_shared_boundary("right")

    def _prune_introns_too_long(self):

        introns_too_long = list()

        for intron in self._intron_objs.values():
            if intron.has_read_type("ref_transcript"):
                # retain ref introns
                continue

            intron_lend, intron_rend = intron.get_coords()
            intron_len = abs(intron_rend - intron_lend) + 1
            if intron_len > LRAA_Globals.config["max_intron_length"]:
                logger.info(
                    "intron {} exceeds max intron len: {}, targeting removal from splice graph".format(
                        str(intron), LRAA_Globals.config["max_intron_length"]
                    )
                )
                introns_too_long.append(intron)

        if len(introns_too_long) > 0:
            self.purge_introns_from_splice_graph(introns_too_long)

        return

    def _prune_spurious_introns_shared_boundary(self, left_or_right):

        idx = 0 if left_or_right == "left" else 1

        other_idx = 1 if idx == 0 else 0  # the opposite

        introns_shared_coord = defaultdict(list)

        for intron in self._intron_objs.values():
            intron_coord = intron.get_coords()[idx]
            introns_shared_coord[intron_coord].append(intron)

        introns_to_delete = set()

        # see if there are relatively poorly supported junctions
        for intron_list in introns_shared_coord.values():
            intron_list = sorted(intron_list, key=lambda x: x.get_read_support())
            most_supported_intron = intron_list.pop()
            most_supported_intron_abundance = most_supported_intron.get_read_support()
            for alt_intron in intron_list:

                assert (
                    alt_intron != most_supported_intron
                    and alt_intron.get_coords() != most_supported_intron.get_coords()
                )

                assert (
                    alt_intron.get_coords()[idx]
                    == most_supported_intron.get_coords()[idx]
                )

                assert (
                    alt_intron.get_coords()[other_idx]
                    != most_supported_intron.get_coords()[other_idx]
                )

                assert (
                    alt_intron.get_read_support()
                    <= most_supported_intron.get_read_support()
                )

                if alt_intron.has_read_type("ref_transcript"):
                    # retaining ref introns
                    continue

                alt_intron_abundance = alt_intron.get_read_support()
                alt_intron_relative_freq = (
                    alt_intron_abundance / most_supported_intron_abundance
                )

                delta_other_boundary = abs(
                    most_supported_intron.get_coords()[other_idx]
                    - alt_intron.get_coords()[other_idx]
                )

                if alt_intron_relative_freq < Splice_graph._min_alt_splice_freq:
                    logger.debug(
                        "alt intron: {}".format(alt_intron)
                        + " has rel freq: {} and will be purged".format(
                            alt_intron_relative_freq
                        )
                    )
                    introns_to_delete.add(alt_intron)

                # check for splice site aggregation if low per_id alignments allowed
                elif (
                    LRAA_Globals.config["aggregate_adjacent_splice_boundaries"] is True
                    and delta_other_boundary
                    <= LRAA_Globals.config["aggregate_splice_boundary_dist"]
                ):
                    logger.debug(
                        "alt intron: {} is within aggregation distance ({} < {})  of {} and will be purged".format(
                            alt_intron,
                            delta_other_boundary,
                            LRAA_Globals.config["aggregate_splice_boundary_dist"],
                            most_supported_intron,
                        )
                    )
                    introns_to_delete.add(alt_intron)

        logger.info(
            "removing {} low frequency introns with shared {} coord".format(
                len(introns_to_delete), left_or_right
            )
        )

        self.purge_introns_from_splice_graph(introns_to_delete)

        return

    def purge_introns_from_splice_graph(self, introns_to_delete):

        introns_remove_from_splice_graph = set()

        for intron in introns_to_delete:
            intron_coords = intron.get_coords()
            intron_key = "{}:{}".format(intron_coords[0], intron_coords[1])
            intron_obj = self._intron_objs[intron_key]
            if intron_obj.has_read_type("ref_transcript"):
                logger.debug(
                    "-retaining intron {} as having ref_transcript support".format(
                        str(intron_obj)
                    )
                )
            else:
                logger.debug(
                    "removing intron: {} {}".format(
                        intron_key, self._intron_objs[intron_key]
                    )
                )
                del self._intron_objs[intron_key]
                introns_remove_from_splice_graph.add(intron)

        # remove from splice graph altogether.
        if self._splice_graph is not None:
            for intron in introns_remove_from_splice_graph:
                if intron in self._splice_graph:
                    self._splice_graph.remove_node(intron)

        return

    def _prune_weak_splice_neighboring_segments(self):

        #                    I
        #              _____________
        #             /             \
        #   *********/****       ****\************
        #       A       B         C      D
        #

        coverage_window_length = 10  # maybe make this configurable.
        pseudocount = 1

        introns_to_delete = set()

        ofh = None
        if LRAA_Globals.DEBUG:
            ofh = open("__splice_neighbor_cov_ratios.dat", "a")

        for intron in self._intron_objs.values():
            lend, rend = intron.get_coords()
            intron_abundance = intron.get_read_support()

            A_mean_cov = self._get_mean_coverage(
                lend - coverage_window_length, lend - 1
            )
            B_mean_cov = self._get_mean_coverage(lend, lend + coverage_window_length)

            ratio_B_A = (B_mean_cov + pseudocount) / (A_mean_cov + pseudocount)

            C_mean_cov = self._get_mean_coverage(rend - coverage_window_length, rend)
            D_mean_cov = self._get_mean_coverage(
                rend + 1, rend + coverage_window_length
            )

            ratio_C_D = (C_mean_cov + pseudocount) / (D_mean_cov + pseudocount)

            if LRAA_Globals.DEBUG:
                ofh.write(
                    "{}".format(lend)
                    + "\tI:{}".format(intron_abundance)
                    + "\tA:{:.3f}".format(A_mean_cov)
                    + "\tB:{:.3f}".format(B_mean_cov)
                    + "\tB/A:{:.3f}".format(ratio_B_A)
                    + "\n{}".format(rend)
                    + "\tC:{:.3f}".format(C_mean_cov)
                    + "\tD:{:.3f}".format(D_mean_cov)
                    + "\tC/D:{:.3f}".format(ratio_C_D)
                )

        return

    def _get_mean_coverage(self, start, end):

        cov_sum = 0
        cov_len = 0

        for i in range(start, end + 1):
            if i >= 0 and i < len(self._contig_base_cov):
                cov_len += 1
                cov_sum += self._contig_base_cov[i]

        if cov_len > 0:
            mean_cov = cov_sum / cov_len
            return mean_cov
        else:
            logger.warning(
                "coverage requested to position {} extends beyond length of contig {}".format(
                    end, self._contig_acc
                )
            )
            return 0

    def _segment_exon_by_coverage_n_splicing(self):

        left_splice_sites = set()
        right_splice_sites = set()

        for intron in self._intron_objs.values():
            lend, rend = intron.get_coords()
            left_splice_sites.add(lend)
            right_splice_sites.add(rend)

        exon_segments = list()
        exon_seg_start = None
        for i in range(1, len(self._contig_base_cov)):
            if exon_seg_start is None:
                if self._contig_base_cov[i]:
                    # start new exon segment
                    exon_seg_start = i
            else:
                # handle current segment
                if self._contig_base_cov[i] == 0:
                    # stop segment, add to seg list
                    exon_segments.append([exon_seg_start, i - 1])
                    exon_seg_start = None

            # splice breakpoint logic
            if (
                i + 1 in left_splice_sites
                or i + 1 in self._input_transcript_lend_boundaries
            ):
                if exon_seg_start is not None:
                    exon_segments.append([exon_seg_start, i])
                exon_seg_start = None

            if i in right_splice_sites or i in self._input_transcript_rend_boundaries:
                if exon_seg_start is not None:
                    exon_segments.append([exon_seg_start, i])
                exon_seg_start = None

        # get last one if it runs off to the end of the contig
        if exon_seg_start is not None:
            exon_segments.append([exon_seg_start, len(self._contig_base_cov)])

        if LRAA_Globals.DEBUG:
            # write exon list to file
            with open("__exon_regions.init.bed", "a") as ofh:
                for segment in exon_segments:
                    ofh.write(
                        "\t".join(
                            [self._contig_acc, str(segment[0] - 1), str(segment[1] + 1)]
                        )
                        + "\n"
                    )

            with open("__introns.init.bed", "a") as ofh:
                for intron in self._intron_objs.values():
                    intron_lend, intron_rend = intron.get_coords()
                    intron_support = intron.get_read_support()

                    ofh.write(
                        "\t".join(
                            [
                                self._contig_acc,
                                str(intron_lend),
                                str(intron_rend),
                                "{}-{}".format(intron_lend, intron_rend),
                                str(intron_support),
                            ]
                        )
                        + "\n"
                    )

        return exon_segments

    def _prune_lowly_expressed_intron_overlapping_exon_segments(self):

        draft_splice_graph = self._splice_graph

        intron_objs = list()
        exon_segment_objs = list()

        for node in draft_splice_graph:
            if type(node) == Intron:
                intron_objs.append(node)
            elif type(node) == Exon:
                exon_segment_objs.append(node)
            else:
                raise RuntimeError(
                    "Error, not identifying node: {} as Exon or Intron type - instead {} ".format(
                        node, type(node)
                    )
                )

        # store splice junction coordinates.
        splice_coords = set()
        for intron_obj in intron_objs:
            lend, rend = intron_obj.get_coords()
            splice_coords.add(lend)
            splice_coords.add(rend)

        # build interval tree for exon segments.

        exon_itree = itree.IntervalTree()
        for exon_seg in exon_segment_objs:
            exon_lend, exon_rend = exon_seg.get_coords()
            exon_itree[exon_lend : exon_rend + 1] = exon_seg

        exons_to_purge = set()

        min_intron_cov_for_filtering = 1 / Splice_graph._min_alt_unspliced_freq + 1

        for intron in intron_objs:

            if intron.get_read_support() < min_intron_cov_for_filtering:
                continue
            if (
                intron.get_feature_length()
                > Splice_graph._max_intron_length_for_exon_segment_filtering
            ):
                continue

            intron_lend, intron_rend = intron.get_coords()
            overlapping_exon_segs = exon_itree[intron_lend : intron_rend + 1]
            # print("Intron: {}".format(intron))

            for overlapping_exon_seg_iv in overlapping_exon_segs:
                # print("\toverlaps: {}".format(overlapping_exon_seg))

                overlapping_exon_seg = overlapping_exon_seg_iv.data

                # see if connected by an intron
                exon_lend, exon_rend = overlapping_exon_seg.get_coords()
                if exon_lend - 1 in splice_coords or exon_rend + 1 in splice_coords:
                    continue

                if (
                    float(overlapping_exon_seg.get_read_support())
                    / float(intron.get_read_support())
                    < Splice_graph._min_alt_unspliced_freq
                ):
                    logger.debug(
                        "-pruning {} as has exon_read_support:{} / intron_read_support:{} < min_alt_unspliced_freq: {}".format(
                            overlapping_exon_seg,
                            overlapping_exon_seg.get_read_support(),
                            intron.get_read_support(),
                            Splice_graph._min_alt_unspliced_freq,
                        )
                    )
                    exons_to_purge.add(overlapping_exon_seg)

        # retain exons of input transcripts
        exons_to_retain = set()
        for exon in exons_to_purge:
            lend, rend = exon.get_coords()
            if len(self._input_transcript_exon_coords_itree[lend : rend + 1]) > 0:
                logger.debug(
                    "Retaining exon segment {}-{} due to overlap with input transcript exon segment".format(
                        lend, rend
                    )
                )
                exons_to_retain.add(exon)

        if exons_to_retain:
            exons_to_purge = exons_to_purge - exons_to_retain

        logger.info(
            "-removing {} lowly expressed exon segments based on intron overlap".format(
                len(exons_to_purge)
            )
        )
        if exons_to_purge:
            draft_splice_graph.remove_nodes_from(exons_to_purge)

        if LRAA_Globals.DEBUG:
            exons_to_purge = list(exons_to_purge)
            exons_to_purge = sorted(exons_to_purge, key=lambda x: x._lend)
            with open("__exon_segments_to_purge.bed", "a") as ofh:
                for exon in exons_to_purge:
                    ofh.write(exon.get_bed_row(pad=1) + "\n")

        return

    def _prune_disconnected_introns(self):

        draft_splice_graph = self._splice_graph

        introns_to_remove = list()
        for node in draft_splice_graph:
            if type(node) == Intron:
                ## check it has at least one parent and one child
                if (
                    len(list(draft_splice_graph.predecessors(node))) == 0
                    or len(list(draft_splice_graph.successors(node))) == 0
                ):

                    introns_to_remove.append(node)

        logger.info(
            "-pruning {} now disconnected introns".format(len(introns_to_remove))
        )

        if LRAA_Globals.DEBUG:
            with open("__pruned_disconnected_introns.bed", "a") as ofh:
                for intron in introns_to_remove:
                    ofh.write(intron.get_bed_row(pad=1) + "\n")

        # remove introns
        for intron_to_remove in introns_to_remove:
            intron_lend, intron_rend = intron_to_remove.get_coords()
            intron_token = "{}:{}".format(intron_lend, intron_rend)
            del self._intron_objs[intron_token]

        draft_splice_graph.remove_nodes_from(introns_to_remove)

        return

    def get_exon_predecessor(self, exon_node):
        assert type(exon_node) == Exon

        predecessors = self._splice_graph.predecessors(exon_node)

        exon_predecessors = list()

        for pred in predecessors:
            if type(pred) == Exon:
                exon_predecessors.append(Exon)

        assert len(exon_predecessors) <= 1

        if len(exon_predecessors) == 1:
            return exon_predecessors[0]
        else:
            return None

    def get_exon_successor(self, exon_node):
        assert type(exon_node) == Exon

        successors = self._splice_graph.successors(exon_node)

        exon_successors = list()

        for succ in successors:
            if type(succ) == Exon:
                exon_successors.append(succ)

        assert len(exon_successors) <= 1

        if len(exon_successors) == 1:
            return exon_successors[0]
        else:
            return None

    def get_introns_attached_to_exon(self, exon_node):
        assert type(exon_node) == Exon

        attached_introns = set()
        for node in self._splice_graph.predecessors(exon_node):
            if type(node) == Intron:
                attached_introns.add(node)

        for node in self._splice_graph.successors(exon_node):
            if type(node) == Intron:
                attached_introns.add(node)

        return attached_introns

    def _prune_low_support_introns(self):

        draft_splice_graph = self._splice_graph

        introns_to_prune = set()

        for node in draft_splice_graph:
            # walk each exon island, capture all linked introns

            if type(node) == Exon and self.get_exon_predecessor(node) is None:

                logger.debug("INIT EXON FOUND: {}".format(node))

                incidental_introns = self.get_introns_attached_to_exon(node)

                logger.debug(
                    "INIT EXON incidental introns: {}".format(incidental_introns)
                )

                next_exon = self.get_exon_successor(node)
                while next_exon is not None:

                    logger.debug("NEXT Exon found: {}".format(next_exon))

                    next_incidental_introns = self.get_introns_attached_to_exon(
                        next_exon
                    )
                    if len(next_incidental_introns) > 0:
                        logger.debug(
                            "NEXT incidental introns: {}".format(
                                next_incidental_introns
                            )
                        )
                        incidental_introns.update(next_incidental_introns)

                    next_exon = self.get_exon_successor(next_exon)

                # prune introns in group according to min alt splice freq
                if len(incidental_introns) > 1:
                    logger.debug(
                        "-evaluating {} introns in exon island for intron pruning".format(
                            len(incidental_introns)
                        )
                    )
                    incidental_introns = sorted(
                        incidental_introns,
                        key=lambda x: x.get_read_support(),
                        reverse=True,
                    )
                    top_support = incidental_introns[0].get_read_support()
                    logger.debug(
                        "-top supported intron in exon island has {} read support".format(
                            top_support
                        )
                    )
                    for i in range(1, len(incidental_introns)):
                        intron_i_read_support = incidental_introns[i].get_read_support()
                        frac_i_read_support = intron_i_read_support / top_support
                        logger.debug(
                            "-island neighboring incidental intron[{}]={} has {} read support = {} frac of dominant".format(
                                i,
                                incidental_introns[i],
                                intron_i_read_support,
                                frac_i_read_support,
                            )
                        )
                        if (
                            frac_i_read_support
                            < LRAA_Globals.config["min_alt_splice_freq"]
                        ):
                            # insufficient support, prune all introns with that or lower support.
                            for intron in incidental_introns[i:]:
                                logger.debug(
                                    "-pruning lowly supported intron: {}".format(intron)
                                )
                                introns_to_prune.add(intron)
                            break

        logger.info(
            "-removing {} low supported introns in exon islands according to min alt splice freq.".format(
                len(introns_to_prune)
            )
        )

        self.purge_introns_from_splice_graph(introns_to_prune)

        return

    def write_intron_exon_splice_graph_bed_files(self, output_prefix, pad=0):

        exons_bed_file = "{}.exons.bed".format(output_prefix)
        introns_bed_file = "{}.introns.bed".format(output_prefix)

        exons_ofh = open(exons_bed_file, "a")
        introns_ofh = open(introns_bed_file, "a")

        exons_list = list()
        introns_list = list()

        for node in self._splice_graph:
            if type(node) == Exon:
                exons_list.append(node)
            elif type(node) == Intron:
                introns_list.append(node)
            else:
                # could be TSS or PolyA
                pass
                # raise RuntimeError("not intron or exon object... bug... ")

        exons_list = sorted(exons_list, key=lambda x: x._lend)
        introns_list = sorted(introns_list, key=lambda x: x._lend)

        for exon in exons_list:
            exons_ofh.write(exon.get_bed_row(pad=pad) + "\n")

        for intron in introns_list:
            introns_ofh.write(intron.get_bed_row(pad=pad) + "\n")

        exons_ofh.close()
        introns_ofh.close()

        return

    def describe_graph(self, outputfilename):

        if True:
            return  # disabling for now - files can be too big

        ofh = open(outputfilename, "a")

        nodes = list(self._splice_graph.nodes)

        nodes = sorted(nodes, key=lambda x: x._lend)

        for node in nodes:

            node_descr = self.describe_node(node)
            ofh.write(node_descr + "\n")

        ofh.close()

        return

    def describe_node(self, node):

        node_descr = ""
        preds = list(self._splice_graph.predecessors(node))
        if preds:
            pred_strs = list()
            for pred in preds:
                # print(pred)
                pred_strs.append(str(pred))
            node_descr += ">;<".join(pred_strs)
        else:
            node_descr += "."

        node_descr += "\t<" + str(node) + ">\t"

        succs = list(self._splice_graph.successors(node))
        if succs:
            succs_strs = list()
            for succ in succs:
                succs_strs.append(str(succ))
            node_descr += ">;<".join(succs_strs)
        else:
            node_descr += "."

        return node_descr

    def _merge_neighboring_proximal_unbranched_exon_segments(self):

        merged_node_ids = list()

        #
        #       \/                           \/
        # ------ ------ --------- -------------------
        #       /\                           /\
        #
        #       |~~~~~~~~~~~~~~~~~~~~~~~~~~~~|

        # start at a left-branched exon segment or no predecessors.

        ## identify all exon segments that are not preceded by exon segments
        def branched_left(node):
            assert type(node) == Exon

            lend, rend = node.get_coords()
            if (
                lend in self._input_transcript_lend_boundaries
                or (lend - 1) in self._input_transcript_rend_boundaries
            ):
                return True

            if len(list(self._splice_graph.predecessors(node))) > 1:
                return True

            return False

        def branched_right(node):
            assert type(node) == Exon

            lend, rend = node.get_coords()

            if (
                rend in self._input_transcript_rend_boundaries
                or (rend + 1) in self._input_transcript_lend_boundaries
            ):
                return True

            if len(list(self._splice_graph.successors(node))) > 1:
                return True

            return False

        def branched_or_nothing_left(node):

            if branched_left(node):
                return True

            pred_nodes = list(self._splice_graph.predecessors(node))
            if len(pred_nodes) == 0:
                return True

            for pred_node in pred_nodes:

                if type(pred_node) == Exon:
                    if branched_right(pred_node):
                        return True

                else:
                    # non-exon pred
                    return True

            return False

        def branched_or_nothing_right(node):

            if branched_right(node):
                return True

            succ_nodes = list(self._splice_graph.successors(node))
            if len(succ_nodes) == 0:
                return True

            for succ_node in succ_nodes:

                if type(succ_node) == Exon:
                    if branched_left(succ_node):
                        return True

                else:
                    # non-exon pred
                    return True

            return False

        exon_segment_objs, intron_objs = self._get_exon_and_intron_nodes()

        init_exons = list()
        for exon_segment in exon_segment_objs:
            assert type(exon_segment) == Exon
            if branched_or_nothing_left(exon_segment) and not branched_or_nothing_right(
                exon_segment
            ):
                init_exons.append(exon_segment)

        def get_unbranched_exon_segments(init_node):
            exon_seg_list = [init_node]

            node = init_node
            while not branched_or_nothing_right(node):
                node = next(self._splice_graph.successors(node))
                assert type(node) == Exon
                exon_seg_list.append(node)

            return exon_seg_list

        for init_exon in init_exons:
            logger.debug(
                "Init exon candidate: {}".format(self.describe_node(init_exon))
            )
            assert type(init_exon) == Exon

            exons_to_merge_list = get_unbranched_exon_segments(init_exon)
            if len(exons_to_merge_list) < 2:
                continue

            # do merge (keep first exon, update attributes, then delete the others.
            exons_to_merge_list[0]._rend = exons_to_merge_list[-1]._rend
            exons_to_merge_list[0]._mean_coverage = self._get_mean_coverage(
                exons_to_merge_list[0]._lend, exons_to_merge_list[0]._rend
            )

            logger.debug(
                "Exons to merge: {}".format(
                    "\n".join(
                        [self.describe_node(exon) for exon in exons_to_merge_list]
                    )
                )
            )

            # direct first node to successors of last node.
            for succ_node in self._splice_graph.successors(exons_to_merge_list[-1]):
                self._splice_graph.add_edge(exons_to_merge_list[0], succ_node)

            # prune the intervening nodes.
            self._splice_graph.remove_nodes_from(exons_to_merge_list[1:])

        """
            
        for i in range(1, len(exon_segment_objs)):
            next_node = exon_segment_objs[i]

            if ( (not self._node_has_successors(prev_node))
                and
                (not self._node_has_predecessors(next_node))
                and next_node._lend - prev_node._rend - 1 < Splice_graph._inter_exon_segment_merge_dist):

                # merge next node into the prev node
                prev_node._rend = next_node._rend
                prev_node._mean_coverage = self._get_mean_coverage(prev_node._lend, prev_node._rend)

                for next_node_successor in self._splice_graph.successors(next_node):
                    self._splice_graph.add_edge(prev_node, next_node_successor)

                # remove next node
                merged_node_ids.append(next_node.get_id())
                self._splice_graph.remove_node(next_node)
            else:
                prev_node = next_node
                
        if LRAA_Globals.DEBUG:
            with open("__merged_exons.list", "a") as ofh:
                print("\n".join(merged_node_ids), file=ofh)

        """

    def _prune_exon_spurs_at_introns(self):

        logger.info("checking for exon spurs at introns")

        exon_segment_objs, intron_objs = self._get_exon_and_intron_nodes()

        #######           ---------------------
        #                /       ^intron^      \
        #    -----------/===                 ===\-------------------
        #                R_spur              L_spur

        ## //TODO: incorporate read coverage checks

        def is_R_spur(exon_node):

            logger.debug("Evaluating {} as potential R_spur".format(exon_node))

            if (
                exon_node.get_feature_length()
                > LRAA_Globals.config["max_exon_spur_length"]
            ):
                logger.debug(
                    "{} not R spur, feature length: {} < max_exon_spur_length: {}".format(
                        exon_node,
                        exon_node.get_feature_length(),
                        LRAA_Globals.config["max_exon_spur_length"],
                    )
                )
                return False

            has_successor = self._node_has_successors(exon_node)

            has_intron_predecessor = False
            has_alt_intron = False

            for predecessor in self._splice_graph.predecessors(exon_node):
                if type(predecessor) == Intron:
                    has_intron_predecessor = True

                elif type(predecessor) == Exon:
                    for successor in self._splice_graph.successors(predecessor):
                        if type(successor) == Intron:
                            has_alt_intron = True

            # return (not has_intron_successor) and (not has_intron_predecessor) and has_alt_intron
            R_spur_boolean = (
                (not has_successor) and (not has_intron_predecessor) and has_alt_intron
            )

            logger.debug(
                "{} is R spur == {}, (not has_successor = {}), (not has_intron_predecessor = {}), and (has_alt_intron = {})".format(
                    exon_node,
                    R_spur_boolean,
                    not has_successor,
                    not has_intron_predecessor,
                    has_alt_intron,
                )
            )

            return R_spur_boolean

        #######           ---------------------
        #                /       ^intron^      \
        #    -----------/===                 ===\-------------------
        #                R_spur              L_spur

        def is_L_spur(exon_node):

            logger.debug("Evaluating {} as potential L_spur".format(exon_node))

            if (
                exon_node.get_feature_length()
                > LRAA_Globals.config["max_exon_spur_length"]
            ):
                logger.debug(
                    "{} not L spur, feature length: {} < max_exon_spur_length: {}".format(
                        exon_node,
                        exon_node.get_feature_length(),
                        LRAA_Globals.config["max_exon_spur_length"],
                    )
                )
                return False

            has_predecessor = self._node_has_predecessors(exon_node)

            if has_predecessor:
                logger.debug(
                    "{} has predecessors: {}".format(
                        exon_node, list(self._splice_graph.predecessors(exon_node))
                    )
                )

            has_intron_successor = False
            has_alt_intron = False

            for successor in self._splice_graph.successors(exon_node):
                if type(successor) == Intron:
                    has_intron_successor = True

                elif type(successor) == Exon:
                    for predecessor in self._splice_graph.predecessors(successor):
                        logger.debug(
                            "evaluating if predecessor {} of successor {} to {} is an intron".format(
                                predecessor, successor, exon_node
                            )
                        )
                        if type(predecessor) == Intron:
                            has_alt_intron = True

            L_spur_boolean = (
                (not has_predecessor) and (not has_intron_successor) and has_alt_intron
            )

            logger.debug(
                "{} is L spur == {}, (not has_predecessor = {}), (not has_intron_successor = {}), and (has_alt_intron = {})".format(
                    exon_node,
                    L_spur_boolean,
                    not has_predecessor,
                    not has_intron_successor,
                    has_alt_intron,
                )
            )

            return L_spur_boolean

        exons_to_prune = list()

        for exon in exon_segment_objs:
            if (
                exon.get_feature_length()
                >= Splice_graph._min_terminal_splice_exon_anchor_length
            ):
                logger.debug(
                    "Retaining {} as spur since feature length {} exceeds min_terminal_splice_exon_anchor_length {}".format(
                        exon,
                        exon.get_feature_length(),
                        Splice_graph._min_terminal_splice_exon_anchor_length,
                    )
                )
                # long enough, we'll keep it for now.
                continue

            # ok shortie, must see if it's a spur
            if is_L_spur(exon) or is_R_spur(exon):
                exons_to_prune.append(exon)

        # no pruning if overlaps an input reference transcript
        if len(exons_to_prune) > 0:
            exons_no_overlap_with_ref = list()
            num_retained = 0
            for exon in exons_to_prune:
                exon_lend, exon_rend = exon.get_coords()
                if (
                    len(
                        self._input_transcript_exon_coords_itree[
                            exon_lend : exon_rend + 1
                        ]
                    )
                    == 0
                ):
                    exons_no_overlap_with_ref.append(exon)
                else:
                    num_retained += 1
                    logger.debug(
                        f"-retaining exon spur {exon} {exon_lend}-{exon_rend} due to overlap with ref transcript exon"
                    )

            logger.info(
                f"-retaining {num_retained} exon spurs due to overlap with input transcripts"
            )
            exons_to_prune = exons_no_overlap_with_ref

        if exons_to_prune:
            logger.info("-removing {} exon spurs".format(len(exons_to_prune)))

            if LRAA_Globals.DEBUG:
                with open("__pruned_exon_spurs.list", "a") as ofh:
                    for exon in exons_to_prune:
                        print(exon.get_id(), file=ofh)

            self._splice_graph.remove_nodes_from(exons_to_prune)
        else:
            logger.debug("-no spurs to prune.")

        return

    def _finalize_splice_graph(self):

        self._itree_exon_segments = itree.IntervalTree()

        ## store node ID to node object
        for node in self._splice_graph:
            self._node_id_to_node[node.get_id()] = node
            if type(node) == Exon:
                # store exon, TSS, and PolyA features in itree for overlap queries
                lend, rend = node.get_coords()
                self._itree_exon_segments[lend : rend + 1] = node

            elif type(node) == TSS:
                TSS_coord, _ = node.get_coords()
                max_dist_between_alt_TSS_sites = LRAA_Globals.config[
                    "max_dist_between_alt_TSS_sites"
                ]
                half_dist = int(max_dist_between_alt_TSS_sites / 2)
                self._itree_exon_segments[
                    TSS_coord - half_dist : TSS_coord + half_dist + 1
                ] = node

            elif type(node) == PolyAsite:
                polyAsite_coord, _ = node.get_coords()
                max_dist_between_alt_polyA_sites = LRAA_Globals.config[
                    "max_dist_between_alt_polyA_sites"
                ]
                half_dist = int(max_dist_between_alt_polyA_sites / 2)
                self._itree_exon_segments[
                    polyAsite_coord - half_dist : polyAsite_coord + half_dist + 1
                ] = node

        self._validate_itree()

        ################################
        ## add introns to separate itree

        self._itree_introns = itree.IntervalTree()

        for intron_obj in self._intron_objs.values():
            lend, rend = intron_obj.get_coords()
            self._itree_introns[lend : rend + 1] = intron_obj

        return

    def _validate_itree(self):

        itree = self._itree_exon_segments

        if LRAA_Globals.DEBUG:
            with open("__itree_contents.txt", "wt") as ofh:
                for interval in itree:
                    print(str(interval), file=ofh)

        for node in self._splice_graph:
            if type(node) != Intron:
                # ensure we find it when querying the itree
                lend, rend = node.get_coords()

                overlapping_intervals = itree[lend : rend + 1]
                overlapping_interval_nodes = [i.data for i in overlapping_intervals]
                if node not in overlapping_interval_nodes:
                    raise RuntimeError(
                        "Error, node {} not found among overlapping intervals in itree: {}".format(
                            node, overlapping_intervals
                        )
                    )

        logger.info("itree validates.")

    def _is_unspliced_exon_segment_artifact(self, exon, intron):

        # ensure not overlapping with a reference exon
        exon_lend, exon_rend = exon.get_coords()
        if (
            self._input_transcript_exon_coords_itree[exon_lend : exon_rend + 1]
            is not None
        ):
            # must retain it if it overlaps with a ref exon segment
            return False

        if (
            exon.get_read_support()
            < Splice_graph._min_alt_unspliced_freq * intron.get_read_support()
        ):
            return True

        intron_lend, intron_rend = intron.get_coords()

        # prune singleton exon segments falling in sufficiently expressed introns:
        if (not self._node_has_predecessors(exon)) and (
            not self._node_has_successors(exon)
        ):
            return True

        return False

    def _prune_unspliced_introns(self):

        draft_splice_graph = self._splice_graph

        intron_objs = list()
        exon_segment_objs = list()

        for node in draft_splice_graph:
            if type(node) == Intron:
                intron_objs.append(node)
            elif type(node) == Exon:
                exon_segment_objs.append(node)
            else:
                raise RuntimeError(
                    "Error, not identifying node: {} as Exon or Intron type - instead {} ".format(
                        node, type(node)
                    )
                )

        # build interval tree for exon segments.

        exon_itree = itree.IntervalTree()
        for exon_seg in exon_segment_objs:
            exon_lend, exon_rend = exon_seg.get_coords()
            exon_itree[exon_lend : exon_rend + 1] = exon_seg

        exons_to_purge = set()

        ## should we restrict to certain introns here? min coverage? YES!!!

        for intron in intron_objs:
            intron_lend, intron_rend = intron.get_coords()
            overlapping_exon_segs = exon_itree[intron_lend : intron_rend + 1]

            for overlapping_exon_seg in overlapping_exon_segs:
                exon_seg = overlapping_exon_seg.data
                if self._is_unspliced_exon_segment_artifact(
                    exon=exon_seg, intron=intron
                ):
                    exons_to_purge.add(exon_seg)

        logger.info(
            "-removing {} likely unspliced exon segments based on intron overlap".format(
                len(exons_to_purge)
            )
        )
        if exons_to_purge:
            draft_splice_graph.remove_nodes_from(exons_to_purge)

        if LRAA_Globals.DEBUG:
            exons_to_purge = list(exons_to_purge)
            exons_to_purge = sorted(exons_to_purge, key=lambda x: x._lend)
            with open(
                "__exon_segments_to_purge.bed", "a"
            ) as ofh:  # file should be already created based on earlier low expressed exon segments overlapping introns removal step
                for exon in exons_to_purge:
                    ofh.write(exon.get_bed_row(pad=1) + "\n")

        return

    def is_empty(self):
        return len(self._splice_graph) == 0

    def _eliminate_low_support_TSS(self, node_list):

        logger.debug("Eliminating low support TSS")

        TSS_list = list()

        for node in node_list:
            if type(node) == TSS:
                TSS_list.append(node)

        if TSS_list:

            sum_TSS_read_support = 0
            for TSS_obj in TSS_list:
                sum_TSS_read_support += TSS_obj.get_read_support()

            min_TSS_iso_fraction = LRAA_Globals.config["min_TSS_iso_fraction"]

            TSS_to_purge = list()

            TSS_list = sorted(
                TSS_list, key=lambda x: x.get_read_support(), reverse=True
            )

            for i, TSS_obj in enumerate(TSS_list):
                frac_read_support = TSS_obj.get_read_support() / sum_TSS_read_support
                if frac_read_support < min_TSS_iso_fraction:
                    TSS_to_purge = TSS_list[i:]
                    TSS_list = TSS_list[0:i]
                    break

            if TSS_to_purge:
                logger.debug(
                    "Purging TSSs due to min isoform fraction requirements: {}".format(
                        TSS_to_purge
                    )
                )
                for TSS_obj in TSS_to_purge:
                    self._TSS_objs.remove(TSS_obj)
                self._splice_graph.remove_nodes_from(TSS_to_purge)
                TSS_to_purge.clear()

            # remove remaining potential degradation products
            # walk the splice graph along linear exon connections and prune alt TSSs that have lower than the frac dominant support
            TSS_to_purge = set()  # reinit as set
            max_frac_TSS_is_degradation = LRAA_Globals.config[
                "max_frac_alt_TSS_from_degradation"
            ]  # if neighboring TSS has this frac or less, gets purged as degradation product
            for TSS_obj in TSS_list:
                if TSS_obj in TSS_to_purge:
                    continue

                logger.debug(
                    "Evaluationg TSS for purge as degradation TSS: {}".format(TSS_obj)
                )

                # check if not connected
                if (
                    len(list(self._splice_graph.successors(TSS_obj))) == 0
                    and len(list(self._splice_graph.predecessors(TSS_obj))) == 0
                ):
                    TSS_to_purge.add(TSS_obj)
                    logger.warning(
                        "TSS_obj wasnt connected in the graph... removing it. {}".format(
                            TSS_obj
                        )
                    )
                    continue

                this_TSS_support = TSS_obj.get_read_support()

                connected_exons = (
                    self._splice_graph.successors(TSS_obj)
                    if self._contig_strand == "+"
                    else self._splice_graph.predecessors(TSS_obj)
                )
                connected_exons = [x for x in connected_exons if type(x) == Exon]

                assert (
                    len(connected_exons) == 1
                ), "Error, TSS_obj is not connected to a single exon: {}, connected_exons: {} ".format(
                    TSS_obj, connected_exons
                )
                connected_exon = connected_exons[0]
                have_connected = True
                while have_connected:
                    logger.debug(
                        "Walking exon segments looking for alt TSSs from {}".format(
                            connected_exon
                        )
                    )
                    have_connected = False
                    connected_exons = (
                        self._splice_graph.successors(connected_exon)
                        if self._contig_strand == "+"
                        else self._splice_graph.predecessors(connected_exon)
                    )
                    connected_exons = [x for x in connected_exons if type(x) == Exon]
                    if len(connected_exons) == 1:
                        have_connected = True
                        connected_exon = connected_exons[0]
                        # examine potential alt TSS candidate
                        alt_TSS_candidates = (
                            self._splice_graph.predecessors(connected_exon)
                            if self._contig_strand == "+"
                            else self._splice_graph.successors(connected_exon)
                        )
                        alt_TSS_candidates = [
                            x for x in alt_TSS_candidates if type(x) == TSS
                        ]
                        if len(alt_TSS_candidates) == 1:
                            alt_TSS_candidate = alt_TSS_candidates[0]
                            if (
                                alt_TSS_candidate not in TSS_to_purge
                                and alt_TSS_candidate.get_read_support()
                                / this_TSS_support
                                <= max_frac_TSS_is_degradation
                            ):
                                TSS_to_purge.add(alt_TSS_candidate)
                                logger.debug(
                                    "-purging degradation TSS: {}".format(
                                        alt_TSS_candidate
                                    )
                                )

            if TSS_to_purge:
                logger.debug(
                    "Purging TSSs due to max frac degradation requirements: {}".format(
                        TSS_to_purge
                    )
                )
                for TSS_obj in TSS_to_purge:
                    self._TSS_objs.remove(TSS_obj)
                self._splice_graph.remove_nodes_from(TSS_to_purge)

        return

    def _eliminate_low_support_PolyA(self, node_list):

        logger.debug("Eliminating low support PolyA")

        PolyA_list = list()

        for node in node_list:
            if type(node) == PolyAsite:
                PolyA_list.append(node)

        if PolyA_list:

            sum_PolyA_read_support = 0
            for PolyA_obj in PolyA_list:
                sum_PolyA_read_support += PolyA_obj.get_read_support()

            min_PolyA_iso_fraction = LRAA_Globals.config["min_PolyA_iso_fraction"]

            PolyA_to_purge = list()

            PolyA_list = sorted(
                PolyA_list, key=lambda x: x.get_read_support(), reverse=True
            )

            for i, PolyA_obj in enumerate(PolyA_list):
                frac_read_support = (
                    PolyA_obj.get_read_support() / sum_PolyA_read_support
                )
                if frac_read_support < min_PolyA_iso_fraction:
                    PolyA_to_purge = PolyA_list[i:]
                    PolyA_list = PolyA_list[0:i]
                    break

            if PolyA_to_purge:
                logger.debug(
                    "Purging PolyAs due to min isoform fraction requirements: {}".format(
                        PolyA_to_purge
                    )
                )
                for PolyA_obj in PolyA_to_purge:
                    self._PolyA_objs.remove(PolyA_obj)
                self._splice_graph.remove_nodes_from(PolyA_to_purge)

        return

    def reset_exon_coverage_via_pretty_alignments(self, pretty_alignments):
        self._initialize_contig_coverage()

        # recompute base coverage
        for pretty_alignment in pretty_alignments:
            for pretty_segment in pretty_alignment.get_pretty_alignment_segments():
                lend, rend = pretty_segment
                for i in range(lend, rend + 1):
                    if i < self._contig_seq_len:
                        self._contig_base_cov[i] += 1
                    else:
                        break

        # reassign exon coverage values.
        exon_segment_objs, intron_objs = self._get_exon_and_intron_nodes()
        for exon_seg in exon_segment_objs:
            exon_lend, exon_rend = exon_seg.get_coords()
            exon_seg._mean_coverage = self._get_mean_coverage(exon_lend, exon_rend)

        if LRAA_Globals.DEBUG:
            self.write_intron_exon_splice_graph_bed_files(
                "__final_graph.pretty." + LRAA_Globals.LRAA_MODE, pad=0
            )

        return


## general utility functions used above.


def filter_non_peaky_positions(
    grouped_position_counts, position_counter, contig_acc, contig_strand
):

    grouped_position_counts_kept = list()

    grouped_position_counts = sorted(grouped_position_counts, key=lambda x: x[0])

    if LRAA_Globals.DEBUG:
        ofh = open("__TSS_filter_non_peaky_positions.tsv", "at")

    pseudocount = 1

    window_len = LRAA_Globals.config["TSS_window_read_enrich_len"]
    window_enrichment_factor = LRAA_Globals.config["TSS_window_read_enrich_factor"]

    for grouped_position in grouped_position_counts:
        position, count = grouped_position
        position_real_count = position_counter[position]

        adjacent_counts = [0, 0]
        for i in range(position - window_len, position + window_len):
            if i != position:
                pos_count = position_counter[i]
                if pos_count > 0:
                    adjacent_counts.append(pos_count)

        median_adjacent_count = (
            statistics.median(adjacent_counts) if len(adjacent_counts) > 0 else 0
        )  # adjacent_counts/(2*window_len)
        pos_frac_counts = (position_real_count + pseudocount) / (
            median_adjacent_count + pseudocount
        )

        kept = pos_frac_counts >= window_enrichment_factor

        if LRAA_Globals.DEBUG:
            print(
                "\t".join(
                    [
                        contig_acc,
                        contig_strand,
                        str(position),
                        str(position_real_count),
                        str(median_adjacent_count),
                        str(pos_frac_counts),
                        str(kept),
                    ]
                ),
                file=ofh,
            )

        if kept:
            grouped_position_counts_kept.append(grouped_position)

    return grouped_position_counts_kept


def aggregate_sites_within_window(
    pos_counter, max_distance_between_aggregated_sites, min_count_aggregated_site
):

    position_count_structs = list()

    for position, count in pos_counter.items():

        position_count_struct = {
            "position": int(position),
            "count": count,
            "selected": False,
            "aggregated_count": count,
            "index": -1,  # updated below after sorting
        }

        position_count_structs.append(position_count_struct)

    # first sort by position and identify max within window distance
    position_count_structs = sorted(position_count_structs, key=lambda x: x["position"])

    # aggregate counts within distance from each candidate site
    num_structs = len(position_count_structs)
    for i, i_struct in enumerate(position_count_structs):
        i_struct["index"] = i

        # look left in window
        j = i - 1
        while j >= 0:
            j_struct = position_count_structs[j]
            if (
                i_struct["position"] - j_struct["position"]
                > max_distance_between_aggregated_sites
            ):
                break
            i_struct["aggregated_count"] += j_struct["count"]
            j = j - 1

        # look right in window
        k = i + 1
        while k < num_structs:
            k_struct = position_count_structs[k]
            if (
                k_struct["position"] - i_struct["position"]
                > max_distance_between_aggregated_sites
            ):
                break
            i_struct["aggregated_count"] += k_struct["count"]
            k += 1

    # now capture the peaks
    peak_sites = list()
    agg_count_sorted_position_count_structs = sorted(
        position_count_structs,
        key=lambda x: (x["count"], x["aggregated_count"]),
        reverse=True,
    )

    ## reset aggregated counts
    for i_struct in agg_count_sorted_position_count_structs:
        if i_struct["selected"]:
            # already part of another defined peak
            continue

        i_struct["aggregated_count"] = i_struct["count"]  # reset
        i = i_struct["index"]

        # look left in window
        j = i - 1
        while j >= 0:
            j_struct = position_count_structs[j]
            assert i_struct["position"] > j_struct["position"]
            if (
                i_struct["position"] - j_struct["position"]
                > max_distance_between_aggregated_sites
            ):
                break
            if not j_struct["selected"]:
                i_struct["aggregated_count"] += j_struct["count"]
                j_struct["selected"] = True
            j -= 1

        # look right in window
        k = i + 1
        while k < num_structs:
            k_struct = position_count_structs[k]
            assert k_struct["position"] > i_struct["position"]
            if (
                k_struct["position"] - i_struct["position"]
                > max_distance_between_aggregated_sites
            ):
                break
            if not k_struct["selected"]:
                i_struct["aggregated_count"] += k_struct["count"]
                k_struct["selected"] = True
            k += 1

        if i_struct["aggregated_count"] >= min_count_aggregated_site:
            peak_sites.append([i_struct["position"], i_struct["aggregated_count"]])

    return peak_sites


def append_log_file(filename, genome_features_list):

    with open(filename, "a") as ofh:
        for feature in genome_features_list:
            print(feature.get_bed_row(), file=ofh)


def write_pos_counter_info(filename, position_counter, contig_acc, contig_strand):

    position_counts = list(position_counter.items())
    position_counts = sorted(position_counts, key=lambda x: x[0])

    with open(filename, "at") as ofh:
        for position, count in position_counts:
            print(
                "\t".join([contig_acc, str(position), contig_strand, str(count)]),
                file=ofh,
            )

    return


#############
## unit tests
#############


def test_aggregate_sites_within_window():

    pos_counter = {10: 1, 15: 2, 20: 3, 40: 1, 45: 5, 50: 2, 75: 3}

    peaks = aggregate_sites_within_window(pos_counter, 5, 3)

    expected_peaks = [[45, 8], [15, 6], [75, 3]]

    assert (
        peaks == expected_peaks
    ), "Error, peaks {} differs from expected peaks {}".format(peaks, expected_peaks)
