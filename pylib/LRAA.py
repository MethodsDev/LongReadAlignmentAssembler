#!/usr/bin/env python
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
from MultiPath import MultiPath
from MultiPathCounter import MultiPathCounter
from LRAA_Globals import SPACER
import LRAA_Globals
from MultiPathGraph import MultiPathGraph
from Vertex import Vertex
from Transcript import Transcript
import Simple_path_utils
from Pretty_alignment import *
from Pretty_alignment_manager import Pretty_alignment_manager
from Scored_path import Scored_path
from shutil import rmtree
import time
from multiprocessing import Process, Queue
import traceback
from MultiProcessManager import MultiProcessManager
from collections import defaultdict
import Util_funcs
import Simple_path_utils as SPU
import math


logger = logging.getLogger(__name__)


## BIG TODO:// move functionality related to mapping features to the splicegraph to the actual splicegraph class instead of using all the wrapper functions herein.

ITER = 0


class LRAA:

    min_transcript_length = LRAA_Globals.config["min_transcript_length"]
    min_mpgn_read_count = 1

    max_contained_to_be_a_pasa_vertex = 10

    def __init__(self, splice_graph, num_parallel_processes=1):

        self._splice_graph = splice_graph

        self._multipath_graph = None  # set under build_multipath_graph()
        self._mp_counter = None  # set under build_multipath_graph()

        self._contig_acc = None  # set under build_multipath_graph()
        self._contig_strand = None  # set under build_multipath_graph()

        self._num_parallel_processes = num_parallel_processes

        # Map read_name -> (lend, rend) genomic span of the chosen alignment for that read
        # Populated during _populate_read_multi_paths and used to refine transcript terminal bounds.
        self._read_name_to_span = dict()

        return

    def build_multipath_graph(
        self,
        contig_acc,
        contig_strand,
        contig_seq,
        bam_file,
        allow_spacers=False,
        input_transcripts=None,
        restrict_splice_type=None,
        SE_read_encapsulation_mask=None,
    ):

        if SE_read_encapsulation_mask is not None and restrict_splice_type != "SE":
            raise RuntimeError(
                "Error, can only apply SE_read_encapsulation_mask when restricting to SE read alignments"
            )

        # Ensure external read-tracking stores are initialized so MultiPath.get_read_names()
        # can stream names and debug/graph-building logic that depends on names/spans works in assembly mode too.
        try:
            import LRAA_Globals as LG
            if getattr(LG, "READ_NAME_STORE", None) is None or getattr(LG, "MP_READ_ID_STORE", None) is None:
                from ReadNameStore import ReadNameStore  # type: ignore
                from MpReadIdStore import MpReadIdStore  # type: ignore
                # Prefer structured temp under LRAA_TMP_DIR when available
                base_root = os.environ.get("LRAA_TMP_DIR", os.getcwd())
                base_dir = os.path.join(base_root, "__read_tracking")
                os.makedirs(base_dir, exist_ok=True)
                token = f"{contig_acc}.{contig_strand}"
                if getattr(LG, "READ_NAME_STORE", None) is None:
                    LG.READ_NAME_STORE = ReadNameStore(os.path.join(base_dir, f"read_names.{token}"))
                if getattr(LG, "MP_READ_ID_STORE", None) is None:
                    LG.MP_READ_ID_STORE = MpReadIdStore(os.path.join(base_dir, f"mp_reads.{token}"))
        except Exception:
            pass

        start_time = time.time()
        mp_counter = self._populate_read_multi_paths(
            contig_acc,
            contig_strand,
            contig_seq,
            bam_file,
            allow_spacers,
            restrict_splice_type,
            SE_read_encapsulation_mask,
        )
        self._mp_counter = mp_counter

        if input_transcripts is not None:
            logger.info("-incorporating input transcripts into multpath graph")
            self._incorporate_transcripts_into_mp_counter(
                mp_counter, input_transcripts, bam_file
            )

        # Ensure any batched LMDB writes are visible before downstream phases that may read from stores
        try:
            import LRAA_Globals as LG
            for _store in (getattr(LG, "READ_NAME_STORE", None), getattr(LG, "MP_READ_ID_STORE", None)):
                if _store is not None and hasattr(_store, "flush"):
                    try:
                        _store.flush()
                    except Exception:
                        pass
        except Exception:
            pass

        multipath_graph = MultiPathGraph(
            mp_counter,
            self._splice_graph,
            contig_acc,
            contig_strand,
            LRAA.min_mpgn_read_count,
            allow_spacers,
            read_name_to_span=self._read_name_to_span,
        )
        self._multipath_graph = multipath_graph
        self._contig_acc = contig_acc
        self._contig_strand = contig_strand

        if LRAA_Globals.DEBUG:
            ## debugging info
            debug_multipath_graph_filename = "__multipath_graph.{}.{}.dat".format(
                contig_acc, contig_strand
            )
            logger.info("writing {}".format(debug_multipath_graph_filename))
            multipath_graph.describe_graph(debug_multipath_graph_filename)

        build_time = time.time() - start_time
        logger.info("-multipath graph building took {:.1f} seconds.".format(build_time))

        return mp_counter

    def reconstruct_isoforms(self, single_best_only=False):

        # define disjoint graph components.
        mpg = self._multipath_graph

        mpg_components = (
            mpg.define_disjoint_graph_components_via_shared_splice_graph_vertex()
        )

        num_mpg_components = len(mpg_components)

        logger.info("{} connected components identified".format(num_mpg_components))

        if LRAA_Globals.DEBUG:
            mpg.write_mp_graph_nodes_to_gtf("__mpgns.pre.gtf")

        mpg_components = mpg.remove_small_components(
            mpg_components, LRAA.min_transcript_length
        )
        logger.info(
            "{} components surviving the min length {} criterion.".format(
                len(mpg_components), LRAA.min_transcript_length
            )
        )

        if LRAA_Globals.DEBUG:
            mpg.write_mp_graph_nodes_to_gtf("__mpgns.post_length_filter.gtf")

        all_reconstructed_transcripts = list()

        USE_MULTIPROCESSOR = self._num_parallel_processes > 1

        q = None
        mpm = None
        if USE_MULTIPROCESSOR:
            logger.info("-Running assembly jobs with multiprocessing")
            q = Queue()
            mpm = MultiProcessManager(self._num_parallel_processes, q)
        else:
            logger.info(
                "-Running using single thread, so multiprocessing disabled here."
            )  # easier for debugging sometimes

        def get_mpgn_list_coord_span(mpgn_list):

            coords = list()
            for mpgn in mpgn_list:
                mpgn_coords = mpgn.get_coords()
                coords.extend(mpgn_coords)

            coords = sorted(coords)

            return coords[0], coords[-1]

        mpg_component_debug_dir = "__mpg_components"
        if LRAA_Globals.DEBUG:
            if not os.path.exists(mpg_component_debug_dir):
                os.makedirs(mpg_component_debug_dir)

        def write_mpg_component_debug_file(mpgn_list, filename):
            logger.debug("-writing mpgn description file: {}".format(filename))
            with open(filename, "a") as ofh:
                for mpgn in mpgn_list:
                    print(str(mpgn), file=ofh)

        all_reconstructed_transcripts = list()
        component_counter = 0
        for mpg_component in mpg_components:

            mpg_component_size = len(mpg_component)

            component_counter += 1
            coord_span = get_mpgn_list_coord_span(mpg_component)
            logger.info(
                "LRAA - assembly of component {} size {} region: {}{}:{}-{}".format(
                    component_counter,
                    mpg_component_size,
                    self._contig_acc,
                    self._contig_strand,
                    coord_span[0],
                    coord_span[1],
                )
            )

            mpg_token = "{}{}:{}-{}".format(
                self._contig_acc, self._contig_strand, coord_span[0], coord_span[1]
            )
            if LRAA_Globals.DEBUG:
                mpgn_description_filename = "{}/{}.mpgns.txt".format(
                    mpg_component_debug_dir, mpg_token
                )
                write_mpg_component_debug_file(mpg_component, mpgn_description_filename)

            if (
                USE_MULTIPROCESSOR
                and mpg_component_size
                >= LRAA_Globals.config["min_mpgn_component_size_for_spawn"]
            ):
                p = Process(
                    target=self._reconstruct_isoforms_single_component,
                    name=mpg_token,
                    args=(
                        q,
                        mpg_component,
                        component_counter,
                        mpg_token,
                        single_best_only,
                    ),
                )

                mpm.launch_process(p)

            else:

                try:

                    # run directly, not through multiprocessing
                    reconstructed_transcripts = (
                        self._reconstruct_isoforms_single_component(
                            None,
                            mpg_component,
                            component_counter,
                            mpg_token,
                            single_best_only,
                        )
                    )
                    all_reconstructed_transcripts.extend(reconstructed_transcripts)

                except Exception as e:
                    if mpm is not None:
                        mpm.terminate_all_processes()
                    raise (e)

                fraction_jobs_complete = component_counter / num_mpg_components
                logger.info(
                    "progress monitor for {} {} : {:.2f}% of components processed.".format(
                        self._contig_acc, self._contig_strand, fraction_jobs_complete * 100
                    )
                )

        if USE_MULTIPROCESSOR:
            logger.info("WAITING ON REMAINING MULTIPROCESSING JOBS")
            num_failures = mpm.wait_for_remaining_processes()
            logger.info(mpm.summarize_status())
            if num_failures:
                raise RuntimeError(
                    "Error, {} component failures encountered".format(num_failures)
                )

            queue_contents = mpm.retrieve_queue_contents()
            for entry in queue_contents:
                all_reconstructed_transcripts.extend(entry)

        logger.info(
            "Finished round of isoform reconstruction for {} {}".format(
                self._contig_acc, self._contig_strand
            )
        )

        return all_reconstructed_transcripts

    def prune_ref_transcripts_as_evidence(self, transcripts):
        for transcript in transcripts:
            transcript.prune_reftranscript_as_evidence()

    ##################
    ## Private methods
    ##################

    def _reconstruct_isoforms_single_component(
        self, q, mpg_component, component_counter, mpg_token, single_best_only=False
    ):

        start_time = time.time()

        using_multiprocessing = q is not None
        component_size = len(mpg_component)
        logger.info(
            f"ISOFORM RECONSTRUCTION. multiprocessing[{using_multiprocessing}] {mpg_token} component_size: {component_size}"
        )

        contig_acc = self._splice_graph.get_contig_acc()
        contig_strand = self._splice_graph.get_contig_strand()

        gene_id_use = ":".join(
            ["g", contig_acc, contig_strand, "comp-" + str(component_counter)]
        )

        ## exclude those mpgn's that are contained by many transcripts to reduce the trellis size.
        contained_mpg_counter = defaultdict(int)
        for mpgn in mpg_component:
            for mpgn_contained in mpgn.get_containments():
                contained_mpg_counter[mpgn_contained] += 1

        mpg_components_for_trellis = list()

        for mpgn in mpg_component:
            if contained_mpg_counter[mpgn] <= LRAA.max_contained_to_be_a_pasa_vertex:
                mpg_components_for_trellis.append(mpgn)

        mpg_component = mpg_components_for_trellis  # replace for trellis building
        logger.info("-num vertices for trellis: {}".format(len(mpg_component)))

        MIN_SCORE = LRAA_Globals.config["min_path_score"]

        best_transcript_paths = list()

        paths_seen = set()

        def reinit_weights(mpgn_list):
            for mpgn in mpgn_list:
                mpgn.set_reweighted_flag(False)
            return

        round_iter = 0

        mpgns_require_representation = set(mpg_component)

        ###################
        ## Building trellis
        pasa_vertices = self._build_trellis(mpg_component, mpg_token)
        ###################

        if logger.getEffectiveLevel() == logging.DEBUG:  ## for debugging info only
            if round_iter == 1:
                for pasa_vertex in pasa_vertices:
                    logger.debug(pasa_vertex.describe_pasa_vertex())

            self._write_all_scored_paths_to_file(
                component_counter, round_iter, mpg_token, pasa_vertices
            )

        all_scored_paths = self._retrieve_all_scored_paths(pasa_vertices)

        all_scored_paths = sorted(all_scored_paths, key=lambda x: x.get_score())

        all_represented_read_ids = set()

        while len(mpgns_require_representation) > 0 and len(all_scored_paths) > 0:

            round_iter += 1

            # reinit_weights(mpg_component)
            top_scored_path = all_scored_paths.pop()
            assert type(top_scored_path) == Scored_path

            if top_scored_path.get_score() < MIN_SCORE:
                break

            mpgns_represented = top_scored_path.get_all_represented_mpgns(
                additional_mpgns_to_check=mpgns_require_representation
            )

            found_prev_unrepresented_mpgn = False
            for mpgn_represented in mpgns_represented:
                if mpgn_represented in mpgns_require_representation:
                    found_prev_unrepresented_mpgn = True
                    mpgns_require_representation.remove(mpgn_represented)

            if found_prev_unrepresented_mpgn:
                logger.debug(
                    "Retrieved best (Score={})  transcript path for mpg {} : {}".format(
                        top_scored_path.get_score(), mpg_token, top_scored_path
                    )
                )
                if LRAA_Globals.DEBUG:
                    self._write_best_score_path_info_to_file(
                        top_scored_path, round_iter, mpg_token
                    )

                best_transcript_paths.append(top_scored_path)

                all_represented_read_ids.update(
                    top_scored_path.get_all_represented_read_ids()
                )

                # adjust weights
                # self._decrement_transcript_path_vertices(top_scored_path, pasa_vertices)

                for path in all_scored_paths:
                    path.rescore(all_represented_read_ids)

                # reprioritize
                all_scored_paths = sorted(all_scored_paths, key=lambda x: x.get_score())

        # from the best transcript paths, reconstruct the actual transcripts themselves:

        self._validate_pairwise_incompatibilities(best_transcript_paths)

        transcripts = list()

        transcript_counter = 0

        for transcript_path in best_transcript_paths:
            assert type(transcript_path) == Scored_path

            transcript_counter += 1

            transcript_obj = transcript_path.toTranscript()
            transcript_obj.set_gene_id(gene_id_use)

            transcript_id_use = ":".join(
                [
                    "t",
                    contig_acc,
                    contig_strand,
                    "comp-" + str(component_counter),
                    "iso-" + str(transcript_counter),
                ]
            )
            transcript_obj.set_transcript_id(transcript_id_use)

            if transcript_obj is not None:
                transcripts.append(transcript_obj)
                logger.debug("-assembled: {}".format(str(transcript_obj)))

        if LRAA_Globals.config["collapse_alt_TSS_and_PolyA"]:
            transcripts = self._collapse_identical_intron_isoforms(transcripts)

        logger.info(
            "-reconstructed {} transcripts from component {}".format(
                len(transcripts), component_counter
            )
        )

        ###################
        #
        #  Adjust UTRs based on contained path nodes
        #
        #

        def _adjust_transcript_terminals_from_reads(transcript):
            # Gather all read names that supported this transcript
            read_names = set()
            for mp in transcript.get_multipaths_evidence_assigned():
                read_names.update(mp.get_read_names())

            # Collect spans for reads we have recorded; exclude reference/fake reads
            spans = []
            for rn in read_names:
                if rn.startswith("reftranscript:") or rn.startswith("fake_for_merge"):
                    continue
                if rn in self._read_name_to_span:
                    spans.append(self._read_name_to_span[rn])
            if not spans:
                return  # nothing to adjust

            min_lend = min(l for l, _ in spans)
            max_rend = max(r for _, r in spans)

            exons = transcript.get_exon_segments()
            if not exons:
                return

            # Determine allowable terminal exon bounds from splice graph nodes
            sg = self.get_splice_graph()
            sp = transcript.get_simple_path()

            # find first and last exon node objects
            first_exon_bounds = None
            last_exon_bounds = None
            for nid in sp:
                if isinstance(nid, str) and nid.startswith("E:"):
                    node = sg.get_node_obj_via_id(nid)
                    if node is not None:
                        first_exon_bounds = list(node.get_coords())
                        break
            for nid in reversed(sp):
                if isinstance(nid, str) and nid.startswith("E:"):
                    node = sg.get_node_obj_via_id(nid)
                    if node is not None:
                        last_exon_bounds = list(node.get_coords())
                        break

            # Fallback to current exon bounds if node lookup fails
            first_lend, first_rend = exons[0]
            last_lend, last_rend = exons[-1]
            if first_exon_bounds is None:
                first_exon_bounds = [first_lend, first_rend]
            if last_exon_bounds is None:
                last_exon_bounds = [last_lend, last_rend]

            # Adjust toward farthest read bounds but clip to exon node extents
            exon_left_limit = first_exon_bounds[0]
            exon_right_limit = last_exon_bounds[1]
            candidate_first_lend = max(exon_left_limit, min_lend)
            candidate_last_rend = min(exon_right_limit, max_rend)

            # Only adjust ends lacking explicit annotation: TSS and PolyA
            has_TSS = transcript.has_TSS()
            has_PolyA = transcript.has_PolyA()
            strand = transcript.get_strand()

            # Map which side to adjust depending on strand
            # left boundary: '+' => TSS, '-' => PolyA
            # right boundary: '+' => PolyA, '-' => TSS
            adjust_left = (strand == "+" and not has_TSS) or (strand == "-" and not has_PolyA)
            adjust_right = (strand == "+" and not has_PolyA) or (strand == "-" and not has_TSS)

            new_first_lend = candidate_first_lend if adjust_left else first_lend
            new_last_rend = candidate_last_rend if adjust_right else last_rend

            # Sanity: ensure we don't invert exons
            if new_first_lend > first_rend or new_last_rend < last_lend:
                return

            if new_first_lend != first_lend or new_last_rend != last_rend:
                # Apply updates
                exons[0] = [new_first_lend, first_rend]
                exons[-1] = [last_lend, new_last_rend]

                # Write back into transcript and update transcript bounds and cDNA length
                transcript._exon_segments = exons
                transcript._lend = exons[0][0]
                transcript._rend = exons[-1][1]
                transcript._cdna_len = sum(seg[1] - seg[0] + 1 for seg in exons)

        for transcript in transcripts:
            _adjust_transcript_terminals_from_reads(transcript)

        #
        ###################

        # lighten the transcripts before returning them
        for transcript in transcripts:
            transcript.lighten()

        end_time = time.time()
        asm_time_min = (end_time - start_time) / 60

        logger.info(f"Assembly time for {mpg_token} is {asm_time_min:.2f} min.")

        if q is not None:
            # using MultiProcessing Queue
            q.put(transcripts)
        else:
            return transcripts

    def _populate_read_multi_paths(
        self,
        contig_acc,
        contig_strand,
        contig_seq,
        bam_file,
        allow_spacers,
        restrict_splice_type,
        SE_read_encapsulation_mask=None,
    ):
        """
        Reads the alignments from the BAM and for each read traces it
        through a path in the graph.
        The path is stored as a multipath object with a count associated with the number of reads assigned to it.
        """

        global ITER
        ITER += 1

        # distill read alignments into unique multipaths (so if 10k alignments yield the same structure, there's one multipath with 10k count associated)
        mp_counter = MultiPathCounter()

        if bam_file is None:
            return mp_counter  # nothing to do here.

        pretty_alignment_manager = Pretty_alignment_manager(self._splice_graph)

        pretty_alignments = pretty_alignment_manager.retrieve_pretty_alignments(
            contig_acc,
            contig_strand,
            contig_seq,
            bam_file,
            region_lend=self._splice_graph._region_lend,
            region_rend=self._splice_graph._region_rend,
            use_cache=True,
            restrict_splice_type=restrict_splice_type,
            try_correct_alignments=LRAA_Globals.config["try_correct_alignments"],
            SE_read_encapsulation_mask=SE_read_encapsulation_mask,
        )

        # must redo base coverage and exon coverage assignments
        self._splice_graph.reset_exon_coverage_via_pretty_alignments(pretty_alignments)

        # grouping read alignments according to read pairings (for illumina PE data):
        # group alignments:  grouped_alignments['read_name'] = list(read1_pretty_alignment, read2_pretty_alignment, ...)
        grouped_alignments = self._group_alignments_by_read_name(pretty_alignments)

        logger.info(
            "-got {} pretty alignments grouped into {} alignment groups.".format(
                len(pretty_alignments), len(grouped_alignments)
            )
        )

        # capture the read->path assignments:
        if LRAA_Globals.DEBUG:
            read_graph_mappings_ofh = open(f"__read_graph_mappings.dat.{ITER}", "a")

        logger.info("-start: mapping read alignments to the graph")
        num_alignments = len(grouped_alignments)

        # progress config
        show_progress = LRAA_Globals.config.get("show_progress_mapping", True)
        prefer_tqdm = LRAA_Globals.config.get("use_tqdm_progress", True)
        update_every_n = LRAA_Globals.config.get("mapping_update_every_n", 10000)
        update_interval = LRAA_Globals.config.get("mapping_update_interval_sec", 2.0)

        pbar = None
        if show_progress and prefer_tqdm and num_alignments > 0:
            try:
                import importlib

                _mod = importlib.import_module("tqdm")
                _tqdm = getattr(_mod, "tqdm", None)
                if _tqdm is not None:
                    pbar = _tqdm(
                        total=num_alignments,
                        desc="map-reads",
                        unit="reads",
                        leave=False,
                        dynamic_ncols=True,
                    )
            except Exception:
                pbar = None

        prev_time = time.time()
        last_time_update = prev_time
        # periodic INFO logging (separate from stderr progress)
        log_interval = LRAA_Globals.config.get("mapping_log_progress_interval_sec", 30.0)
        last_log_emit = prev_time
        # mapping counters for summary and periodic logging
        num_processed = 0
        num_kept = 0
        num_discard_spacer = 0
        num_no_path = 0
        # optional memory usage helper
        def _rss_mb():
            try:
                import psutil  # type: ignore
                return psutil.Process(os.getpid()).memory_info().rss / (1024.0 * 1024.0)
            except Exception:
                return None
        for i, read_name in enumerate(grouped_alignments):
            num_processed += 1
            if pbar is not None:
                try:
                    pbar.update(1)
                except Exception:
                    pass
            elif show_progress and num_alignments > 0:
                should_update = False
                if update_every_n and update_every_n > 0 and i % update_every_n == 0:
                    should_update = True
                now = time.time()
                if (
                    update_interval
                    and update_interval > 0
                    and now - last_time_update >= update_interval
                ):
                    should_update = True
                if should_update:
                    frac_done = i / num_alignments * 100
                    time_delta = now - prev_time
                    proc = max(i, 1)
                    rate = proc / max(now - prev_time, 1e-6)
                    sys.stderr.write(
                        f"\r[{i} / {num_alignments} =  {frac_done:.2f} rate={rate:.3e} reads/sec    "
                    )
                    sys.stderr.flush()
                    last_time_update = now

            # print("{}\t{}".format(read_name, len(grouped_alignments[read_name])))
            paths_list = list()
            path_candidates = list()  # (path, (lend, rend)) per pretty_alignment
            read_type = None
            for pretty_alignment in grouped_alignments[read_name]:
                if read_type is None:
                    read_type = pretty_alignment.get_read_type()

                path = self._map_read_to_graph(
                    pretty_alignment.get_pretty_alignment_segments()
                )

                logger.debug(
                    "pretty_alignment: {} maps to graph path: {}".format(
                        pretty_alignment, path
                    )
                )
                if path and path != SPACER:
                    assert path[0] != SPACER, "path[0] is SPACER, not allowed"
                    assert path[-1] != SPACER, "path[-1] is SPACER, not allowed"
                    paths_list.append(path)
                    try:
                        span = pretty_alignment.get_alignment_span()
                    except Exception:
                        span = None
                    path_candidates.append((path, span))

            ## not allowing spacers in paths
            paths_list_no_spacers = list()
            candidates_no_spacers = list()
            for idx, path in enumerate(paths_list):
                if SPACER not in path:
                    paths_list_no_spacers.append(path)
                    candidates_no_spacers.append(path_candidates[idx])
                else:
                    num_discard_spacer += 1
                    if LRAA_Globals.DEBUG:
                        read_graph_mappings_ofh.write(
                            "\t".join(
                                [
                                    read_name,
                                    str(grouped_alignments[read_name]),
                                    str(paths_list),
                                    "DISCARDED-SPACER",
                                ]
                            )
                            + "\n"
                        )

            if paths_list_no_spacers:
                # take first one to represent this read
                chosen_path, chosen_span = candidates_no_spacers[0]
                paths_list = [chosen_path]
                num_kept += 1
                if chosen_span is not None:
                    self._read_name_to_span[read_name] = chosen_span
            else:
                num_no_path += 1
                continue

            # Construct MultiPath with optional name storage to reduce memory
            # Always avoid storing read names in memory; rely on external store for names
            store_names = False
            mp = MultiPath(
                self._splice_graph,
                paths_list,
                read_types={read_type},
                read_names={read_name} if store_names else set(),
                read_count=1,
            )

            logger.debug("paths_list: {} -> mp: {}".format(paths_list, mp))

            # Add to multipath counter and capture canonical mp for attaching read ID
            mp_pair = mp_counter.add(mp)
            canonical_mp, _ = mp_pair.get_multipath_and_count()

            # If external read tracking stores are available, persist mp_id -> read_id
            try:
                name_store = getattr(LRAA_Globals, "READ_NAME_STORE", None)
                mp_store = getattr(LRAA_Globals, "MP_READ_ID_STORE", None)
                if name_store is not None and mp_store is not None:
                    rid = name_store.get_or_add(read_name)
                    mp_store.append(canonical_mp.get_id(), rid)
            except Exception:
                pass

            if LRAA_Globals.DEBUG:
                # Report the canonical multipath so that rnames can be resolved via external stores
                read_graph_mappings_ofh.write(
                    "\t".join([read_name, str(pretty_alignment), str(canonical_mp)])
                    + "\n"
                )

            # periodic INFO logs (besides stderr progress bar) for better observability
            if log_interval and log_interval > 0:
                now = time.time()
                if now - last_log_emit >= float(log_interval):
                    frac = (i + 1) / max(num_alignments, 1)
                    rate = (i + 1) / max(now - prev_time, 1e-6)
                    rss = _rss_mb()
                    mem_txt = f" rss={rss:.1f}MB" if rss is not None else ""
                    logger.info(
                        f"[map-reads] {i+1}/{num_alignments} ({frac*100:.2f}%) kept={num_kept} discard_spacer={num_discard_spacer} no_path={num_no_path} rate={rate:.2f}/s{mem_txt}"
                    )
                    last_log_emit = now

        if LRAA_Globals.DEBUG:
            read_graph_mappings_ofh.close()

            with open(f"__mp_counter.{ITER}", "a") as mp_counter_ofh:
                print(str(mp_counter), file=mp_counter_ofh)

        if pbar is not None:
            try:
                pbar.close()
            except Exception:
                pass
        elif show_progress and num_alignments > 0:
            try:
                sys.stderr.write("\n")
                sys.stderr.flush()
            except Exception:
                pass

        # Final summary and completion log
        logger.info(
            f"-done: mapping read alignments to the graph; total={num_processed} kept={num_kept} discard_spacer={num_discard_spacer} no_path={num_no_path}"
        )

        return mp_counter

    def assign_transcripts_paths_in_graph(self, transcripts):

        for transcript in transcripts:
            logger.debug("-mapping transcript to graph: {}".format(transcript))
            segments = transcript.get_exon_segments()
            path = self._map_read_to_graph(
                segments, refine_TSS_simple_path=True, refine_PolyA_simple_path=True
            )
            logger.debug(str(transcript) + " maps to graph as " + str(path))
            assert (
                path is not None
            ), "Error, input transcript {} has no path in graph.".format(
                transcript.get_transcript_id()
            )

            assert (
                SPACER not in path
            ), "Error, found SPACER in input transcript {} with path {}".format(
                transcript.get_transcript_id(), str(path)
            )

            transcript.set_simple_path(path)

        return

    def _group_alignments_by_read_name(self, pretty_alignments):

        grouped_alignments = defaultdict(list)

        for pretty_alignment in pretty_alignments:
            read_name = pretty_alignment.get_read_name()
            grouped_alignments[read_name].append(pretty_alignment)

        return grouped_alignments

    def _map_read_to_graph(
        self,
        alignment_segments,
        refine_TSS_simple_path=True,
        refine_PolyA_simple_path=True,
    ):

        logger.debug("incoming segments: {}".format(alignment_segments))

        path = list()

        num_segments = len(alignment_segments)

        for i in range(num_segments):

            segment = alignment_segments[i]

            path_part = None

            ## determine type of segment
            if i == 0 and num_segments == 1:
                # single exon segment type
                path_part = self._map_segment_to_graph_SINGLE(segment)
            elif i == 0:
                # initial segment
                path_part = self._map_segment_to_graph_INITIAL(segment)
            elif i == num_segments - 1:
                # terminal segment
                path_part = self._get_intron_node_id(alignment_segments[i - 1], segment)
                if not path_part:
                    path_part = [SPACER]
                terminal_segment = self._map_segment_to_graph_TERMINAL(segment)
                if not terminal_segment and path_part != [SPACER]:
                    terminal_segment = [SPACER]
                if terminal_segment:
                    path_part.extend(terminal_segment)
            else:
                # internal segment
                #   first, get preceding intron
                path_part = self._get_intron_node_id(alignment_segments[i - 1], segment)
                if not path_part:
                    path_part = [SPACER]
                internal_segment = self._map_segment_to_graph_INTERNAL(segment)
                if not internal_segment and path_part != [SPACER]:
                    internal_segment = [SPACER]
                if internal_segment:
                    path_part.extend(internal_segment)

            logger.debug("segment: {}  mapped to {}".format(segment, path_part))

            if path_part:
                path.extend(path_part)
                # print("\tpath extended to: {}".format(path))
            else:
                if len(path) == 0 or path[-1] != SPACER:
                    path.append(SPACER)  # spacer

        # refine path

        if SPACER in path:
            path = Simple_path_utils.trim_terminal_spacers(path)

        if refine_TSS_simple_path:
            path = Simple_path_utils.refine_TSS_simple_path(
                self.get_splice_graph(), path
            )

        if refine_PolyA_simple_path:
            path = Simple_path_utils.refine_PolyA_simple_path(
                self.get_splice_graph(), path
            )

        if Simple_path_utils.count_exons_in_simple_path(path) == 0:
            logger.debug("path {} has no exons. Ignoring path.".format(path))
            return None

        if SPACER in path:
            path = self._remove_stutters(path)

        path = Simple_path_utils.add_spacers_between_disconnected_nodes(
            self.get_splice_graph(), path
        )

        # if SPACER in path:
        #    path = self._try_easy_fill_spacers(path)

        return path

    def _get_intron_node_id(self, prev_segment, next_segment):

        intron_lend = prev_segment[1] + 1
        intron_rend = next_segment[0] - 1

        intron_obj = self._splice_graph.get_intron_node_obj(intron_lend, intron_rend)
        if intron_obj:
            return [intron_obj.get_id()]
        else:
            return None

    def _map_segment_to_graph_SINGLE(self, segment):

        contig_strand = self._splice_graph.get_contig_strand()

        overlapping_segments = self._splice_graph.get_overlapping_exon_segments(
            segment[0],
            segment[1],
            min_frac_feature_overlap=LRAA_Globals.config["min_feature_frac_overlap"],
        )

        if contig_strand == "+":
            overlapping_segments = sorted(
                overlapping_segments, key=lambda x: (x._lend, x._weight)
            )
        else:
            overlapping_segments = sorted(
                overlapping_segments, key=lambda x: (x._rend, x._weight)
            )

        path = list()
        for exon_segment in overlapping_segments:
            id = exon_segment.get_id()
            path.append(id)

        return path

    def _map_segment_to_graph_INITIAL(self, segment):

        path = list()

        contig_strand = self._splice_graph.get_contig_strand()

        overlapping_segments = self._splice_graph.get_overlapping_exon_segments(
            segment[0], segment[1]
        )

        if contig_strand == "+":
            overlapping_segments = sorted(
                overlapping_segments, key=lambda x: (x._lend, x._weight)
            )
        else:
            overlapping_segments = sorted(
                overlapping_segments, key=lambda x: (x._rend, x._weight)
            )

        # , min_frac_feature_overlap=LRAA_Globals.config['min_feature_frac_overlap'])

        for i, exon_segment in enumerate(overlapping_segments):
            # check for overlap and not extending beyond feature rend
            if not (
                (
                    segment[0] <= exon_segment._rend
                    and segment[1] >= exon_segment._lend
                    and exon_segment._rend <= segment[1]
                )
            ):

                continue

            # check amount of feature overlap if not splice-adjacent
            if exon_segment._rend == segment[1]:
                path.append(exon_segment.get_id())
            else:
                num_overlap_bases = Util_funcs.get_num_overlapping_bases(
                    [exon_segment._lend, exon_segment._rend], segment
                )
                feature_len = exon_segment._rend - exon_segment._lend + 1
                frac_feature_overlap = num_overlap_bases / feature_len
                if (
                    frac_feature_overlap
                    >= LRAA_Globals.config["min_feature_frac_overlap"]
                ):
                    path.append(exon_segment.get_id())

        return path

    def _map_segment_to_graph_TERMINAL(self, segment):

        path = list()

        contig_strand = self._splice_graph.get_contig_strand()

        overlapping_segments = self._splice_graph.get_overlapping_exon_segments(
            segment[0], segment[1]
        )

        if contig_strand == "+":
            overlapping_segments = sorted(
                overlapping_segments, key=lambda x: (x._lend, x._weight)
            )
        else:
            overlapping_segments = sorted(
                overlapping_segments, key=lambda x: (x._rend, x._weight)
            )

        for exon_segment in overlapping_segments:
            # check for overlap and not extending beyond feature rend
            if not (
                (
                    segment[0] <= exon_segment._rend
                    and segment[1] >= exon_segment._lend
                    and exon_segment._lend >= segment[0]
                )
            ):

                continue

                # check amount of feature overlap if not splice-adjacent
            if exon_segment._lend == segment[0]:
                path.append(exon_segment.get_id())
            else:
                num_overlap_bases = Util_funcs.get_num_overlapping_bases(
                    [exon_segment._lend, exon_segment._rend], segment
                )
                feature_len = exon_segment._rend - exon_segment._lend + 1
                frac_feature_overlap = num_overlap_bases / feature_len
                if (
                    frac_feature_overlap
                    >= LRAA_Globals.config["min_feature_frac_overlap"]
                ):
                    path.append(exon_segment.get_id())

        # logger.debug("segment {} maps to {}".format(segment, path))

        return path

    def _map_segment_to_graph_INTERNAL(self, segment):

        path = list()

        contig_strand = self._splice_graph.get_contig_strand()

        overlapping_segments = self._splice_graph.get_overlapping_exon_segments(
            segment[0],
            segment[1],
            min_frac_feature_overlap=LRAA_Globals.config["min_feature_frac_overlap"],
        )

        if contig_strand == "+":
            overlapping_segments = sorted(
                overlapping_segments, key=lambda x: (x._lend, x._weight)
            )
        else:
            overlapping_segments = sorted(
                overlapping_segments, key=lambda x: (x._rend, x._weight)
            )

        logger.debug("{} overlaps segments: {}".format(segment, overlapping_segments))

        for exon_segment in overlapping_segments:
            # check for overlap only  (prev: and not extending beyond feature rend)
            if segment[0] <= exon_segment._rend and segment[1] >= exon_segment._lend:
                # and
                # segment[0] <= exon_segment._lend and
                # exon_segment._rend <= segment[1]):

                path.append(exon_segment.get_id())

        logger.debug(
            "\t{} restricted to overlapping segments: {}".format(
                segment, overlapping_segments
            )
        )

        return path

    def _remove_stutters(self, path):

        new_path = list()

        for i, node_id in enumerate(path):
            if node_id != SPACER:
                if len(new_path) == 0:
                    new_path.append(node_id)
                elif new_path[-1] != SPACER:
                    if new_path[-1] != node_id:
                        new_path.append(node_id)
                elif (
                    new_path[-1] == SPACER
                    and len(new_path) > 1
                    and new_path[-2] != node_id
                ):
                    new_path.append(node_id)
            # it is a spacer
            elif len(new_path) > 0:
                if new_path[-1] != SPACER:
                    new_path.append(SPACER)

        new_path = Simple_path_utils.trim_terminal_spacers(new_path)

        if new_path != path:
            logger.debug("{} -> {}".format(path, new_path))

        return new_path

    def _try_easy_fill_spacers(self, path):

        return Simple_path_utils.try_fill_spacers_via_splicegraph(
            self._splice_graph, path
        )

    def _build_trellis(self, mpg_component, mpg_token):

        logger.debug("-building trellis for {}".format(mpg_token))

        mpg = self._multipath_graph

        nodes = sorted(
            mpg_component,
            key=lambda x: (
                x._lend,
                x._rend,
                x.get_left_boundary_sort_weight(),
                x.get_right_boundary_sort_weight(),
            ),
        )

        # init the pasa vertex list

        pasa_vertices = list()

        for node in nodes:
            pasa_vertex = Vertex(node)
            pasa_vertices.append(pasa_vertex)

        for i in range(1, len(pasa_vertices)):
            pv_i = pasa_vertices[i]

            for j in range(i - 1, -1, -1):
                pv_j = pasa_vertices[j]

                if mpg.has_edge(
                    pv_j.get_multipath_graph_node(), pv_i.get_multipath_graph_node()
                ):
                    pv_i.add_highest_scoring_path_extension(pv_j)

        return pasa_vertices

    def _retrieve_best_transcript(self, pasa_vertices):

        best_scoring_path = None
        best_score = 0

        for pasa_vertex in pasa_vertices:

            scored_paths = pasa_vertex.get_fromPaths()
            for scored_path in scored_paths:
                if scored_path.get_score() > best_score:
                    best_score = scored_path.get_score()
                    best_scoring_path = scored_path

        logger.debug(
            "-retrieved best transcript path {} with score: {}".format(
                best_scoring_path, best_score
            )
        )

        return best_scoring_path

    def _retrieve_all_scored_paths(self, pasa_vertices):

        all_scored_paths = list()

        for pasa_vertex in pasa_vertices:
            these_scored_paths = pasa_vertex.get_fromPaths()
            all_scored_paths.extend(these_scored_paths)

        return all_scored_paths

    def _decrement_transcript_path_vertices(self, transcript_path, pasa_vertices):

        logger.debug("_decrement_transcript_path_vertices")

        assert type(transcript_path) == Scored_path

        # examine all pasa vertices that are contained and compatible with the transcript_path
        mpgn_list = list()

        transcript_path_multipath_obj = transcript_path.get_multiPath_obj()

        mpgns_not_compatible = list()

        for pasa_vertex in pasa_vertices:
            mpgn = pasa_vertex.get_multipath_graph_node()
            other_multipath_obj = mpgn.get_multiPathObj()
            if transcript_path_multipath_obj.is_overlapping_contained_and_compatible(
                other_multipath_obj
            ):
                mpgn_list.append(mpgn)
            else:
                mpgns_not_compatible.append(mpgn)

        logger.debug(
            "mpgns found compatible with transcript path: {} include {}".format(
                transcript_path_multipath_obj, mpgn_list
            )
        )
        logger.debug("mpgns found INcomptable are: {}".format(mpgns_not_compatible))

        def recursive_reweight(mpgn):
            if mpgn.get_reweighted_flag() is False:
                mpgn.set_weight(0.000001)
                for contained_mpgn in mpgn.get_containments():
                    recursive_reweight(contained_mpgn)

        for mpgn in mpgn_list:
            logger.debug("_decrement: {}".format(mpgn))
            recursive_reweight(mpgn)

            # if mpgn.get_reweighted_flag() is False:
            #    mpgn.reevaluate_weighting_via_path_compatibilities(transcript_path_multipath_obj)

        ## //TODO: //FIXME  why vertex path incompatible but still part of transcript path?

        for mpgn in transcript_path.get_path_mpgn_list():
            recursive_reweight(mpgn)

        return

    def _rescore_transcript_paths(self, pasa_vertices):

        for pasa_vertex in pasa_vertices:
            pasa_vertex.rescore_fromPaths()

        return

    def _write_all_scored_paths_to_file(
        self, component_counter, round_iter, mpg_token, pasa_vertices
    ):

        outdirname = "__all_scored_paths"

        if not os.path.exists(outdirname):
            os.makedirs(outdirname)

        outputfilename = "{}/scored_paths.{}.R{}.gtf".format(
            outdirname, mpg_token, round_iter
        )

        ofh = open(outputfilename, "a")

        for pasa_vertex in pasa_vertices:
            print("## pasa vertex:", file=ofh)
            print(pasa_vertex.describe_pasa_vertex(), file=ofh)
            from_paths = pasa_vertex.get_fromPaths()
            for from_path in from_paths:
                trans_obj = from_path.toTranscript()
                trans_obj.add_meta("score", from_path.get_score())
                gtf = trans_obj.to_GTF_format(include_TPM=False)
                ofh.write(gtf + "\n")
            print("", file=ofh)  # spacer"

        ofh.close()

        return

    def _write_best_score_path_info_to_file(
        self, transcript_path, round_iter, mpg_token
    ):

        outdirname = "__selected_best_scored_paths"

        if not os.path.exists(outdirname):
            try:
                os.makedirs(outdirname)
                # ignore race condition across multiple threads.
            except:
                pass

        outputfilename = "{}/selected_best_path.{}.R{}.gtf".format(
            outdirname, mpg_token, round_iter
        )

        with open(outputfilename, "a") as ofh:
            print("Transcript path: " + str(transcript_path), file=ofh)

    def get_splice_graph(self):
        return self._splice_graph

    def _validate_pairwise_incompatibilities(self, pasa_scored_paths):

        if len(pasa_scored_paths) < 2:
            return

        # If assembly is restricted to collapse, proactively prune contained/duplicate paths
        # instead of raising. This restores pre-refactor behavior where nested paths were
        # collapsed, not treated as a hard error.
        if LRAA_Globals.config.get("restrict_asm_to_collapse", False):
            # Deduplicate identical simple paths, keeping the higher-scored instance
            token_to_idx = dict()
            keep = [True] * len(pasa_scored_paths)
            for idx, spath in enumerate(pasa_scored_paths):
                sp = spath.get_multiPath_obj().get_simple_path()
                token = "^".join(sp)
                prev_idx = token_to_idx.get(token)
                if prev_idx is None:
                    token_to_idx[token] = idx
                else:
                    # keep the higher score, drop the other
                    if pasa_scored_paths[idx].get_score() > pasa_scored_paths[prev_idx].get_score():
                        keep[prev_idx] = False
                        token_to_idx[token] = idx
                    else:
                        keep[idx] = False

            # Prune strictly contained paths (spacer-aware). Prefer keeping the longer path.
            sg = self.get_splice_graph()
            for i in range(len(pasa_scored_paths)):
                if not keep[i]:
                    continue
                sp_A = pasa_scored_paths[i].get_multiPath_obj().get_simple_path()
                for j in range(len(pasa_scored_paths)):
                    if i == j or not keep[j]:
                        continue
                    sp_B = pasa_scored_paths[j].get_multiPath_obj().get_simple_path()

                    # If A contains B, drop B; if B contains A, drop A.
                    if Simple_path_utils.simple_path_A_contains_and_compatible_with_simple_path_B_spacer_aware_both_paths(
                        sg, sp_A, sp_B
                    ) and sp_A != sp_B:
                        keep[j] = False
                    elif Simple_path_utils.simple_path_A_contains_and_compatible_with_simple_path_B_spacer_aware_both_paths(
                        sg, sp_B, sp_A
                    ) and sp_A != sp_B:
                        keep[i] = False
                        break

            # Apply pruning in-place if anything changed
            if not all(keep):
                new_list = [p for k, p in enumerate(pasa_scored_paths) if keep[k]]
                pasa_scored_paths.clear()
                pasa_scored_paths.extend(new_list)
            return

        # Otherwise, treat overlapping and compatible final paths as a hard error (unchanged)
        for i in range(1, len(pasa_scored_paths)):
            mp_A = pasa_scored_paths[i].get_multiPath_obj()
            sp_A = mp_A.get_simple_path()
            sg = mp_A.get_splice_graph()

            for j in range(i - 1, -1, -1):
                mp_B = pasa_scored_paths[j].get_multiPath_obj()
                sp_B = mp_B.get_simple_path()

                if Simple_path_utils.simple_paths_overlap_and_compatible_spacer_aware_both_paths(
                    sg, sp_A, sp_B
                ):
                    raise RuntimeError(
                        "\n\n\n##****************************************************************************\n"
                        + "Error, final paths:\n{} and \n{}\noverlap and are compatible - should have gotten assembled together".format(
                            pasa_scored_paths[i], pasa_scored_paths[j]
                        )
                        + "\n##***************************************************************************************\n\n"
                    )

        return

    def _collapse_identical_intron_isoforms(self, transcripts):

        def get_intron_token(transcript):
            sp = transcript.get_simple_path()
            intron_ids = [x for x in sp if re.match("I:", x)]
            if len(intron_ids) > 0:
                token = "^".join(intron_ids)
                return token
            else:
                return None

        sg = self.get_splice_graph()
        contig_acc = sg.get_contig_acc()
        contig_strand = sg.get_contig_strand()

        transcripts_ret = list()

        intron_tok_to_transcripts = defaultdict(list)

        for transcript in transcripts:
            intron_tok = get_intron_token(transcript)
            if intron_tok is not None:
                intron_tok_to_transcripts[intron_tok].append(transcript)
            else:
                transcripts_ret.append(transcript)

        for intron_tok, transcript_list in intron_tok_to_transcripts.items():
            if len(transcript_list) == 1:
                transcripts_ret.extend(transcript_list)
            else:
                # merge them.
                merged_transcript = transcript_list.pop()
                # retain the original ids here
                trans_id = merged_transcript.get_transcript_id()
                gene_id = merged_transcript.get_gene_id()

                while len(transcript_list) > 0:

                    next_transcript = transcript_list.pop()

                    logger.debug(
                        "Merging ident intron transcript {} into {}".format(
                            merged_transcript, next_transcript
                        )
                    )

                    all_read_names = set(merged_transcript.get_read_names()) | set(
                        next_transcript.get_read_names()
                    )

                    new_transcript_mp = MultiPath(
                        sg,
                        [
                            merged_transcript.get_simple_path(),
                            next_transcript.get_simple_path(),
                        ],
                        read_names=all_read_names,
                    )

                    merged_transcript = new_transcript_mp.toTranscript()

                merged_transcript.set_gene_id(gene_id)
                merged_transcript.set_transcript_id(trans_id)
                transcripts_ret.append(merged_transcript)

        return transcripts_ret

    def _incorporate_transcripts_into_mp_counter(
        self, mp_counter, input_transcripts, bam_file
    ):

        ## if bam_file is None, then working in LRAA transcript merge-only mode.
        ## and should treat each input transcript like it's supported by reads in par with TPM value.

        for input_transcript in input_transcripts:
            simple_path = input_transcript.get_simple_path()
            if bam_file is not None:
                mp = MultiPath(
                    self._splice_graph,
                    [simple_path],
                    read_types={"reftranscript"},
                    read_names={
                        "reftranscript:" + input_transcript.get_transcript_id()
                    },
                )
            else:
                # fake read, transcript-merge mode.
                fake_read_prefix = (
                    input_transcript.get_transcript_id() + "-" + str(time.time())
                )
                if (
                    input_transcript.has_annotated_TPM()
                    and LRAA_Globals.LRAA_MODE != "MERGE"
                ):
                    num_fake_reads = math.ceil(input_transcript.get_TPM())
                else:
                    num_fake_reads = LRAA_Globals.config["min_reads_novel_isoform"]

                fake_read_names = set(
                    [f"{fake_read_prefix}.{i}" for i in range(num_fake_reads)]
                )

                mp = MultiPath(
                    self._splice_graph,
                    [simple_path],
                    read_types={"fake_for_merge"},
                    read_names=fake_read_names,
                )

            mp_counter.add(mp)

        return

    @classmethod
    def differentiate_known_vs_novel_isoforms(
        cls, transcripts, reference_splice_hashes=None
    ):
        """Set novelty flag on transcripts based on splice pattern (introns) only.

        Rationale: We want only transcripts whose splice pattern exactly matches one
        provided among the reference input transcripts (and that are spliced) to be
        treated as known. All others, including ME/SE-derived isoforms (even if they
        may share some multipath evidence with references) are considered novel and
        subjected to min_reads_novel_isoform filtering. Mono-exonic transcripts are
        always classified as novel here per user guidance (restriction to spliced).
        """

        if reference_splice_hashes is not None and not isinstance(
            reference_splice_hashes, set
        ):
            reference_splice_hashes = set(reference_splice_hashes)

        for transcript in transcripts:
            if transcript.has_introns() and reference_splice_hashes is not None:
                splice_hash = Util_funcs.get_hash_code(transcript.get_introns_string())
                if splice_hash in reference_splice_hashes:
                    transcript.set_is_novel_isoform(False)
                    continue
            # default novel
            transcript.set_is_novel_isoform(True)
