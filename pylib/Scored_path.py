import sys, os, re
import MultiPath
import MultiPathGraph
import MultiPathGraphNode
import GenomeFeature
import Simple_path_utils
import Transcript
import LRAA_Globals

import math
import logging

logger = logging.getLogger(__name__)


class Scored_path:

    def __init__(self, path_list_of_multipath_graph_nodes):

        self._all_represented_mpgns = (
            set()
        )  # stores all input mpgns and their contained mpgns

        def recursively_capture_nodes(mpgn):
            assert type(mpgn) == MultiPathGraphNode.MultiPathGraphNode
            if mpgn not in self._all_represented_mpgns:
                self._all_represented_mpgns.add(mpgn)
                for contained_mpgn in mpgn.get_containments():
                    recursively_capture_nodes(contained_mpgn)

        for mpgn in path_list_of_multipath_graph_nodes:
            recursively_capture_nodes(mpgn)

        self._all_represented_read_names = set()
        for mpgn in self._all_represented_mpgns:
            self._all_represented_read_names.update(mpgn.get_read_names())

        self._mpgn_list_path = path_list_of_multipath_graph_nodes

        self._multiPath_obj = MultiPath.MultiPath.multiPath_from_mpgn_list(
            self._mpgn_list_path
        )

        self._cdna_len = self._multiPath_obj.get_cdna_length()

        span_lend, span_rend = self._multiPath_obj.get_coords()

        self._contig_span_len = span_rend - span_lend + 1

        # init
        self._score = -1
        self._initial_score = -1

        score = self.compute_path_score()

        # set
        self._score = score
        self._initial_score = score

        self._validate_compatible_containments()

    def __repr__(self):
        txt = "Scored_path: {} (Score={:.5f}, InitScore={:.5f})\nmpgns:\n".format(
            self.get_multiPath_obj(), self.get_score(), self.get_initial_score()
        )

        for mpgn in self.get_path_mpgn_list():
            txt += str(mpgn) + "\n"

        return txt

    def get_score(self):
        return self._score

    def get_initial_score(self):
        return self._initial_score

    def get_path_mpgn_list(self):
        return list(self._mpgn_list_path)

    def get_simple_path(self):
        return self.get_multiPath_obj().get_simple_path()

    def get_all_represented_mpgns(self, additional_mpgns_to_check=None):

        represented_mpgns = set(self._all_represented_mpgns)

        if additional_mpgns_to_check:

            scored_simple_path = self.get_multiPath_obj().get_simple_path()

            for mpgn in additional_mpgns_to_check:
                sg = mpgn.get_splice_graph()
                mpgn_simple_path = mpgn.get_simple_path()
                if Simple_path_utils.simple_path_A_contains_and_compatible_with_simple_path_B_spacer_aware_both_paths(
                    sg, scored_simple_path, mpgn_simple_path
                ):
                    represented_mpgns.add(mpgn)

        return represented_mpgns

    def get_multiPath_obj(self):
        return self._multiPath_obj

    def incompatibility_detected(self, extension_mpgn):

        mpg = extension_mpgn.get_multipath_graph()

        for mpgn in self.get_path_mpgn_list():
            if mpg.incompatible_mpgn_pair(mpgn, extension_mpgn):
                return True

        return False

    def create_scored_path_extension(self, extension_mpgn):

        path_list = self.get_path_mpgn_list() + [extension_mpgn]

        extension_scored_path = Scored_path(path_list)

        return extension_scored_path

    def rescore(self, exclude_read_names=set()):
        self._score = self.compute_path_score(exclude_read_names)

        if self._score > self._initial_score:
            raise RuntimeError(
                "Error, rescored path exceeds initial score for path: " + str(self)
            )

        return

    def get_all_represented_reads(self):

        read_names = set()

        for mpgn in self._all_represented_mpgns:
            for read_name in mpgn.get_read_names():
                read_names.add(read_name)
            # has_containments is a method; call it to avoid always-True evaluation
            if mpgn.has_containments():
                for contained_mpgn in mpgn.get_containments():
                    for read_name in contained_mpgn.get_read_names():
                        read_names.add(read_name)

        return read_names

    def toTranscript(self):

        mpgn_list = self.get_path_mpgn_list()

        mpgn_list = sorted(mpgn_list, key=lambda x: x._rend)

        splice_graph = mpgn_list[0].get_splice_graph()

        orient = splice_graph.get_contig_strand()

        read_names = self.get_all_represented_reads()

        # merge to a single multipath object
        simple_path_list = list()
        for mpgn in mpgn_list:
            mp = mpgn.get_multiPathObj()
            simple_path_list.append(mp.get_simple_path())

        transcript_mp = MultiPath.MultiPath(
            splice_graph, simple_path_list, read_names=read_names
        )

        transcript_obj = transcript_mp.toTranscript()

        return transcript_obj

    def get_all_represented_read_names(self):
        # returns copy of set
        return self._all_represented_read_names.copy()

    def compute_path_score(self, exclude_read_names=set()):

        assert self._cdna_len > 0 and self._contig_span_len > 0

        score = 0

        # Primary scoring: unique supporting read names across represented nodes
        for read_name in self._all_represented_read_names:
            if read_name not in exclude_read_names:
                score += 1

        # Fallback (initial scoring only): if no names are available at all (e.g., external read-name
        # store not populated or inaccessible) and we are NOT rescoring with exclusions, approximate
        # support by summed node counts across represented mpgns. Critically, do not apply this fallback
        # during rescoring with exclude_read_names, otherwise overlapping candidates may retain positive
        # scores after their supporting reads were already assigned, defeating the greedy exclusion logic.
        if score == 0 and not exclude_read_names:
            try:
                approx = 0
                for mpgn in self._all_represented_mpgns:
                    c = 0
                    try:
                        c = int(mpgn.get_count())
                    except Exception:
                        c = 0
                    approx += max(0, c)
                score = approx
            except Exception:
                # keep score at 0 if anything goes wrong
                pass

        logger.debug(
            str(self.get_simple_path()) + "\n^^^ computed with score = {}".format(score)
        )

        return score

    def _validate_compatible_containments(self):

        scored_simple_path = self.get_multiPath_obj().get_simple_path()

        for mpgn in self.get_path_mpgn_list():
            sg = mpgn.get_splice_graph()
            mpgn_simple_path = mpgn.get_simple_path()

            if Simple_path_utils.simple_path_A_contains_and_compatible_with_simple_path_B_spacer_aware_both_paths(
                sg, scored_simple_path, mpgn_simple_path
            ):
                logger.debug(
                    "Validated scored path: {}\ncontains path {}\n".format(
                        scored_simple_path, mpgn_simple_path
                    )
                )
            else:
                raise RuntimeError(
                    "Error, scored path: {}\ndoes not contain path {}\n".format(
                        scored_simple_path, mpgn_simple_path
                    )
                )
