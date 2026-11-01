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

        # Aggregate compact read IDs for internal scoring
        self._all_represented_read_ids = set()
        for mpgn in self._all_represented_mpgns:
            self._all_represented_read_ids.update(mpgn.get_read_ids())

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

    def rescore(self, exclude_read_ids=set()):
        self._score = self.compute_path_score(exclude_read_ids)

        if self._score > self._initial_score:
            raise RuntimeError(
                "Error, rescored path exceeds initial score for path: " + str(self)
            )

        return

    def get_all_represented_reads(self):
        """
        Backward-compatible alias returning read NAMES for display purposes.
        Internally, IDs are preferred. This resolves IDs -> names using READ_NAME_STORE.
        """
        return self.get_all_represented_read_names()

    def get_all_represented_read_ids(self):
        # returns copy of set of IDs
        return set(self._all_represented_read_ids)

    def get_all_represented_read_names(self):
        # Resolve IDs -> names for display/output
        names = set()
        try:
            import LRAA_Globals as LG
            name_store = getattr(LG, "READ_NAME_STORE", None)
            if name_store is not None:
                for rid in self._all_represented_read_ids:
                    nm = name_store.get_name(int(rid))
                    if nm is not None:
                        names.add(nm)
                return names
        except Exception:
            pass
        # Fallback: stringify IDs
        return set(str(rid) for rid in self._all_represented_read_ids)

    def toTranscript(self):

        mpgn_list = self.get_path_mpgn_list()

        mpgn_list = sorted(mpgn_list, key=lambda x: x._rend)

        splice_graph = mpgn_list[0].get_splice_graph()

        orient = splice_graph.get_contig_strand()

        # Prefer passing IDs into MultiPath; it will store IDs in-memory
        read_ids = self.get_all_represented_read_ids()

        # merge to a single multipath object
        simple_path_list = list()
        for mpgn in mpgn_list:
            mp = mpgn.get_multiPathObj()
            simple_path_list.append(mp.get_simple_path())

        transcript_mp = MultiPath.MultiPath(
            splice_graph, simple_path_list, read_names=read_ids
        )

        transcript_obj = transcript_mp.toTranscript()

        return transcript_obj

    # get_all_represented_read_names provided above as resolver from IDs

    def compute_path_score(self, exclude_read_ids=set()):

        assert self._cdna_len > 0 and self._contig_span_len > 0

        score = 0

        # Primary scoring: unique supporting read IDs across represented nodes.
        #
        # Synthetic reads injected for input transcripts are templates, not
        # observations, and two different failure modes hang on how they are counted:
        #
        #   a) a structure NO read supports must never be reconstructed. Its only ids
        #      are synthetic, and it has to score zero.
        #   b) a supplied model whose real reads were all claimed by a dominant
        #      overlapping isoform must still reach quantification. Scoring it zero
        #      drops it at selection, before anything can judge whether it was
        #      expressed; on chr20 HG002 that cost 110 supplied models, only 7 of
        #      which had no expression once quantified.
        #
        # These are distinguishable: (a) has no real read id at all, (b) has real read
        # ids that happen to be excluded. So the template counts only for a path that
        # carries at least one real read, giving (b) a floor of 1 that greedy exclusion
        # cannot erode while leaving (a) at zero.
        #
        # Selection is not reporting. A model selected on that floor is still dropped
        # after EM by filter_multiexonic_isoforms_by_TPM_threshold if quantification
        # gave it nothing.
        synthetic_read_ids = LRAA_Globals.SYNTHETIC_READ_IDS
        has_real_read = any(
            rid not in synthetic_read_ids for rid in self._all_represented_read_ids
        )
        for rid in self._all_represented_read_ids:
            if rid in exclude_read_ids:
                continue
            if rid in synthetic_read_ids and not has_real_read:
                continue
            score += 1

        # Fallback (initial scoring only): if no ids are available at all (e.g., external read-id
        # store not populated or inaccessible) and we are NOT rescoring with exclusions, approximate
        # support by summed node counts across represented mpgns. Critically, do not apply this fallback
        # during rescoring with exclude_read_ids, otherwise overlapping candidates may retain positive
        # scores after their supporting reads were already assigned, defeating the greedy exclusion logic.
        # The test is whether any ids exist, not whether the score came out zero: a score of zero from
        # ids that were all excluded or all synthetic is a real answer, not missing data.
        if not self._all_represented_read_ids and not exclude_read_ids:
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
