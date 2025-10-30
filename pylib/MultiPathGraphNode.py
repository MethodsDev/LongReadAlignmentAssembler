import sys, os, re
import Splice_graph
import MultiPath
import MultiPathGraph
import MultiPathCounter
import networkx as nx
import Simple_path_utils as Simple_path_utils
import LRAA_Globals
from LRAA_Globals import SPACER
from GenomeFeature import *
import math

import logging

logger = logging.getLogger(__name__)


NORMALIZE_SCORES_BY_SEQ_LENGTH = False
MIN_WEIGHT = 0.01


class MultiPathGraphNode:

    mp_id_counter = 0

    def __init__(self, multiPathObj, count, lend, rend, mpg):

        # ensure we have positive support
        assert count > 0, "error, multipath obj lacks read support - shouldnt be possible: " + str(
            multiPathObj
        )

        self._multiPath = multiPathObj
        self._count = count

        self._lend = lend
        self._rend = rend

        # set these if bounded by polyA or TSS
        self._left_boundary = 0
        self._right_boundary = 0

        self._mpg = mpg

        self._containments = set()  # other MPG nodes fully contained by this node.

        MultiPathGraphNode.mp_id_counter += 1
        self._id = "mp{}x".format(MultiPathGraphNode.mp_id_counter)

        self._seq_length = self._compute_seq_length()

        self._weight = 1.0  # contribution of read content towards scores

        self._prev_weight = self._weight  # retain previous weight setting

        self._component_id = (
            0  # will be set to a component ID after components are defined
        )

        self._reweighted_flag = (
            False  # should be initialized to False before each reconstruction round.
        )

        # Optional: min/max genomic span across reads supporting this node (if provided by caller)
        self._min_read_lend = None
        self._max_read_rend = None

    def get_id(self):
        return self._id

    def get_multiPathObj(self):
        return self._multiPath

    def get_multipath_graph(self):
        return self._mpg

    def get_simple_path(self):
        return self._multiPath.get_simple_path()

    def get_read_names(self):

        # return all read names in this mpgn and all contained mpgns

        all_read_names = set()
        all_read_names.update(self._multiPath.get_read_names())

        for contained_mpgn in self.get_containments():
            all_read_names.update(contained_mpgn._multiPath.get_read_names())

        return all_read_names

    def get_read_ids(self):
        """
        Returns the set of compact read IDs across this node and any contained nodes.
        Prefer this in algorithmic paths; resolve names only for display/output.
        """
        all_ids = set()
        try:
            all_ids.update(self._multiPath.get_read_ids())
        except Exception:
            pass
        for contained_mpgn in self.get_containments():
            try:
                all_ids.update(contained_mpgn._multiPath.get_read_ids())
            except Exception:
                continue
        return all_ids

    def set_reweighted_flag(self, flag_setting):
        self._reweighted_flag = flag_setting
        return

    def get_reweighted_flag(self):
        return self._reweighted_flag

    def get_coords(self):
        return (self._lend, self._rend)

    def set_read_span_bounds(self, lend, rend):
        self._min_read_lend = lend
        self._max_read_rend = rend

    def get_read_span_bounds(self):
        if self._min_read_lend is None or self._max_read_rend is None:
            return None
        return (self._min_read_lend, self._max_read_rend)

    def get_count(self):
        return self._count

    def set_count(self, count_val):
        self._count = count_val
        return

    def get_weight(self):
        return self._weight

    def set_weight(self, weight):
        if self.get_reweighted_flag() is True:
            raise RuntimeError(
                "Error, cant set weight on mpgn when reweighted flag is already set"
            )

        self._prev_weight = self._weight

        self._weight = weight
        self.set_reweighted_flag(True)

        logger.debug(
            "changed mpgn weight from {} to {}".format(self._prev_weight, self._weight)
        )
        # sys.stderr.write("changed mpgn weight from {} to {}".format(self._prev_weight, self._weight))

        return

    def get_left_boundary_sort_weight(self):
        if self._left_boundary == 1:
            return 1
        else:
            return 0

    def get_right_boundary_sort_weight(self):
        if self._right_boundary == 1:
            return 1
        else:
            return 0

    def get_prev_weight(self):
        return self._prev_weight

    def get_score_EXCLUDE_containments(self, use_prev_weight=False):
        weight = self._prev_weight if use_prev_weight else self._weight
        score = self._count * weight
        if NORMALIZE_SCORES_BY_SEQ_LENGTH:
            score = score / self._seq_length  # normalize read count by feature length

        return score

    def get_score_INCLUDE_containments(self, use_prev_weight=False, mpgn_ignore=set()):

        seen = set(mpgn_ignore)

        total_counts = 0

        all_relevant_nodes = [self]
        contained_nodes = self.get_containments()
        if contained_nodes:
            all_relevant_nodes.extend(contained_nodes)

        for node in all_relevant_nodes:
            assert type(node) == MultiPathGraphNode
            if node not in seen:
                weight = node._prev_weight if use_prev_weight else node._weight
                total_counts += node._count * weight

        score = total_counts

        if NORMALIZE_SCORES_BY_SEQ_LENGTH:
            score = score / self._seq_length  # normalize read count by feature length

        return score

    def get_seq_length(self):
        return self._seq_length

    def get_component_id(self):
        return self._component_id

    def set_component_id(self, component_id):
        self._component_id = component_id

    def __repr__(self):
        return self.toString(recursive=True)  # to get first round of containments

    def toString(self, recursive=False):
        containments = self.get_containments()

        # Abbreviate read names for logging to avoid huge outputs
        try:
            names = self.get_read_names()
        except Exception:
            names = set()
        if len(names) > 10:
            read_names_show = str(list(names)[0:10]) + "....{} num reads".format(
                len(names)
            )
        else:
            read_names_show = str(names)

        text = "<{}:{} {}-{} Count:{} W:{:0.8f} Containments:{}, ScoreExcCont:{:.4f} ScoreInclCon:{:.4f} len:{} read_names:{}>".format(
            self.get_id(),
            self.get_simple_path(),
            self._lend,
            self._rend,
            self._count,
            self._weight,
            len(containments),
            self.get_score_EXCLUDE_containments(use_prev_weight=False),
            self.get_score_INCLUDE_containments(use_prev_weight=False),
            self._seq_length,
            read_names_show,
        )

        if recursive:
            for containment in containments:
                text += "\n\tcontained: " + containment.toString(recursive=False)

        return text

    def has_successors(self):
        if len(list(self._mpg.successors(self))) > 0:
            return True
        else:
            return False

    def get_successors(self):
        return list(self._mpg.successors(self))

    def has_predecessors(self):
        if len(list(self._mpg.predecessors(self))) > 0:
            return True
        else:
            return False

    def get_predecessors(self):
        return list(self._mpg.predecessors(self))

    def coords_overlap(self, other_node):
        my_lend, my_rend = self.get_coords()
        other_lend, other_rend = other_node.get_coords()

        if my_lend <= other_rend and my_rend >= other_lend:
            return True
        else:
            return False

    def add_containment(self, other_node):
        self._containments.add(other_node)

    def has_containments(self):
        if len(self._containments) > 0:
            return True
        else:
            return False

    def get_containments(self):
        return list(self._containments)

    def contains_other_node(self, other_node):

        self_path = self.get_simple_path()
        other_path = other_node.get_simple_path()

        # this is spacer-aware
        if Simple_path_utils.simple_path_A_contains_and_compatible_with_simple_path_B_spacefree_region_path_A(
            self.get_splice_graph(), self_path, other_path
        ):
            # logger.debug("Path {}\ncontains Path {}".format(self.get_simple_path(), other_node.get_simple_path()))
            assert len(self_path) >= len(
                other_path
            ), "Error, len({}) not >= len({})".format(self_path, other_path)
            return True
        else:
            return False

    def compatible(self, other_node):
        if self._multiPath.is_overlapping_and_compatible(other_node._multiPath):
            return True
        else:
            return False

    def get_mpgn_pair_token(mpgn_A, mpgn_B):
        assert type(mpgn_A) == MultiPathGraphNode
        assert type(mpgn_B) == MultiPathGraphNode

        token = "^^".join(sorted((mpgn_A._id, mpgn_B._id)))
        return token

    def get_splice_graph(self):
        return self.get_multipath_graph().get_splice_graph()

    def get_splice_graph_node_objs_for_path(self):
        simple_path = self.get_simple_path()
        splice_graph_node_objs = list()

        sg = self.get_splice_graph()

        for node_id in simple_path:
            if node_id == SPACER:
                splice_graph_node_objs.append(None)
            else:
                splice_graph_node_obj = sg.get_node_obj_via_id(node_id)
                splice_graph_node_objs.append(splice_graph_node_obj)

        return splice_graph_node_objs

    def _compute_seq_length(self):
        sg_nodes = self.get_splice_graph_node_objs_for_path()

        seq_length = 0
        for sg_node in sg_nodes:
            if sg_node is not None:
                lend, rend = sg_node.get_coords()
                feature_len = rend - lend + 1
                if type(sg_node) == Exon:
                    seq_length += feature_len
                # elif type(sg_node) == Intron:
                #    seq_length += min(feature_len, 10000)  #//FIXME: determine max allowable length based on intron length distribution.

        return seq_length

    def reevaluate_weighting_via_path_compatibilities(self, transcript_multiPath):

        logger.debug("reevaluating weights for {}".format(self))
        # sys.stderr.write("reevaluating weights for {}".format(self))

        assert type(transcript_multiPath) == MultiPath.MultiPath

        self.set_weight(0.001)

        """
        transcript_simple_path = transcript_multiPath.get_simple_path()
        
        compatible_score = self.get_score_INCLUDE_containments(use_prev_weight=True)
        incompatible_score = 0

        connected_mpgns = self.get_predecessors() + self.get_successors()

        for node in connected_mpgns:
            node_simple_path = node.get_simple_path()
            if Simple_path_utils.path_A_contains_path_B(transcript_simple_path, node_simple_path):
                compatible_score += node.get_score_INCLUDE_containments(use_prev_weight=True)
            else:
                incompatible_score += node.get_score_INCLUDE_containments(use_prev_weight=True)

        pseudocount = 1e-3 
        
        fraction_compatible = (compatible_score + pseudocount) / (compatible_score + incompatible_score + pseudocount)
        
        current_weight = self.get_weight()

        adjusted_weight = current_weight - (fraction_compatible * current_weight)

        adjusted_weight = max(MIN_WEIGHT, adjusted_weight) # to avoid going too low.
        
        self.set_weight(adjusted_weight)

        logger.debug("fraction compatibile: {}, adjusted weight -> {}".format(fraction_compatible, adjusted_weight))
        """

        return
