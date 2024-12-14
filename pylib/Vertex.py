import sys, os, re
from MultiPath import MultiPath
from MultiPathGraph import MultiPathGraphNode
from Scored_path import Scored_path
import LRAA_Globals

import logging

logger = logging.getLogger(__name__)


class Vertex:

    def __init__(self, multipath_graph_node):

        assert type(multipath_graph_node) == MultiPathGraphNode

        self._multipath_graph_node = multipath_graph_node

        self._fromPaths = list()  # hold scored paths.

        # add unextended current path node as initial path
        initial_scored_path = Scored_path([self._multipath_graph_node])
        self._fromPaths.append(initial_scored_path)

        return

    def __repr__(self):
        return "Vertex for {}".format(self._multipath_graph_node.get_simple_path())

    def get_multipath_graph_node(self):
        return self._multipath_graph_node

    def get_mpgn(self):
        return self.get_multipath_graph_node()

    def get_multipath_graph(self):
        return self._multipath_graph_node.get_multipath_graph()

    def get_fromPaths(self):
        return list(self._fromPaths)

    def add_highest_scoring_path_extension(self, prev_pasa_vertex):

        best_score = 0
        best_prev_scored_path = None

        self_mpgn = self.get_mpgn()
        assert type(self_mpgn) == MultiPathGraphNode

        if LRAA_Globals.DEBUG:
            path_ext_dir = "__path_extension_audits"
            if not os.path.exists(path_ext_dir):
                os.makedirs(path_ext_dir)
            mp_id = self_mpgn.get_id()
            path_ext_audit_file = os.path.join(
                path_ext_dir, "{}.extension_audit".format(mp_id)
            )
            extension_audit_ofh = open(path_ext_audit_file, "at")
            print(
                "############################################################################",
                file=extension_audit_ofh,
            )
            print(
                "# Evaluating comparison to prev pasa vertex: {} with from paths:".format(
                    prev_pasa_vertex
                ),
                file=extension_audit_ofh,
            )
            for i, from_path in enumerate(prev_pasa_vertex.get_fromPaths()):
                print(
                    "[From path {}: {}".format(i, from_path), file=extension_audit_ofh
                )
            print("##### Done from path listing.", file=extension_audit_ofh)

        for prev_scored_path in prev_pasa_vertex.get_fromPaths():
            if prev_scored_path.incompatibility_detected(self_mpgn):

                if LRAA_Globals.DEBUG:
                    print(
                        "-incompatible extension of {}\nto prev scored path: {}".format(
                            self_mpgn, prev_scored_path
                        ),
                        file=extension_audit_ofh,
                    )

            else:

                extension_path_candidate = (
                    prev_scored_path.create_scored_path_extension(self_mpgn)
                )

                if LRAA_Globals.DEBUG:
                    print(
                        "-extension path candidate: {}".format(
                            extension_path_candidate
                        ),
                        file=extension_audit_ofh,
                    )

                if extension_path_candidate.get_score() > best_score:
                    best_score = extension_path_candidate.get_score()
                    best_prev_scored_path = extension_path_candidate

                    if LRAA_Globals.DEBUG:
                        print("\t** best so far.", file=extension_audit_ofh)

        if best_prev_scored_path is not None:
            self._fromPaths.append(best_prev_scored_path)
            logger.debug(
                "Added best extension path: {} to {}".format(
                    best_prev_scored_path.get_simple_path(), self
                )
            )

            if LRAA_Globals.DEBUG:
                print(
                    "**  best selected extension path: {}".format(
                        best_prev_scored_path
                    ),
                    file=extension_audit_ofh,
                )

        return

    def rescore_fromPaths(self):

        for scored_path in self.get_fromPaths():
            scored_path.rescore()

        return

    def describe_pasa_vertex(self):

        pv_text = str(self) + "\n"

        from_paths = self.get_fromPaths()

        from_paths = sorted(from_paths, key=lambda x: x._score, reverse=True)

        from_path_counter = 0
        for path in from_paths:
            from_path_counter += 1
            pv_text += "\n\t[fromPath:{}] ".format(from_path_counter) + str(path) + "\n"

        return pv_text
