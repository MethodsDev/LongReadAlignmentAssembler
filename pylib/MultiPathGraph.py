import sys, os, re
import Splice_graph
import MultiPath
import MultiPathCounter
import networkx as nx
import Simple_path_utils as Simple_path_utils
import LRAA_Globals
from LRAA_Globals import SPACER
from GenomeFeature import *
from MultiPathGraphNode import MultiPathGraphNode
from collections import defaultdict

import logging

logger = logging.getLogger(__name__)

ITER = 0


class MultiPathGraph:

    def __init__(
        self,
        multiPathCounter,
        splice_graph,
        contig_acc,
        contig_strand,
        min_mpgn_read_count=1,
        allow_spacers=False,
        read_name_to_span=None,
    ):

        local_debug = False

        logger.info(f"START building MultiPathGraph for {contig_acc}")

        assert type(multiPathCounter) == MultiPathCounter.MultiPathCounter
        assert type(splice_graph) == Splice_graph.Splice_graph

        self._splice_graph = splice_graph
        self._contig_acc = contig_acc
        self._contig_strand = contig_strand

        mp_graph = nx.DiGraph()
        self._mp_graph = mp_graph
        self._read_name_to_span = (
            read_name_to_span if read_name_to_span is not None else {}
        )

        self._mp_graph_nodes_list = list()

        self._incompatible_pairs = set()  # store sets of pairs of incompatible nodes.

        self._mp_id_to_node = dict()  # mp_id -> mpg_node

        sg_component_to_mp_id = defaultdict(set)

        multiPathCountPairs = multiPathCounter.get_all_MultiPathCountPairs()
        for mpCountPair in multiPathCountPairs:
            mp, count = mpCountPair.get_multipath_and_count()

            if count < min_mpgn_read_count:
                continue

            # print("mp: {}, count {}".format(mp, count))

            first_node_id = mp[0]
            last_node_id = mp[-1]

            first_node_obj = splice_graph.get_node_obj_via_id(first_node_id)
            last_node_obj = splice_graph.get_node_obj_via_id(last_node_id)

            lend_coord = first_node_obj.get_coords()[0]
            rend_coord = last_node_obj.get_coords()[1]

            mp_graph_node = MultiPathGraphNode(
                mp, count, lend_coord, rend_coord, mpg=self
            )
            mp_graph.add_node(mp_graph_node)

            # Optional: set read span bounds from supporting reads
            if self._read_name_to_span:
                all_reads = mp_graph_node.get_read_names()
                spans = [self._read_name_to_span[r] for r in all_reads if r in self._read_name_to_span]
                if spans:
                    min_l = min(x[0] for x in spans)
                    max_r = max(x[1] for x in spans)
                    mp_graph_node.set_read_span_bounds(min_l, max_r)

            if type(first_node_obj) in (TSS, PolyAsite):
                mp_graph_node._left_boundary = 1

            if type(last_node_obj) in (TSS, PolyAsite):
                mp_graph_node._right_boundary = 1

            self._mp_graph_nodes_list.append(mp_graph_node)

            # assign mp to splice graph component
            mp_graph_node_id = mp_graph_node.get_id()
            self._mp_id_to_node[mp_graph_node_id] = mp_graph_node

            component_id = self._splice_graph._node_id_to_component[first_node_id]
            # print(f"{mp_graph_node_id} first node is {first_node_id} and assigned to component_id {component_id}")
            sg_component_to_mp_id[component_id].add(mp_graph_node_id)

        ## sort
        self._mp_graph_nodes_list = sorted(
            self._mp_graph_nodes_list,
            key=lambda x: (
                x._lend,
                x._rend,
                x.get_left_boundary_sort_weight(),
                x.get_right_boundary_sort_weight(),
            ),
        )

        if LRAA_Globals.DEBUG:
            mpg_build_dir = f"__mpg_building/{self._contig_acc}"
            if not os.path.exists(mpg_build_dir):
                os.makedirs(mpg_build_dir)
            # build_file = os.path.join(mpg_build_dir, "build-{}.txt".format(self._contig_acc))
            # build_ofh = open(build_file, "wt")

        if LRAA_Globals.DEBUG:
            global ITER
            ITER += 1
            component_descr_file = f"__MPGN_components_described.{ITER}.bed"
            if os.path.exists(component_descr_file):
                raise RuntimeError(
                    "Error! - delete existing file {}".format(component_descr_file)
                )
            component_descr_ofh = open(component_descr_file, "a")

        sorted_component_ids = sorted(
            sg_component_to_mp_id.keys(),
            key=lambda x: len(sg_component_to_mp_id[x]),
            reverse=True,
        )
        for component_id in sorted_component_ids:
            mp_node_set = sg_component_to_mp_id[component_id]
            # print(mp_node_set)
            num_paths = len(mp_node_set)
            logger.info(f"Component {component_id} has {num_paths} paths assigned.")
            if LRAA_Globals.DEBUG:
                component_description_file = os.path.join(
                    mpg_build_dir, f"{component_id}.comp.descr.tsv"
                )
                mp_nodes = [
                    self._mp_id_to_node[mp_node_id] for mp_node_id in mp_node_set
                ]
                mp_nodes = sorted(mp_nodes, key=lambda x: x._lend)
                with open(component_description_file, "a") as ofh:
                    for mp_node in mp_nodes:
                        print(
                            "\t".join(
                                [
                                    contig_acc,
                                    str(mp_node._lend),
                                    str(mp_node._rend),
                                    str(component_id),
                                    str(mp_node.get_count()),
                                    contig_strand,
                                    str(mp_node.get_simple_path()),
                                    ",".join(mp_node.get_read_names()),
                                ]
                            ),
                            file=ofh,
                        )

                        print(
                            "\t".join(
                                [
                                    contig_acc,
                                    str(mp_node._lend),
                                    str(mp_node._rend),
                                    "Comp:"
                                    + str(component_id)
                                    + ", count: "
                                    + str(mp_node.get_count()),
                                    contig_strand,
                                    str(mp_node.get_simple_path()),
                                ]
                            ),
                            file=component_descr_ofh,
                        )

        for component_id in sorted_component_ids:
            mp_node_set = sg_component_to_mp_id[component_id]

            ordered_nodes = list()
            for mp_node_id in mp_node_set:
                mp_node_obj = self._mp_id_to_node[mp_node_id]
                ordered_nodes.append(mp_node_obj)
                # print(f"component_id: {component_id}\tmp_node_obj: {mp_node_obj}")

            max_nodes = LRAA_Globals.config["max_path_nodes_per_component"]
            if len(ordered_nodes) > max_nodes:
                logger.info(
                    "Size of component node set is too large... shrinking to {}".format(
                        max_nodes
                    )
                )
                ordered_nodes = sorted(
                    ordered_nodes,
                    key=lambda x: (x.get_count(), x._seq_length),
                    reverse=True,
                )

                # prune the excess:
                nodes_to_prune = ordered_nodes[max_nodes:]
                # retain max nodes
                ordered_nodes = ordered_nodes[0:max_nodes]

                self._prune_nodes(nodes_to_prune)
                sg_component_to_mp_id[component_id] = ordered_nodes

            ordered_nodes = sorted(
                ordered_nodes,
                key=lambda x: (
                    x._lend,
                    x._rend,
                    x.get_left_boundary_sort_weight(),
                    x.get_right_boundary_sort_weight(),
                ),
            )

            num_ordered_nodes = len(ordered_nodes)
            region = "{}{}:{}-{}".format(
                contig_acc,
                contig_strand,
                ordered_nodes[0]._lend,
                ordered_nodes[-1]._rend,
            )
            logger.info(
                f"Building MP Graph for component_id {component_id} spanning {region} with {num_ordered_nodes} multipaths"
            )

            if LRAA_Globals.DEBUG and local_debug:
                component_build_file = os.path.join(
                    mpg_build_dir, f"{component_id}.comp.buildgraph.tsv"
                )
                build_ofh = open(component_build_file, "a")

            ## define edges, containments, and incompatibilities
            for i in range(0, len(ordered_nodes)):
                node_i = ordered_nodes[i]

                for j in range(i - 1, -1, -1):
                    node_j = ordered_nodes[j]

                    if LRAA_Globals.DEBUG and local_debug:
                        print(
                            "\n\n# comparing prev_j\n{}\nto_i\n{}".format(
                                node_j, node_i
                            ),
                            file=build_ofh,
                        )

                    # nope - need more clever logic tracking prev max rend in ordered list.
                    # if node_j._rend < node_i._lend:
                    #    if LRAA_Globals.DEBUG:
                    #        print("-non-overlapping, short-circuiting", file=build_ofh)
                    #    break # all earlier node j's will also be non-overlapping

                    if node_j._rend < node_i._lend:
                        # they do not overlap and so cannot be contained or overlapping/incompatible
                        continue

                    if node_i.contains_other_node(node_j):
                        # i contains j
                        if LRAA_Globals.DEBUG and local_debug:
                            print("i-contains-j", file=build_ofh)
                        node_i.add_containment(node_j)

                    elif node_j.contains_other_node(node_i):
                        # j contains i
                        if LRAA_Globals.DEBUG and local_debug:
                            print("j-contains-i", file=build_ofh)
                        node_j.add_containment(node_i)

                    elif node_i.compatible(node_j):
                        # draw edge between overlapping and compatible nodes.
                        if LRAA_Globals.DEBUG and local_debug:
                            print("i-COMPATIBLE-j", file=build_ofh)
                            # logger.debug("adding edge: {},{}".format(node_j, node_i))

                        if not LRAA_Globals.config["restrict_asm_to_collapse"]:
                            self._mp_graph.add_edge(node_j, node_i)

                    else:
                        # incompatible pairs
                        if LRAA_Globals.DEBUG and local_debug:
                            print("i-NOTcompatible-j", file=build_ofh)
                        incompatible_pair_token = (
                            MultiPathGraphNode.get_mpgn_pair_token(node_i, node_j)
                        )
                        self._incompatible_pairs.add(incompatible_pair_token)

        if LRAA_Globals.DEBUG:
            component_descr_ofh.close()

        logger.info(f"DONE building MultiPathGraph for {contig_acc}")

        return

    def get_ordered_nodes(self):
        # these are sorted by rend
        return list(self._mp_graph_nodes_list)

    def init_mpgn_reweighting_flags(self):
        mpgn_nodes = self.get_ordered_nodes()
        for mpgn in mpgn_nodes:
            mpgn.set_reweighted_flag(False)
        return

    def has_edge(self, multiPath_before, multiPath_after):
        return self._mp_graph.has_edge(multiPath_before, multiPath_after)

    def successors(self, mpgn):
        return self._mp_graph.successors(mpgn)

    def predecessors(self, mpgn):
        return self._mp_graph.predecessors(mpgn)

    def incompatible_mpgn_pair(self, mpgn_A, mpgn_B):
        mpgn_pair_token = MultiPathGraphNode.get_mpgn_pair_token(mpgn_A, mpgn_B)
        if mpgn_pair_token in self._incompatible_pairs:
            return True
        else:
            return False

    def get_splice_graph(self):
        return self._splice_graph

    def define_disjoint_graph_components_via_graph_traversal(self):

        mpgn_list = self.get_ordered_nodes()
        mpgn_seen = set()

        component_list = list()

        while len(mpgn_list) > 0:

            queue = list()
            seed_node = mpgn_list.pop(0)

            if seed_node in mpgn_seen:
                continue

            # start a new component.
            component = list()
            queue.append(seed_node)

            while len(queue) > 0:
                node = queue.pop(0)
                if node not in mpgn_seen:
                    mpgn_seen.add(node)
                    component.append(node)
                    # add predecessors and successors to queue
                    if node.has_predecessors():
                        queue.extend(node.get_predecessors())
                    if node.has_successors():
                        queue.extend(node.get_successors())
                    if node.has_containments():
                        queue.extend(node.get_containments())

            if len(component) > 0:
                component_list.append(component)

        logger.info(
            "identified {} disjoint graph components".format(len(component_list))
        )

        ## assign the component ids
        component_counter = 0
        for component in component_list:
            component_counter += 1
            # print("component: {}, counter: {}".format(component, component_counter))
            for mpgn in component:
                mpgn.set_component_id(component_counter)

        return component_list

    def define_disjoint_graph_components_via_shared_splice_graph_vertex(self):

        all_splice_graph_nodes = set()
        splice_graph_node_to_mpgns = defaultdict(set)

        # populate splice graph node to mpgn data structures
        mpgn_list = self.get_ordered_nodes()
        for mpgn in mpgn_list:
            splice_graph_nodes = mpgn.get_splice_graph_node_objs_for_path()
            for splice_graph_node in splice_graph_nodes:
                all_splice_graph_nodes.add(splice_graph_node)
                splice_graph_node_to_mpgns[splice_graph_node].add(mpgn)

        # build disjoint components
        splice_graph_nodes_visited = set()

        component_list = list()
        all_splice_graph_nodes = list(
            all_splice_graph_nodes
        )  # convert earlier set to list
        while len(all_splice_graph_nodes) > 0:

            queue = list()
            seed_splice_graph_node = all_splice_graph_nodes.pop()

            if seed_splice_graph_node in splice_graph_nodes_visited:
                continue

            # start a new component
            component = list()
            queue.append(seed_splice_graph_node)

            while len(queue) > 0:
                splice_graph_node = queue.pop(0)
                mpgns = list(splice_graph_node_to_mpgns[splice_graph_node])
                splice_graph_nodes_visited.add(splice_graph_node)
                for mpgn in mpgns:
                    if mpgn not in component:
                        component.append(mpgn)

                    mpgn_splice_graph_nodes = mpgn.get_splice_graph_node_objs_for_path()
                    for splice_graph_node in mpgn_splice_graph_nodes:
                        if (splice_graph_node not in splice_graph_nodes_visited) and (
                            splice_graph_node not in queue
                        ):
                            queue.append(splice_graph_node)

            component_list.append(component)

        logger.info(
            "identified {} disjoint graph components".format(len(component_list))
        )

        ## assign the component ids
        component_counter = 0
        for component in component_list:
            component_counter += 1
            # print("component: {}, counter: {}".format(component, component_counter))
            for mpgn in component:
                mpgn.set_component_id(component_counter)

        return component_list

    def describe_graph(self, output_filename):

        if True:
            return  # files can be too large, disabling for now

        ofh = open(output_filename, "a")

        ofh2 = open(output_filename + ".tsv", "a")  # tabulated output

        mpgn_list = self.get_ordered_nodes()
        for mpgn in mpgn_list:
            ofh.write(str(mpgn) + "\n")
            all_mpgn_read_names = mpgn.get_read_names()
            all_contained_read_names = set()

            for containment in mpgn.get_containments():
                ofh.write("\t" + str(containment) + "\n")
                containment_read_names = containment.get_read_names()
                all_contained_read_names.update(containment_read_names)
            ofh.write("\n")

            read_names_not_from_containments = (
                all_mpgn_read_names - all_contained_read_names
            )
            mpgn_id = mpgn.get_id()
            for read_name in read_names_not_from_containments:
                print("\t".join([mpgn_id, "primary", read_name]), file=ofh2)
            for read_name in all_contained_read_names:
                print("\t".join([mpgn_id, "contained", read_name]), file=ofh2)

        ofh.close()
        ofh2.close()

        return

    def write_mp_graph_nodes_to_gtf(self, gtf_output_filename):

        contig_acc = self._splice_graph.get_contig_acc()

        ofh = open(gtf_output_filename, "a")

        mpgn_list = self.get_ordered_nodes()
        for mpgn in mpgn_list:
            splice_graph_nodes = mpgn.get_splice_graph_node_objs_for_path()

            mpgn_read_count = mpgn.get_count()
            mpgn_component_id = mpgn.get_component_id()

            mpgn_id = mpgn.get_id()
            trans_id = (
                "t__count={}_Comp={}_.".format(mpgn_read_count, mpgn_component_id)
                + mpgn_id
            )
            gene_id = (
                "g__count={}_Comp={}_".format(mpgn_read_count, mpgn_component_id)
                + mpgn_id
            )

            for splice_graph_node in splice_graph_nodes:
                if splice_graph_node is not None and type(splice_graph_node) == Exon:

                    coords = splice_graph_node.get_coords()

                    ofh.write(
                        "\t".join(
                            [
                                contig_acc,
                                "MPGN",
                                "exon",
                                str(coords[0]),
                                str(coords[1]),
                                ".",
                                "?",
                                ".",
                                'gene_id "{}"; transcript_id "{}";'.format(
                                    gene_id, trans_id
                                ),
                            ]
                        )
                        + "\n"
                    )

        ofh.close()

        return

    def remove_small_components(self, mpg_components, min_transcript_length):

        surviving_components = list()

        for mpgn_list in mpg_components:
            max_seq_len = 0
            for mpgn in mpgn_list:
                max_seq_len += mpgn.get_seq_length()
            if max_seq_len < min_transcript_length:
                # component is too small to generate a sufficiently large transcript
                self._mp_graph.remove_nodes_from(mpgn_list)
            else:
                surviving_components.append(mpgn_list)

        return surviving_components

    def _prune_nodes(self, nodes_to_prune):

        logger.info(
            "-pruning {} nodes to reduce component size".format(len(nodes_to_prune))
        )
        for node in nodes_to_prune:
            del self._mp_id_to_node[node.get_id()]
            self._mp_graph.remove_node(node)
            self._mp_graph_nodes_list.remove(node)

        return
