#!/usr/bin/env python
# encoding: utf-8

import sys, os, re
import hashlib
from collections import defaultdict
from LRAA_Globals import SPACER
import Simple_path_utils
from Util_funcs import coordpairs_overlap
import logging
from GenomeFeature import Exon, Intron
from unittest.mock import Mock
from Splice_graph import Splice_graph
import Transcript

logger = logging.getLogger(__name__)


##
## The MultiPath is a list of node IDs that can be interrupted by SPACERs.
## The MultiPath stores the splice graph so that the nodes corresponding to those IDs can be retrieved as needed.
##


class MultiPath:
        # Memory model for read tracking:
        # - In RAM: store compact integer read identifiers in _read_ids.
        # - On disk: ReadNameStore (SQLite/LMDB) maps long read-name strings <-> integer IDs.
        # - If _read_ids is empty, get_read_ids() can stream IDs from MP_READ_ID_STORE.
        # Prefer using IDs for algorithmic decisions; resolve names only for reporting.

    MP_id_counter = 0

    def __init__(self, splice_graph, paths_list, read_types=set(), read_names=set(), read_count=None):

        self._splice_graph = splice_graph

        self._simple_path = self._merge_paths_to_simple_multi_path(paths_list)
        self._simple_path_trimmed_tuple = tuple(
            Simple_path_utils.trim_terminal_spacers(self._simple_path.copy())
        )

        # print(str(self._simple_path))

        MultiPath.MP_id_counter += 1

        self._id = "MP{}".format(MultiPath.MP_id_counter)

        # determine span
        coords = list()
        exon_segments = list()
        for node_id in self._simple_path:
            if node_id != SPACER:
                node_obj = self._splice_graph.get_node_obj_via_id(node_id)
                coordpair = node_obj.get_coords()
                coords.extend(coordpair)
                if type(node_obj) == Exon:
                    exon_segments.append(coordpair)

        coords = sorted(coords)

        self._lend = coords[0]
        self._rend = coords[-1]

        self._read_types = read_types
        # Store compact read identifiers (ints) in-memory; resolve to long names only on demand
        self._read_ids = set()
        if read_names is not None:
            for rn in read_names:
                rid = self._coerce_read_identifier(rn)
                if rid is not None:
                    self._read_ids.add(int(rid))
        # maintain count independent of whether names are stored
        if read_count is not None:
            self._read_count = int(read_count)
        else:
            self._read_count = len(self._read_ids)
            if self._read_count == 0 and read_names is not None:
                try:
                    self._read_count = len(read_names)
                except Exception:
                    self._read_count = 0

        self._exon_segments = Simple_path_utils.merge_adjacent_segments(exon_segments)

        return

    def get_id(self):
        return self._id

    def get_simple_path(self):
        if self._simple_path_trimmed_tuple is None:
            self._simple_path_trimmed_tuple = tuple(
                Simple_path_utils.trim_terminal_spacers(self._simple_path.copy())
            )
        return list(self._simple_path_trimmed_tuple)

    def get_simple_path_tuple(self):
        if self._simple_path_trimmed_tuple is None:
            self._simple_path_trimmed_tuple = tuple(
                Simple_path_utils.trim_terminal_spacers(self._simple_path.copy())
            )
        return self._simple_path_trimmed_tuple

    def get_exon_segments(self):
        return self._exon_segments.copy()

    # read type handling

    def include_read_type(self, read_type):
        if type(read_type) in [list, set]:
            for r in read_type:
                self._read_types.add(r)
        else:
            self._read_types.add(read_type)

    def has_read_type(self, read_type):
        bool_has_read_type = read_type in self._read_types
        return bool_has_read_type

    def get_read_types(self):
        return self._read_types.copy()

    # read name handling

    def get_read_names(self):
        """
        Returns the set of read names supporting this multipath.

        Behavior:
        - MultiPath stores compact read IDs in-memory. This method resolves those IDs to long string names
          using LRAA_Globals.READ_NAME_STORE as needed.
        - If no IDs are retained in-memory, attempts to stream read IDs from LRAA_Globals.MP_READ_ID_STORE and
          resolve to names via READ_NAME_STORE.

        Notes:
        - Resolving names can be relatively expensive. Prefer get_read_count() or get_read_ids() when possible.
        """
        ids = set(self._read_ids)
        if not ids:
            # try to stream IDs from external store if configured
            try:
                import LRAA_Globals  # local import to avoid cycles
                mp_store = getattr(LRAA_Globals, "MP_READ_ID_STORE", None)
                if mp_store is not None:
                    ids = set(int(rid) for rid in mp_store.iter_read_ids(self.get_id()))
            except Exception:
                ids = set()
        # resolve ids -> names
        try:
            import LRAA_Globals
            name_store = getattr(LRAA_Globals, "READ_NAME_STORE", None)
            if name_store is not None:
                names = set()
                for rid in ids:
                    nm = name_store.get_name(int(rid))
                    if nm is not None:
                        names.add(nm)
                return names
        except Exception:
            pass
        # Fallback if name store not available: return stringified IDs
        return set(str(rid) for rid in ids)

    def get_read_ids(self):
        """
        Returns the set of compact read IDs supporting this multipath. If not retained in-memory,
        will attempt to stream from MP_READ_ID_STORE when available.
        """
        if len(self._read_ids) > 0:
            return set(self._read_ids)
        try:
            import LRAA_Globals
            mp_store = getattr(LRAA_Globals, "MP_READ_ID_STORE", None)
            if mp_store is not None:
                return set(int(rid) for rid in mp_store.iter_read_ids(self.get_id()))
        except Exception:
            pass
        return set()

    def get_read_count(self):
        return self._read_count

    def get_read_names_count(self):
        # count of IDs stored or externally associated
        if len(self._read_ids) > 0:
            return len(self._read_ids)
        try:
            import LRAA_Globals
            mp_store = getattr(LRAA_Globals, "MP_READ_ID_STORE", None)
            if mp_store is not None:
                return mp_store.count(self.get_id())
        except Exception:
            pass
        return 0

    # Alias reflecting ID-based semantics
    def get_read_id_count(self):
        return self.get_read_names_count()

    def include_read_name(self, read_name):
        # Accept a single read identifier or a collection; coerce to ID when possible
        if type(read_name) in [list, set]:
            added = 0
            for r in read_name:
                rid = self._coerce_read_identifier(r)
                if rid is not None:
                    self._read_ids.add(int(rid))
                    added += 1
            self._read_count += added if added > 0 else len(read_name)
        else:
            rid = self._coerce_read_identifier(read_name)
            if rid is not None:
                self._read_ids.add(int(rid))
            self._read_count += 1

    # New preferred alias for clarity
    def include_read_id(self, read_id):
        if type(read_id) in [list, set]:
            added = 0
            for r in read_id:
                try:
                    self._read_ids.add(int(r))
                    added += 1
                except Exception:
                    continue
            self._read_count += added if added > 0 else len(read_id)
        else:
            try:
                self._read_ids.add(int(read_id))
            except Exception:
                pass
            self._read_count += 1

    def include_read_count(self, increment=1):
        # count-only change; IDs unchanged
        self._read_count += int(increment)

    def remove_read_name(self, read_name):
        rid = self._coerce_read_identifier(read_name)
        if rid is not None and int(rid) in self._read_ids:
            self._read_ids.discard(int(rid))
            self._read_count = max(0, self._read_count - 1)

    # Alias for clarity
    def remove_read_id(self, read_id):
        try:
            rid = int(read_id)
            if rid in self._read_ids:
                self._read_ids.discard(rid)
                self._read_count = max(0, self._read_count - 1)
        except Exception:
            pass

    def prune_reftranscript_as_evidence(self):
        # Remove any reads whose resolved names start with "reftranscript:"
        try:
            import LRAA_Globals
            name_store = getattr(LRAA_Globals, "READ_NAME_STORE", None)
        except Exception:
            name_store = None

        if name_store is None:
            # Without a name store we cannot resolve, so nothing to prune safely
            return

        ids_to_keep = set()
        pruned = 0
        # Prefer in-memory IDs if present; otherwise try streaming IDs
        ids = set(self._read_ids)
        if not ids:
            try:
                import LRAA_Globals
                mp_store = getattr(LRAA_Globals, "MP_READ_ID_STORE", None)
                if mp_store is not None:
                    ids = set(int(rid) for rid in mp_store.iter_read_ids(self.get_id()))
            except Exception:
                ids = set()

        for rid in ids:
            nm = name_store.get_name(int(rid))
            if nm is None or not nm.startswith("reftranscript:"):
                ids_to_keep.add(int(rid))
            else:
                pruned += 1

        # Replace internal set with kept IDs and adjust count
        self._read_ids = ids_to_keep
        if pruned > 0:
            self._read_count = max(0, self._read_count - pruned)

    # splice graph operations

    def get_splice_graph(self):
        return self._splice_graph

    def get_ordered_exons_and_introns(self):
        simple_path = self.get_simple_path()

        sg = self.get_splice_graph()

        # spacers become None

        exons_and_introns = list()
        for node_id in simple_path:
            if node_id == SPACER:
                exons_and_introns.append(None)
            else:
                obj = sg.get_node_obj_via_id(node_id)
                exons_and_introns.append(obj)

        return exons_and_introns

    def get_introns(self):
        simple_path = self.get_simple_path()

        sg = self.get_splice_graph()

        introns = list()
        for node_id in simple_path:
            if node_id != SPACER:
                obj = sg.get_node_obj_via_id(node_id)
                if type(obj) == Intron:
                    introns.append(obj)

        return introns

    def __len__(self):
        return len(self._simple_path)

    def __getitem__(self, i):
        length = len(self)
        if i < 0:
            i += length
        if 0 <= i < length:
            return self._simple_path[i]

        raise IndexError("Index out of range: {}".format(i))

    def get_coords(self):
        return (self._lend, self._rend)

    def get_cdna_length(self):

        cdna_len = 0

        exons_and_introns = self.get_ordered_exons_and_introns()
        for feature in exons_and_introns:
            if type(feature) == Exon:
                lend, rend = feature.get_coords()
                exon_len = rend - lend + 1
                cdna_len += exon_len

        return cdna_len

    def _merge_paths_to_simple_multi_path(self, paths_list):

        paths_to_asm = [path for path in paths_list]  # copy incoming list

        sg = self._splice_graph

        def sort_func(simple_path):
            return self._splice_graph.get_node_obj_via_id(simple_path[0]).get_coords()[
                0
            ]

        paths_to_asm = sorted(paths_to_asm, key=sort_func)

        ## perform cycles of merge attempts, retaining relative ordering.
        assembled_paths = list()
        # seed it with the first entry
        assembled_paths.append(paths_to_asm.pop(0))

        while True:

            unmerged_paths = list()

            seed_asm = assembled_paths[-1]

            merged_flag = False

            for other_path in paths_to_asm:
                # check for merge
                if Simple_path_utils.are_overlapping_and_compatible_NO_gaps_in_overlap(
                    seed_asm, other_path
                ):
                    # mergeable, so lets merge them
                    merged_asm = Simple_path_utils.merge_simple_paths(
                        seed_asm, other_path
                    )
                    # update seed asm
                    seed_asm = assembled_paths[-1] = merged_asm
                    merged_flag = True

                elif Simple_path_utils.simple_paths_overlap_and_compatible_spacer_aware_both_paths(
                    sg, seed_asm, other_path
                ):
                    seed_asm = assembled_paths[-1] = (
                        Simple_path_utils.merge_simple_paths_containing_spacers(
                            sg, seed_asm, other_path
                        )
                    )
                    merged_flag = True
                elif Simple_path_utils.simple_paths_overlap(sg, seed_asm, other_path):
                    # logger.warning("-warning: multipath subpaths overlap but contain spacers or are incompatible")
                    # logger.warning("A: {}\nB:{}".format(seed_asm, other_path))
                    # use the longer one.
                    seed_asm = assembled_paths[-1] = (
                        seed_asm if len(seed_asm) > len(other_path) else other_path
                    )
                    # //FIXME: should use better criteria than this above.

                else:
                    unmerged_paths.append(other_path)

            if unmerged_paths:
                if not merged_flag:
                    # must make a new seed.
                    assembled_paths.append(unmerged_paths.pop(0))

                paths_to_asm = unmerged_paths

            else:
                break  # done assembling

        ## build multipath for post-assembly products
        simple_multipath = []
        for i, path in enumerate(assembled_paths):
            simple_multipath += path
            if i != len(assembled_paths) - 1:
                simple_multipath.append(SPACER)

        simple_multipath = Simple_path_utils.trim_terminal_spacers(simple_multipath)

        if SPACER in simple_multipath:
            simple_multipath = Simple_path_utils.try_fill_spacers_via_splicegraph(
                sg, simple_multipath
            )

        return simple_multipath

    def is_overlapping_contained_and_compatible(self, other_multipath):

        my_lend, my_rend = self.get_coords()
        other_lend, other_rend = other_multipath.get_coords()
        if my_lend <= other_lend and my_rend >= other_rend:
            # first check containment
            return self.is_overlapping_and_compatible(other_multipath)
        else:
            return False

    def is_overlapping_and_compatible(self, other_multipath):
        # overlapping parts are required to be identical
        # compatible means no conflict detected.
        # spacer-aware

        assert type(other_multipath) == MultiPath

        if not coordpairs_overlap(self.get_coords(), other_multipath.get_coords()):
            return False

        my_path = self.get_simple_path()
        other_path = other_multipath.get_simple_path()

        if my_path == other_path:
            return True
        else:
            return Simple_path_utils.simple_paths_overlap_and_compatible_spacefree_region_path_A(
                self.get_splice_graph(), my_path, other_path
            )

    def split_multipath_at_spacers(self):
        simple_path = self.get_simple_path()
        if SPACER not in simple_path:
            return [self]

        # need to split at spacers
        logger.info("attempt split path at spacers: {}".format(simple_path))
        split_simple_paths_list = Simple_path_utils.split_path_at_spacers(simple_path)

        split_mps = list()
        for split_simple_path in split_simple_paths_list:
            split_mp = MultiPath(
                self._splice_graph,
                [split_simple_path],
                read_types=self.get_read_types(),
                read_names=self.get_read_ids(),
            )

            split_mps.append(split_mp)

        logger.debug("-SPACER SPLITTING of {}\ninto\n{}".format(self, split_mps))

        return split_mps

    def toTranscript(self):

        sg = self.get_splice_graph()
        contig_acc = sg.get_contig_acc()
        contig_strand = sg.get_contig_strand()

        exon_segments = self.get_exon_segments()

        transcript_obj = Transcript.Transcript(contig_acc, exon_segments, contig_strand)
        transcript_obj._multipath = self
        transcript_obj._simplepath = self.get_simple_path()

        transcript_obj.add_multipaths_evidence_assigned(self)
        transcript_obj.set_multipaths_evidence_weights({self: 1.0})

        return transcript_obj

    def __repr__(self):
        # Resolve names on demand for human-readable representation
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

        return (
            self.get_id()
            + ":"
            + str(self._simple_path)
            + " rtypes:"
            + str(self._read_types)
            + " rnames: "
            + read_names_show
        )

    def toShortDescr(self):
        return self.get_id() + ":" + str(self.get_simple_path())

    @staticmethod
    def multiPath_from_mpgn_list(mpgn_list):

        assert len(mpgn_list) > 0

        path_list = list()

        sg = mpgn_list[0].get_splice_graph()

        for mpgn in mpgn_list:
            path_list.append(mpgn.get_simple_path())

        mp_obj = MultiPath(sg, path_list)

        return mp_obj

    def _coerce_read_identifier(self, value):
        """
        Try to coerce an incoming read identifier to an integer ID. If a string name is provided,
        use READ_NAME_STORE to map to an ID when available. Returns None if coercion is not possible.
        """
        try:
            if isinstance(value, int):
                return int(value)
            if isinstance(value, str):
                # numeric string?
                try:
                    return int(value)
                except Exception:
                    pass
                # look up via name store
                try:
                    import LRAA_Globals
                    name_store = getattr(LRAA_Globals, "READ_NAME_STORE", None)
                    if name_store is not None:
                        return int(name_store.get_or_add(value))
                except Exception:
                    return None
                # fall back to deterministic hash-based identifier when no store is available
                digest = hashlib.sha1(value.encode("utf-8", "ignore")).hexdigest()
                return int(digest[:16], 16)
        except Exception:
            return None
        return None

    # Merge helpers (do not change counts)

    def merge_read_ids(self, read_ids):
        """
        Merge a collection of compact read IDs into this MultiPath without adjusting read_count.
        Returns the number of *new* unique IDs added so the caller can keep external
        counters in sync without having to materialize the existing ID set (which can be large).
        """
        if read_ids is None:
            return 0

        added = 0
        try:
            for rid in read_ids:
                try:
                    rid_int = int(rid)
                except Exception:
                    continue

                if rid_int not in self._read_ids:
                    self._read_ids.add(rid_int)
                    added += 1
        except Exception:
            return added

        return added


def __get_dummy_splice_graph():

    sg = Splice_graph()

    Exon.reset_counter()

    #   E1:100-200   E2:300-400      E3:500-600          E4:700-800    E5:900-1000
    #    [-----]     [--------]      [---------]         [--------]    [---------]
    #

    e1 = Exon("contig", 100, 200, "+", 1)
    e1_ID = e1.get_id()
    sg._node_id_to_node[e1_ID] = e1

    e2 = Exon("contig", 300, 400, "+", 1)
    e2_ID = e2.get_id()
    sg._node_id_to_node[e2_ID] = e2

    e3 = Exon("contig", 500, 600, "+", 1)
    e3_ID = e3.get_id()
    sg._node_id_to_node[e3_ID] = e3

    e4 = Exon("contig", 700, 800, "+", 1)
    e4_ID = e4.get_id()
    sg._node_id_to_node[e4_ID] = e4

    e5 = Exon("contig", 900, 1000, "+", 1)
    e5_ID = e5.get_id()
    sg._node_id_to_node[e5_ID] = e5

    # print(str(sg._node_id_to_node))

    return sg


def test_overlapping_n_compatible():

    sg = __get_dummy_splice_graph()

    mp1 = MultiPath(sg, [["E:1", "E:2", "E:3"]])
    mp2 = MultiPath(sg, [["E:2", "E:3", "E:4"]])

    # test compatible paths - no spacers
    assert mp1.is_overlapping_and_compatible(mp2) == True
    assert mp2.is_overlapping_and_compatible(mp1) == True

    # test incompatible paths - no spacers
    mp3 = MultiPath(sg, [["E:1", "E:2", "E:4"]])
    assert mp1.is_overlapping_and_compatible(mp3) == False
    assert mp3.is_overlapping_and_compatible(mp1) == False

    assert mp3.is_overlapping_and_compatible(mp2) == False
    assert mp2.is_overlapping_and_compatible(mp3) == False

    # test compatible paths with spacers
    mp_sp1 = MultiPath(sg, [["E:1", SPACER, "E:4"]])
    mp_sp2 = MultiPath(sg, [["E:1", SPACER, "E:3", "E:4"]])
    assert mp_sp1.is_overlapping_and_compatible(mp_sp2) == False

    mp_sp3 = MultiPath(sg, [["E:1", SPACER, "E:3"]])
    assert mp_sp1.is_overlapping_and_compatible(mp_sp3) == False

    # test incompatible paths with spacers
    mp4 = MultiPath(sg, [["E:2", "E:3", "E:5"]])
    assert mp_sp2.is_overlapping_and_compatible(mp4) == False

    # test multiple spacers
    mp_sp4 = MultiPath(sg, [["E:1", SPACER, "E:2", "E:3", "E:4", "E:5"]])
    mp_sp5 = MultiPath(sg, [["E:1", "E:2", "E:3", "E:4", SPACER, "E:5"]])
    assert mp_sp4.is_overlapping_and_compatible(mp_sp5) == False

    # test incompatible multipel spacers
    mp_sp6 = MultiPath(sg, [["E:1", "E:2", "E:5"]])
    assert mp_sp6.is_overlapping_and_compatible(mp_sp4) == False


def test_merge_paths_to_simple_multi_path():

    sg = __get_dummy_splice_graph()

    paths_list = [["E:1", "E:2"], ["E:2", "E:3"]]
    mp = MultiPath(sg, paths_list)

    assert mp.get_simple_path() == ["E:1", "E:2", "E:3"]

    paths_list = [["E:1", "E:2"], ["E:3", "E:4"]]
    mp = MultiPath(sg, paths_list)
    assert mp.get_simple_path() == ["E:1", "E:2", SPACER, "E:3", "E:4"]
