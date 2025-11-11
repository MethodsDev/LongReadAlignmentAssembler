#!/usr/bin/env python
# encoding: utf-8

import sys, os, re
from collections import defaultdict
import logging
import MultiPath
import LRAA_Globals
from LRAA_Globals import SPACER

logger = logging.getLogger(__name__)


class MultiPathCountPair:

    def __init__(self, multipath):
        self._multipath = multipath
        self.reset_count()

    def increment(self, increment=1):
        # Keep the pair count and the underlying multipath's internal read_count in sync
        try:
            inc = int(increment)
        except Exception:
            inc = 1
        if inc <= 0:
            return
        self._count += inc
        # Also bump the multipath's internal count so downstream users of mp.get_read_count()
        # (e.g., EM, equal assignment) see the aggregated count
        try:
            self._multipath.include_read_count(inc)
        except Exception:
            pass

    def get_multipath_and_count(self):
        return (self._multipath, self._count)

    def include_read_type(self, read_type):
        self._multipath.include_read_type(read_type)

    def include_read_name(self, read_name):
        # retained for backward compatibility; prefer include_read_count in counter
        if type(read_name) in [list, set]:
            self._multipath.include_read_count(len(read_name))
        else:
            self._multipath.include_read_count(1)

    def reset_count(self):
        self._count = self._multipath.get_read_count()

    def __repr__(self):
        ret_text = "{}\t{}".format(str(self._multipath), self._count)
        return ret_text

    def prune_reftranscript_as_evidence(self):
        self._multipath.prune_reftranscript_as_evidence()
        self.reset_count()

    def get_count(self):
        return self._count


class MultiPathCounter:

    def __init__(self):

        self._multipath_counter = dict()

        return

    def add(self, multipath_obj):

        # count must be positive
        assert multipath_obj.get_read_count() > 0

        assert type(multipath_obj) == MultiPath.MultiPath
        multipath_key = multipath_obj.get_simple_path_tuple()

        if multipath_key in self._multipath_counter:
            orig_mp_count_pair = self._multipath_counter[multipath_key]
            orig_mp, _ = orig_mp_count_pair.get_multipath_and_count()

            # Always union read types
            orig_mp_count_pair.include_read_type(multipath_obj.get_read_types())

            # Prefer uniqueness-aware increment when we have explicit IDs
            try:
                incoming_ids = multipath_obj.get_read_ids()
            except Exception:
                incoming_ids = set()

            if incoming_ids:
                # merge provenance IDs first and learn how many were genuinely new
                try:
                    if hasattr(orig_mp, "merge_read_ids"):
                        newly_added = orig_mp.merge_read_ids(incoming_ids)
                    else:
                        newly_added = 0
                except Exception:
                    newly_added = 0

                # increment count only by truly new unique IDs
                if newly_added and newly_added > 0:
                    orig_mp_count_pair.increment(newly_added)

                # DEBUG guardrail: ensure count equals number of unique IDs
                if LRAA_Globals.DEBUG:
                    unique_ct = orig_mp.get_read_id_count()
                    if unique_ct != orig_mp_count_pair.get_count():
                        raise RuntimeError(
                            f"MultiPathCounter invariant violated: count={orig_mp_count_pair.get_count()} != unique_ids={unique_ct} for path {multipath_key}"
                        )
            else:
                # Always require explicit read IDs for aggregation
                raise RuntimeError(
                    "Incoming multipath lacks read IDs; aggregation requires explicit IDs to ensure uniqueness."
                )

            return orig_mp_count_pair

        else:

            if LRAA_Globals.DEBUG:
                # verify this multipath obj is not already stored.
                logger.debug(
                    "-verifying simple paths are unique in the multpath counter and lack spacers"
                )
                for mp_count_pair_obj in self._multipath_counter.values():
                    mp, count = mp_count_pair_obj.get_multipath_and_count()
                    if (
                        mp == multipath_obj
                        or mp.get_simple_path() == multipath_obj.get_simple_path()
                    ):

                        raise RuntimeError(
                            "Encountered multipath as value already.  Incoming: {}, already stored {}".format(
                                str(multipath_obj), str(mp)
                            )
                        )
                    if SPACER in mp.get_simple_path():
                        raise RuntimeError(
                            "Encountered unallowed SPACER in multipath: {}".format(
                                mp.get_simple_path()
                            )
                        )
            self._multipath_counter[multipath_key] = MultiPathCountPair(multipath_obj)

            return self._multipath_counter[multipath_key]

    def get_all_MultiPathCountPairs(self):
        # self.validate_MultiPathCounter()
        return self._multipath_counter.values()

    def __repr__(self):

        return "\n".join([str(x) for x in self._multipath_counter.values()])

    def validate_MultiPathCounter(self):

        # verify this multipath obj is not already stored.
        logger.debug("-verifying simple paths are unique in the multipath counter")

        mp_keys = list()
        mp_count_pairs = list()
        for x, y in self._multipath_counter.items():
            mp_keys.append(x)
            mp_count_pairs.append(y)

        for i in range(0, len(mp_count_pairs) - 1):

            mp_i, count_i = mp_count_pairs[i].get_multipath_and_count()

            for j in range(i + 1, len(mp_count_pairs)):

                mp_j, count_j = mp_count_pairs[j].get_multipath_and_count()

                # logger.debug("Comparing mps: {} vs. {}".format(mp_i, mp_j))

                if mp_i == mp_j:

                    raise RuntimeError(
                        "Encountered multipath stored under two keys:\nkey1: {}\nkey2:{}\nmp:{}".format(
                            mp_keys[i], mp_keys[j], mp_i
                        )
                    )

                if mp_i.get_simple_path() == mp_j.get_simple_path():
                    raise RuntimeError(
                        "Encountered two multipaths stored separately but have the same simple path:\n"
                        + "simple path: {}\n".format(str(mp_i.get_simple_path))
                        + "mp_i: {}\n".format(mp_i)
                        + "mp_j: {}\n".format(mp_j)
                    )

        logger.debug("MultPathCounter validates.")

        return

    def prune_ref_transcripts_as_evidence(self):
        multipaths_to_purge = set()
        for mp_key, mp_count_pair in self._multipath_counter.items():
            mp_count_pair.prune_reftranscript_as_evidence()
            if mp_count_pair.get_count() == 0:
                multipaths_to_purge.add(mp_key)

        if len(multipaths_to_purge) > 0:
            for mp_key in multipaths_to_purge:
                del self._multipath_counter[mp_key]

        return
