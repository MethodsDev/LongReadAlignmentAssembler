#!/usr/bin/env python3

"""Pipeliner.__init__ must not race on creating its own checkpoint_dir.

Found via a genuinely flaky (~1/3) failure in test_chunk_forwards_single_cell.py
once a strandless chunk's two orientations started running as sibling threads:
each spawns a normalize_bam_by_strand.py subprocess, and that script constructs
its own Pipeliner("__chckpts") rooted at the SAME chunk directory. The old
`if not os.path.exists(d): os.makedirs(d)` is a classic check-then-act race --
two callers can both observe "not exists" before either finishes creating it,
and the loser gets FileExistsError. Nothing about this is specific to two
threads in one Python process either: two independent subprocesses hit the
same OS-level race the same way, which is the actual shape that reproduced it.
"""

import os
import shutil
import sys
import threading

REPO = os.path.abspath(os.path.join(os.path.dirname(os.path.realpath(__file__)), ".."))
sys.path.insert(0, os.path.join(REPO, "pylib"))

import Pipeliner  # noqa: E402  (path insert must precede the import)


def test_concurrent_construction_never_raises(tmp_path):
    """A barrier forces every thread to call makedirs at nearly the same
    instant, which is what turns a rare race into a reliable repro -- run
    several rounds against a fresh, nonexistent directory each time, since a
    single lucky pass proves nothing."""

    for round_num in range(20):
        checkpoint_dir = str(tmp_path / "chckpts_{}".format(round_num))
        assert not os.path.exists(checkpoint_dir)

        n = 8
        barrier = threading.Barrier(n)
        errors = []
        lock = threading.Lock()

        def make_one():
            barrier.wait()
            try:
                Pipeliner.Pipeliner(checkpoint_dir)
            except Exception as err:  # noqa: BLE001 - captured for the assertion below
                with lock:
                    errors.append(err)

        threads = [threading.Thread(target=make_one) for _ in range(n)]
        for t in threads:
            t.start()
        for t in threads:
            t.join()

        assert not errors, (
            "round {}: {} of {} concurrent Pipeliner() constructions raised: {}".format(
                round_num, len(errors), n, errors
            )
        )
        assert os.path.isdir(checkpoint_dir)


def test_reusing_an_existing_checkpoint_dir_still_works(tmp_path):
    """The ordinary sequential case -- a resumed run pointing at a directory
    that already exists -- must keep working exactly as before."""

    checkpoint_dir = str(tmp_path / "chckpts")
    os.makedirs(checkpoint_dir)
    sentinel = os.path.join(checkpoint_dir, "already.here")
    open(sentinel, "w").close()

    Pipeliner.Pipeliner(checkpoint_dir)

    assert os.path.exists(sentinel), "an existing directory's contents must survive"
