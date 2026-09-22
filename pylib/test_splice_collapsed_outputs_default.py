#!/usr/bin/env python3

"""The splice-collapsed gtf and site beds default ON, and chunk workers must opt OUT.

A chunk worker is a fresh `LRAA --no_chunk` process holding ONE chunk, so any output
that defaults on reaches it by default too. Its models are a fragment of the run, the
root derives the real artifacts from the merged gtf, and a per-unit collapse would
therefore be both wasted and misleading -- a `<unit>.TSS.bed` describing one chunk sits
in the work directory looking exactly like the run's answer.

`lraa_cmd` builds the worker argv from an explicit allowlist rather than forwarding the
root's, so "the root did not ask for it" is NOT what keeps it off the worker: the
worker resolves its own default. This is the same trap `--no_chunk` documents at that
call site, where a default flipping to True turned an omitted flag into a run-killing
re-entry. Pinned here because the failure is silent: the run still succeeds and still
publishes correct top-level files.
"""

import argparse
import os
import sys

import pytest

REPO_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if os.path.join(REPO_ROOT, "pylib") not in sys.path:
    sys.path.insert(0, os.path.join(REPO_ROOT, "pylib"))

import ChunkedRun

FLAG = "--no_splice_collapsed_outputs"


def _args(**overrides):
    args = argparse.Namespace(
        discovery=True,
        HiFi=False,
        no_rdna_mask=False,
        rdna_mask_fasta=None,
        cell_list=None,
        min_mapping_quality=0,
        min_mapping_quality_for_final_quant=0,
        stream_reads=False,
        stream_reads_rescue_unassigned=False,
        rescue_unassigned_reads_via_transcriptome_alignment=False,
    )
    for key, value in overrides.items():
        setattr(args, key, value)
    return args


def _cmd(**overrides):
    return ChunkedRun.lraa_cmd(
        _args(**overrides),
        bam_for_quant="quant.bam",
        bam_for_sg="sg.bam",
        genome="genome.fa",
        gtf=None,
        out_prefix="unit",
        num_total_reads=1000,
        cpu_budget=1,
    )


def test_chunk_workers_are_told_not_to_write_them():
    assert FLAG in _cmd(), "a chunk worker would collapse its own fragment"


def test_workers_opt_out_in_quant_only_too():
    """Quant-only emits no gtf to collapse, so the flag is inert there -- but the
    argv must not depend on that, or the guarantee moves into another file."""
    assert FLAG in _cmd(discovery=False)


def test_the_flag_is_a_real_opt_out_and_not_a_typo():
    """A misspelled flag is accepted by nothing and would abort every chunk worker.

    Checked against the CLI itself rather than a copy of its option list, so renaming
    the option fails here instead of at the first chunked run.  --show_full_usage_info
    rather than --help: the latter prints a curated short list that omits most options.
    """
    import subprocess

    help_text = subprocess.run(
        [sys.executable, os.path.join(REPO_ROOT, "LRAA"), "--show_full_usage_info"],
        capture_output=True,
        text=True,
    ).stdout

    assert FLAG in help_text
    # An opt-out, not an opt-in: no --include_... remains, so a caller cannot be
    # carrying a flag that silently stopped meaning anything.
    assert "--include_splice_collapsed_outputs" not in help_text
