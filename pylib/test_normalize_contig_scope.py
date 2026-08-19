#!/usr/bin/env python3

"""C3: --contig / --restrict_to_chromosomes must scope normalization, not just the quant
loop.

Before this change, `_normalize_bam_for_splice_graph`'s two call sites always passed
`scope=None` to `Util_funcs.splice_graph_norm_cache_stem` and always handed
`normalize_bam_by_strand.py` the whole input bam, regardless of `--contig` /
`--restrict_to_chromosomes`. Measured on a real fleet run, `--contig chr21` spent 3,324
of 3,474 total seconds normalizing the FULL library while chr21 held 1.56% of its reads
(`BENCHMARKING.chunked_runs.md`, `[HG002-fleet]`). This file tests the fix: contig
restriction is resolved once, early (`_resolve_contig_restriction`), and threaded into
`_normalize_bam_for_splice_graph` as `contig_scope`, which both (a) pre-filters the
normalizer's input to just those contigs and (b) names that restriction in the cache
stem so a restricted run cannot collide with -- or be silently served by -- a whole-bam
cache entry for the same source file.
"""

import importlib.machinery
import importlib.util
import os
import subprocess
import sys
from types import SimpleNamespace

import pysam
import pytest

REPO_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
for _sub in ("pylib", "util"):
    _p = os.path.join(REPO_ROOT, _sub)
    if _p not in sys.path:
        sys.path.insert(0, _p)

LRAA_BIN = os.path.join(REPO_ROOT, "LRAA")


def _load_lraa():
    loader = importlib.machinery.SourceFileLoader(
        "lraa_contig_scope_test_module", LRAA_BIN
    )
    spec = importlib.util.spec_from_loader(loader.name, loader)
    module = importlib.util.module_from_spec(spec)
    loader.exec_module(module)
    return module


lraa = _load_lraa()


def _args(**overrides):
    values = {"contig": None, "region": None, "restrict_to_chromosomes": None}
    values.update(overrides)
    return SimpleNamespace(**values)


# --------------------------------------------------------------- _resolve_contig_restriction


def _genome(tmp_path, contigs=("chr1", "chr2", "chr3")):
    fasta = tmp_path / "genome.fa"
    with open(fasta, "w") as fh:
        for c in contigs:
            fh.write(">{}\n{}\n".format(c, "A" * 200))
    pysam.faidx(str(fasta))
    return str(fasta)


def test_no_restriction_flags_returns_the_whole_genome_contig_listing(tmp_path):
    genome = _genome(tmp_path)
    contigs, restricted, orient = lraa._resolve_contig_restriction(_args(), genome)
    assert contigs == ["chr1", "chr2", "chr3"]
    assert restricted is None
    assert orient is None


def test_contig_flag_yields_a_single_contig_scope(tmp_path):
    genome = _genome(tmp_path)
    contigs, restricted, orient = lraa._resolve_contig_restriction(
        _args(contig="chr2"), genome
    )
    assert contigs == ["chr2"]
    assert restricted == "chr2"
    assert orient is None


def test_contig_flag_strips_the_orientation_suffix_from_the_scope():
    """The scope is WHOLE contigs; +/- selects a strand within one, not a different bam
    region, so it must not leak into what gets passed to samtools view."""
    contigs, restricted, orient = lraa._resolve_contig_restriction(
        _args(contig="chr2+"), "unused.fa"
    )
    assert contigs == ["chr2"]
    assert restricted == "chr2"
    assert orient == "+"


def test_restrict_to_chromosomes_filters_the_genome_contig_listing(tmp_path):
    genome = _genome(tmp_path)
    contigs, restricted, orient = lraa._resolve_contig_restriction(
        _args(restrict_to_chromosomes="chr1,chr3"), genome
    )
    assert contigs == ["chr1", "chr3"]
    assert restricted is None
    assert orient is None


def test_restrict_to_chromosomes_naming_an_absent_contig_is_dropped_silently(tmp_path):
    """Matches the pre-existing behavior this function was extracted from: a name not
    in the genome fasta is filtered out rather than raising."""
    genome = _genome(tmp_path)
    contigs, _, _ = lraa._resolve_contig_restriction(
        _args(restrict_to_chromosomes="chr1,chrNotReal"), genome
    )
    assert contigs == ["chr1"]


def test_contig_flag_takes_precedence_over_restrict_to_chromosomes(tmp_path):
    """The two flags have never combined: --contig alone decides genome_contigs_list
    when both are given, exactly as before this refactor."""
    genome = _genome(tmp_path)
    contigs, restricted, _ = lraa._resolve_contig_restriction(
        _args(contig="chr2", restrict_to_chromosomes="chr1,chr3"), genome
    )
    assert contigs == ["chr2"]
    assert restricted == "chr2"


# --------------------------------------------------------- _normalize_bam_for_splice_graph


CONTIG_A = "chrA"
CONTIG_B = "chrB"
CONTIG_A_LEN = 5000
CONTIG_B_LEN = 5000


def _header():
    return pysam.AlignmentHeader.from_dict(
        {
            "HD": {"VN": "1.6", "SO": "coordinate"},
            "SQ": [
                {"SN": CONTIG_A, "LN": CONTIG_A_LEN},
                {"SN": CONTIG_B, "LN": CONTIG_B_LEN},
            ],
        }
    )


def _read(header, name, contig, start, length=100):
    a = pysam.AlignedSegment(header)
    a.query_name = name
    a.reference_id = header.get_tid(contig)
    a.reference_start = start
    a.mapping_quality = 60
    a.cigarstring = "{}M".format(length)
    a.query_sequence = "A" * length
    a.query_qualities = pysam.qualitystring_to_array("I" * length)
    return a


def _two_contig_bam(path):
    """chrA: 40 overlapping reads at one locus, deep enough that a low
    normalize_max_cov_level really thins it. chrB: 5 reads elsewhere, entirely
    disjoint. Distinguishable by contig alone, which is exactly what --contig /
    --restrict_to_chromosomes scope on."""
    header = _header()
    with pysam.AlignmentFile(str(path), "wb", header=header) as fh:
        for i in range(40):
            fh.write(_read(header, "readA{}".format(i), CONTIG_A, 1000 + i))
        for i in range(5):
            fh.write(_read(header, "readB{}".format(i), CONTIG_B, 2000 + i))
    pysam.index(str(path))
    return path


@pytest.fixture
def two_contig_bam(tmp_path):
    return _two_contig_bam(tmp_path / "two_contig.bam")


@pytest.fixture(autouse=True)
def _norm_config(monkeypatch):
    """Pin the identity/mapq floors these fixtures assume, and restore them after each
    test via monkeypatch so this file cannot leak config state into tests that run
    later in the same pytest session (LRAA_Globals.config is a shared module-level
    dict, not per-test)."""
    monkeypatch.setitem(lraa.LRAA_Globals.config, "min_per_id", 0)
    monkeypatch.setitem(lraa.LRAA_Globals.config, "min_mapping_quality", 0)
    monkeypatch.setitem(
        lraa.LRAA_Globals.config, "min_mapping_quality_for_final_quant", 0
    )
    monkeypatch.setitem(lraa.LRAA_Globals.config, "max_intron_length", 200000)


def _normalize(bam, cache_dir, *, level=10, contig_scope=None, label=None):
    return lraa._normalize_bam_for_splice_graph(
        str(bam), level, str(cache_dir), True, label=label, contig_scope=contig_scope
    )


def _read_xw(bam_path):
    """(name, flag, start, contig) -> XW weight, for every retained record."""
    out = {}
    with pysam.AlignmentFile(str(bam_path), "rb") as fh:
        for read in fh.fetch(until_eof=True):
            key = (read.query_name, read.flag, read.reference_start, read.reference_name)
            out[key] = read.get_tag("XW") if read.has_tag("XW") else None
    return out


def test_cache_stem_is_scopenone_when_unrestricted(tmp_path, two_contig_bam):
    out = _normalize(two_contig_bam, tmp_path / "cache", contig_scope=None)
    assert ".scopenone." in out


def test_cache_stem_names_the_restricted_scope(tmp_path, two_contig_bam):
    out = _normalize(two_contig_bam, tmp_path / "cache", contig_scope=[CONTIG_A])
    assert ".scope{}.".format(CONTIG_A) in out
    assert ".scopenone." not in out


def test_scoped_normalization_never_emits_records_from_another_contig(
    tmp_path, two_contig_bam
):
    """The load-bearing behavior: a --contig run's normalized bam must contain only that
    contig's records. A regression that forgot to pre-filter the input (passing scope
    only to the cache stem, as the pre-fix code effectively did with scope=None always)
    would emit chrB records here too, since chrB is well below the target depth and
    every one of its reads would be kept untouched."""
    out = _normalize(two_contig_bam, tmp_path / "cache", contig_scope=[CONTIG_A])
    with pysam.AlignmentFile(out, "rb") as fh:
        contigs_seen = {r.reference_name for r in fh.fetch(until_eof=True)}
    assert contigs_seen == {CONTIG_A}


def test_scoped_normalization_is_byte_identical_to_whole_bam_normalization_filtered_after(
    tmp_path, two_contig_bam
):
    """The load-bearing correctness claim (C3 acceptance #2): restricting the
    normalizer's INPUT to one contig must not change that contig's retained-read set or
    any of its XW weights, versus normalizing the whole library and keeping only that
    contig's output records. True because normalize_bam_by_strand.py's depth
    accumulation is keyed per contig (window_bases) and this driver always fixes
    window_origin, which anchors window boundaries at absolute coordinates -- so which
    OTHER contigs are present cannot move chrA's own windows or measured depth.
    """
    whole_out = _normalize(
        two_contig_bam, tmp_path / "whole_cache", contig_scope=None, label="whole"
    )
    scoped_out = _normalize(
        two_contig_bam,
        tmp_path / "scoped_cache",
        contig_scope=[CONTIG_A],
        label="whole",
    )

    whole_xw = _read_xw(whole_out)
    whole_chrA = {k: v for k, v in whole_xw.items() if k[3] == CONTIG_A}
    scoped_xw = _read_xw(scoped_out)

    assert len(whole_chrA) > 0
    # confirm real thinning happened, not a trivial pass-through
    assert any(v is not None and v > 1.0 for v in whole_chrA.values())
    assert whole_chrA == scoped_xw


def test_a_contig_absent_from_the_bam_header_yields_an_empty_input_not_the_whole_bam(
    tmp_path, two_contig_bam
):
    """A requested contig that is not in the bam's header at all (as opposed to being in
    the header with zero reads) must not fall back to "no filter" -- that would silently
    normalize the whole, unscoped source under a name that claims it is restricted."""
    out = _normalize(two_contig_bam, tmp_path / "cache", contig_scope=["chrNotInBam"])
    with pysam.AlignmentFile(out, "rb") as fh:
        records = list(fh.fetch(until_eof=True))
    assert records == []


def test_a_cache_hit_does_no_prefiltering_work(tmp_path, two_contig_bam):
    """The prefiltered bam is scoped by the same cache stem as the final output and must
    be gated behind the same checkpoint: a rerun that already has the normalized bam
    cached must not touch the (possibly cleaned-up) prefiltered intermediate at all."""
    cache_dir = tmp_path / "cache"
    out = _normalize(two_contig_bam, cache_dir, contig_scope=[CONTIG_A])

    scoped_intermediate = [
        p for p in cache_dir.iterdir() if p.name.endswith(".contigscope.bam")
    ]
    assert scoped_intermediate
    for p in scoped_intermediate:
        p.unlink()
    for p in cache_dir.glob("*.contigscope.bam.ok"):
        p.unlink()

    # rerun: normalized_path's own checkpoint still exists, so this must return
    # immediately without trying to recreate the intermediate it just deleted
    out2 = _normalize(two_contig_bam, cache_dir, contig_scope=[CONTIG_A])
    assert out2 == out
    assert not any(
        p.name.endswith(".contigscope.bam") for p in cache_dir.iterdir()
    )


def test_restricted_and_unrestricted_scopes_use_different_cache_entries(
    tmp_path, two_contig_bam
):
    """A2's collision guard, exercised end to end: normalizing chrA alone must not reuse
    -- or be reused by -- a whole-library normalization of the same source bam."""
    whole_out = _normalize(two_contig_bam, tmp_path / "cache", contig_scope=None)
    scoped_out = _normalize(
        two_contig_bam, tmp_path / "cache", contig_scope=[CONTIG_A]
    )
    assert whole_out != scoped_out


# --------------------------------------------------------------------- end to end (CLI)


def _fai(tmp_path):
    gtf = tmp_path / "a.gtf"
    gtf.write_text(
        '{}\ttest\texon\t1001\t1150\t.\t+\t.\tgene_id "g1"; transcript_id "t1";\n'.format(
            CONTIG_A
        )
    )
    genome = tmp_path / "g.fa"
    with open(genome, "w") as fh:
        fh.write(">{}\n{}\n".format(CONTIG_A, "A" * CONTIG_A_LEN))
        fh.write(">{}\n{}\n".format(CONTIG_B, "A" * CONTIG_B_LEN))
    return gtf, genome


def test_contig_cli_flag_scopes_normalization_cache_end_to_end(tmp_path):
    """Full CLI invocation: `--contig chrA` must produce a norm-cache file named
    `.scopechrA.`, not `.scopenone.` -- the observable, on-disk proof that C3's fix is
    wired all the way from argv to the normalizer, not just reachable in a unit test."""
    bam = _two_contig_bam(tmp_path / "reads.bam")
    gtf, genome = _fai(tmp_path)
    cmd = [
        sys.executable,
        LRAA_BIN,
        "--bam", str(bam),
        "--genome", str(genome),
        "--gtf", str(gtf),
        "--quant_only",
        "--contig", CONTIG_A,
        "--normalize_max_cov_level", "10",
        "--output_prefix", str(tmp_path / "out"),
    ]
    result = subprocess.run(
        cmd, capture_output=True, text=True, cwd=str(tmp_path), timeout=300
    )
    assert result.returncode == 0, (result.stdout + result.stderr)[-4000:]

    cache_dirs = [
        p
        for p in tmp_path.iterdir()
        if p.is_dir() and p.name.endswith(".norm_cache")
    ]
    assert cache_dirs, "no normalization cache directory was created"
    norm_bams = list(cache_dirs[0].glob("*.bam"))
    assert norm_bams, "no normalized bam was written"
    assert any(".scope{}.".format(CONTIG_A) in p.name for p in norm_bams)
    assert not any(".scopenone." in p.name for p in norm_bams)

    # the normalized bam itself must carry only chrA -- chrB was never read by
    # normalization at all, matching what the quant loop will restrict to anyway
    scoped_bam = next(p for p in norm_bams if p.name.endswith(".bam"))
    with pysam.AlignmentFile(str(scoped_bam), "rb") as fh:
        contigs_seen = {r.reference_name for r in fh.fetch(until_eof=True)}
    assert contigs_seen <= {CONTIG_A}
