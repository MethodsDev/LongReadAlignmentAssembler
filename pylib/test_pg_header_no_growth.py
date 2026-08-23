#!/usr/bin/env python3

"""Tests that the BAM-writing steps do not grow the @PG header chain.

samtools appends one @PG record PER EXISTING CHAIN TIP on every write, and
`samtools merge` first concatenates all of its inputs' chains. On a BAM whose
header already carries many parallel chains this is multiplicative, not
additive. Measured on XP132160.ucsc.bam, which arrives from upstream with
34,976 @PG records across 5,824 parallel chains (one per minimap2 alignment
shard, never collapsed): the cluster-guided path amplified that 78x, to
2,727,296 records = 1,188,154,439 bytes of uncompressed SAM header text. That
matters beyond tidiness because resolving a region NAME to a tid forces a full
header parse, so each per-chromosome region query against the merged BAM cost
~5 minutes of pure parsing, invariant of how many alignments it returned.

`--no-PG` suppresses only the records a given invocation would add; inherited
chains still concatenate, so it caps growth rather than eliminating it. These
tests pin the cap, which is the property a refactor is most likely to lose
silently: dropping the flag reintroduces per-tip growth with no error, no
warning, and no output difference a normal test would notice.

Header @PG counts are read straight out of the BGZF stream rather than via
`samtools view -H`, because samtools appends its own per-tip records to the
header it emits -- a 34,976-record header reads as 40,800 through samtools.
Measuring this with samtools would report the measurement artifact.
"""

import gzip
import re
import struct
import subprocess
import sys
from pathlib import Path

import pysam
import pytest


def _raw_header_text(bam_path):
    """Header text straight from the BGZF stream, without invoking samtools."""
    with gzip.open(bam_path, "rb") as fh:
        magic = fh.read(4)
        assert magic == b"BAM\x01", f"not a BAM: {magic!r}"
        (l_text,) = struct.unpack("<i", fh.read(4))
        return fh.read(l_text).decode("utf-8", "replace")


def _pg_records(bam_path):
    return [
        line
        for line in _raw_header_text(bam_path).split("\n")
        if line.startswith("@PG")
    ]


def _tips(bam_path):
    """@PG IDs no other record points at via PP -- the append fan-out factor."""
    fields = [
        dict(kv.split(":", 1) for kv in rec.split("\t")[1:] if ":" in kv)
        for rec in _pg_records(bam_path)
    ]
    referenced = {f["PP"] for f in fields if "PP" in f}
    return [f["ID"] for f in fields if f.get("ID") not in referenced]


def _write_bam(path, n_chains, records=2):
    """A BAM with n_chains PARALLEL @PG chains of two records each.

    Parallel chains are the point: a single chain would make the append
    behaviour look additive and the test would not detect the regression.
    """
    header = ["@HD\tVN:1.6\tSO:coordinate", "@SQ\tSN:chrT\tLN:1000"]
    for c in range(n_chains):
        header.append(f"@PG\tID:root{c}\tPN:toolA\tVN:1\tCL:toolA shard{c}")
        header.append(f"@PG\tID:leaf{c}\tPN:toolB\tPP:root{c}\tVN:1\tCL:toolB shard{c}")
    body = [
        f"r{i}\t0\tchrT\t{100 + i * 10}\t60\t10M\t*\t0\t0\tACGTACGTAC\tIIIIIIIIII"
        for i in range(records)
    ]
    sam = path.with_suffix(".sam")
    sam.write_text("\n".join(header + body) + "\n")
    # --no-PG on the conversion itself, so the fixture's count is exactly what
    # was asked for rather than that plus one record per chain.
    subprocess.run(
        ["samtools", "view", "-b", "--no-PG", "-o", str(path), str(sam)],
        check=True,
    )
    subprocess.run(["samtools", "index", str(path)], check=True)
    return path


def _record_count(bam_path):
    """Alignment count via a sequential pass, so no index is required."""
    with pysam.AlignmentFile(bam_path, "rb", check_sq=False) as bam:
        return sum(1 for _ in bam.fetch(until_eof=True))


@pytest.fixture
def two_bams(tmp_path):
    a = _write_bam(tmp_path / "a.bam", n_chains=4)
    b = _write_bam(tmp_path / "b.bam", n_chains=4)
    return a, b


def test_fixture_has_parallel_chains(two_bams):
    """Guard the fixture itself: 4 parallel chains, 8 records, 4 tips."""
    a, _ = two_bams
    assert len(_pg_records(a)) == 8
    assert len(_tips(a)) == 4


def test_merge_without_no_pg_grows_per_tip(two_bams):
    """Characterises the defect, so the fix below is not vacuous.

    Without --no-PG a merge inherits both chains AND adds one record per tip of
    each input: 2 x 8 inherited + (4 + 4) added = 24.
    """
    a, b = two_bams
    out = a.parent / "merged_plain.bam"
    subprocess.run(
        ["samtools", "merge", "-f", str(out), str(a), str(b)], check=True
    )
    inherited = len(_pg_records(a)) + len(_pg_records(b))
    added = len(_pg_records(out)) - inherited
    assert added == len(_tips(a)) + len(_tips(b)), (
        "expected one appended @PG per input tip; samtools append behaviour "
        f"may have changed (inherited={inherited}, total={len(_pg_records(out))})"
    )


def test_merge_with_no_pg_adds_nothing(two_bams):
    """The property the pipeline depends on: --no-PG caps growth at inheritance."""
    a, b = two_bams
    out = a.parent / "merged_nopg.bam"
    subprocess.run(
        ["samtools", "merge", "--no-PG", "-f", str(out), str(a), str(b)],
        check=True,
    )
    assert len(_pg_records(out)) == len(_pg_records(a)) + len(_pg_records(b))


def test_pysam_view_with_no_pg_adds_nothing(two_bams):
    """partition_data_by_chromosome.py extracts per chromosome via pysam.view.

    pysam.view is a samtools CLI wrapper, so it appends per tip -- unlike
    pysam.AlignmentFile, which appends nothing and is why the strand splitter
    never contributed to the chain. Asserts the flag reaches samtools rather
    than being silently accepted and dropped.
    """
    a, _ = two_bams
    out = a.parent / "chrT_nopg.bam"
    pysam.view(
        "--no-PG", "-h", "-b", "-o", str(out), str(a), "chrT", catch_stdout=False
    )
    assert len(_pg_records(out)) == len(_pg_records(a))


def test_pysam_view_without_no_pg_grows(two_bams):
    """Positive control for the test above: without the flag it does grow."""
    a, _ = two_bams
    out = a.parent / "chrT_plain.bam"
    pysam.view("-h", "-b", "-o", str(out), str(a), "chrT", catch_stdout=False)
    assert len(_pg_records(out)) == len(_pg_records(a)) + len(_tips(a))


def test_shipped_call_sites_pass_no_pg():
    """The flag is present at every BAM-writing site that appends per tip.

    A refactor that drops it reintroduces multiplicative header growth with no
    functional symptom, so the call sites are pinned by text. pysam.AlignmentFile
    writers are deliberately not required to carry it -- they append nothing.
    """
    repo = Path(__file__).resolve().parents[1]
    expectations = [
        ("util/normalize_bam_by_strand.py", "samtools merge --no-PG"),
        ("WDL/LRAA_quant_by_cluster.wdl", "samtools merge --no-PG"),
        ("util/partition_data_by_chromosome.py", '"--no-PG"'),
    ]
    for relative, needle in expectations:
        text = (repo / relative).read_text()
        assert needle in text, f"{relative} no longer passes --no-PG ({needle!r})"


def test_reheader_collapse_requires_reindexing(two_bams):
    """Why a collapse must precede indexing, pinned as executable documentation.

    Rewriting the header shifts every BGZF virtual offset, so an index built
    before a `samtools reheader -P` is invalid against the result. It fails
    loudly -- non-zero exit, "Invalid BGZF header at offset N" -- which is the
    good outcome, but it does fail, so anyone adding a collapse must index
    after it rather than reusing the pre-collapse .bai. Pinned because the
    ordering constraint is invisible when reviewing either step alone.
    """
    a, _ = two_bams
    sorted_bam = a.parent / "sorted.bam"
    subprocess.run(
        ["samtools", "sort", "--no-PG", "-o", str(sorted_bam), str(a)], check=True
    )
    subprocess.run(["samtools", "index", str(sorted_bam)], check=True)
    expected = subprocess.run(
        ["samtools", "view", "-c", str(sorted_bam), "chrT"],
        capture_output=True, text=True, check=True,
    ).stdout.strip()
    assert int(expected) > 0, "fixture must have records in the queried region"

    minimal = a.parent / "min.sam"
    minimal.write_text("@HD\tVN:1.6\tSO:coordinate\n@SQ\tSN:chrT\tLN:1000\n")
    collapsed = a.parent / "collapsed.bam"
    with open(collapsed, "wb") as fh:
        subprocess.run(
            ["samtools", "reheader", "-P", str(minimal), str(sorted_bam)],
            stdout=fh, check=True,
        )
    assert _pg_records(collapsed) == [], "reheader -P should leave no @PG records"

    # Stale index carried over from before the reheader: must not be trusted.
    stale_bai = collapsed.with_suffix(".bam.bai")
    stale_bai.write_bytes((a.parent / "sorted.bam.bai").read_bytes())
    stale = subprocess.run(
        ["samtools", "view", "-c", str(collapsed), "chrT"],
        capture_output=True, text=True,
    )
    assert stale.returncode != 0, (
        "a pre-reheader .bai must not silently satisfy a region query against "
        f"the reheadered BAM (got rc=0, stdout={stale.stdout.strip()!r})"
    )

    # Rebuilt index: correct answer restored.
    stale_bai.unlink()
    subprocess.run(["samtools", "index", str(collapsed)], check=True)
    fresh = subprocess.run(
        ["samtools", "view", "-c", str(collapsed), "chrT"],
        capture_output=True, text=True, check=True,
    ).stdout.strip()
    assert fresh == expected


def _collapse_script():
    return (
        Path(__file__).resolve().parents[1] / "util" / "misc" / "collapse_bam_pg_header.py"
    )


def test_collapse_removes_inherited_chain(two_bams):
    """The collapse removes accumulated history that --no-PG can only cap.

    --no-PG suppresses the records a write would append but cannot touch chains
    inherited from its inputs; on the measured corpus that left 1,818,752
    records. This is the step that actually empties the chain, and it must
    preserve @HD/@SQ and every alignment while doing so.
    """
    a, b = two_bams
    merged = a.parent / "merged.bam"
    subprocess.run(
        ["samtools", "merge", "--no-PG", "-f", str(merged), str(a), str(b)], check=True
    )
    inherited = len(_pg_records(merged))
    assert inherited > 0, "fixture must inherit a chain for the collapse to remove"
    records_before = _record_count(merged)

    subprocess.run(
        [sys.executable, str(_collapse_script()), "--input_bam", str(merged)],
        check=True, capture_output=True,
    )

    assert _pg_records(merged) == [], "collapse must leave no @PG records"
    assert _record_count(merged) == records_before, "collapse must not drop alignments"
    header = pysam.AlignmentFile(merged, check_sq=False).header.to_dict()
    assert header["HD"]["SO"] == "coordinate", "@HD must survive the collapse"
    assert len(header["SQ"]) > 0, "@SQ must survive the collapse"


def test_collapse_refuses_when_pg_tags_reference_the_chain(two_bams):
    """Refusal, not silent breakage, when an alignment references a @PG record.

    Dropping @PG definitions is only sound if nothing references them. The
    corpus was verified tag-free across 63,154,704 records, but an externally
    supplied BAM may not be, so the check is enforced rather than assumed.
    """
    a, _ = two_bams
    sam = a.parent / "tagged.sam"
    text = subprocess.run(
        ["samtools", "view", "-h", str(a)], capture_output=True, text=True, check=True
    ).stdout
    lines = []
    for line in text.splitlines():
        lines.append(line if line.startswith("@") else line + "\tPG:Z:someprog")
    sam.write_text("\n".join(lines) + "\n")
    tagged = a.parent / "tagged.bam"
    subprocess.run(
        ["samtools", "view", "-b", "--no-PG", "-o", str(tagged), str(sam)], check=True
    )

    before = _pg_records(tagged)
    done = subprocess.run(
        [sys.executable, str(_collapse_script()), "--input_bam", str(tagged)],
        capture_output=True, text=True,
    )
    assert done.returncode == 2, f"expected refusal (2), got {done.returncode}"
    assert "PG:Z:" in done.stderr
    assert _pg_records(tagged) == before, "a refused collapse must not modify the BAM"

    forced = subprocess.run(
        [sys.executable, str(_collapse_script()), "--input_bam", str(tagged), "--force"],
        capture_output=True, text=True,
    )
    assert forced.returncode == 0, "--force must override the refusal"
    assert _pg_records(tagged) == []


def test_wdl_merge_collapses_before_indexing_and_gates_on_pg_tags():
    """The WDL merge task must collapse, gate on PG:Z:, and index last.

    merge_bams is reachable with `pre_normalized_cluster_bams`, which bypasses
    the normalizer that would otherwise have collapsed its inputs, so this task
    has to handle a dirty chain itself. It does so with samtools alone: the task
    runs in a separately pinned image while the .wdl is read from the checkout,
    so calling a newly added helper from here would break on version skew.
    """
    repo = Path(__file__).resolve().parents[1]
    wdl = (repo / "WDL" / "LRAA_quant_by_cluster.wdl").read_text()
    task = wdl[wdl.index("task merge_bams"):]
    command = task[task.index("command <<<"): task.index(">>>", task.index("command <<<"))]

    # Negative assertions must look at executable code, not at the comments
    # explaining why an approach was rejected, nor at the diagnostic strings
    # that name the helper as a remedy for the operator. Both mention these
    # strings deliberately, so strip comments and double-quoted literals.
    code = re.sub(r"#[^\n]*", "", command)
    bare = re.sub(r'"[^"]*"', '""', code)

    assert "-d PG" in code, "must preflight for PG:Z:-tagged records"
    assert "exit 2" in code, "must fail rather than emit an uncollapsed BAM"
    assert "reheader -P -c" in code, "must collapse with reheader, -P to add nothing"
    assert "samtools view -H" not in bare, (
        "must not dump the header with samtools view -H: it appends its own "
        "records and has been seen to hang on headers this size"
    )
    assert "collapse_bam_pg_header.py" not in bare, (
        "must not INVOKE the helper script here: this task runs in a separately "
        "pinned image while the .wdl is read from the checkout. Naming it inside "
        "a diagnostic message is fine; calling it is not."
    )

    # The collapse rewrites the largest artifact the workflow produces, and the
    # tag preflight is a full decode pass over it. Both must be skipped when the
    # header is already empty, which is the normal route. Keyed on the measured
    # header rather than the input route: provenance would wrongly skip a real
    # chain when a newer .wdl runs against an image lacking the upstream
    # collapse.
    assert "samtools view -H" not in bare, "header count must not use samtools"
    assert 'count(b"\\n@PG")' in code, (
        "must count @PG exactly from the raw BGZF header, streaming; samtools "
        "appends one record per chain tip (measured: 1 for an empty chain, 85 "
        "for 84) and materializing the header needs 2.2-5.8x its size"
    )
    assert "gzip.open" in code and "CHUNK" in code, (
        "the count must read the header in bounded chunks, not all at once"
    )
    count_at = code.index('count(b"\\n@PG")')
    for guarded in ("-d PG", "reheader -P -c"):
        assert count_at < code.index(guarded), (
            f"{guarded!r} must be gated behind the @PG count, not run "
            "unconditionally on every merge"
        )

    gate = code.index("-d PG")
    assert gate < code.index("reheader -P -c"), "gate must precede the collapse"
    assert code.index("reheader -P -c") < code.index("samtools index"), (
        "collapse must precede indexing: rewriting the header invalidates a .bai"
    )


def _load_collapse_module():
    import importlib.util
    path = Path(__file__).resolve().parents[1] / "util" / "misc" / "collapse_bam_pg_header.py"
    spec = importlib.util.spec_from_file_location("collapse_bam_pg_header", path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def _bam_with_pg_records(tmp_path, n, name="h.bam"):
    sam = tmp_path / (name + ".sam")
    with open(sam, "w") as fh:
        fh.write("@HD\tVN:1.6\tSO:coordinate\n@SQ\tSN:chrT\tLN:1000\n")
        for i in range(n):
            fh.write("@PG\tID:p%d\tPN:minimap2\tVN:1.21\tCL:mm2 %s\n" % (i, "y" * 200))
        fh.write("r1\t0\tchrT\t100\t60\t10M\t*\t0\t0\tACGTACGTAC\tIIIIIIIIII\n")
    bam = tmp_path / name
    subprocess.run(
        ["samtools", "view", "-b", "--no-PG", "-o", str(bam), str(sam)], check=True
    )
    return bam


@pytest.mark.parametrize("n", [0, 1, 2, 24, 84, 5000])
def test_streaming_pg_count_matches_htslib(tmp_path, n):
    """The streaming counter must agree with htslib at every chain size.

    It reads the header in fixed chunks and counts b"\n@PG" across a 3-byte
    overlap, so an off-by-one at a chunk boundary or at the first line is the
    plausible bug. Cross-checked against pysam, which parses the header
    properly, and against the raw text.
    """
    module = _load_collapse_module()
    bam = _bam_with_pg_records(tmp_path, n)
    assert module.count_pg(bam) == n
    truth = len(pysam.AlignmentFile(bam, check_sq=False).header.to_dict().get("PG", []))
    assert module.count_pg(bam) == truth
    assert len(_pg_records(bam)) == n


def test_streaming_count_allocates_independently_of_header_size(tmp_path):
    """Bounded allocation is the point: the check must not scale with the header.

    A materializing count needs ~2.2x the header text (measured: 191 MB peak RSS
    on an 85 MB header) and pysam's to_dict ~5.8x (469 MB). On the 1,188,154,439
    byte header this code exists to remove, those are multi-GB -- the check would
    fail on exactly the input it is meant to detect.

    Measured with tracemalloc rather than peak RSS on purpose: ru_maxrss is
    inherited across fork(), so a child spawned from a large parent (pytest is
    ~67 MB) reports the parent's high-water mark and the signal disappears
    entirely -- observed returning an identical 68,212 kB for headers 2,511x
    apart. tracemalloc counts this interpreter's own allocations and cannot be
    contaminated that way.
    """
    module = _load_collapse_module()
    small = _bam_with_pg_records(tmp_path, 20, "small.bam")
    large = _bam_with_pg_records(tmp_path, 50000, "large.bam")

    def header_bytes(path):
        with gzip.open(path, "rb") as fh:
            fh.read(4)
            return struct.unpack("<i", fh.read(4))[0]

    growth = header_bytes(large) / header_bytes(small)
    assert growth > 100, "fixture must span a wide header-size range"

    probe = (
        "import sys, tracemalloc, importlib.util;"
        "spec = importlib.util.spec_from_file_location('c', sys.argv[1]);"
        "m = importlib.util.module_from_spec(spec);"
        "spec.loader.exec_module(m);"
        "tracemalloc.start();"
        "n = m.count_pg(sys.argv[2]);"
        "print(tracemalloc.get_traced_memory()[1], n)"
    )
    results = {}
    for label, bam in (("small", small), ("large", large)):
        out = subprocess.run(
            [sys.executable, "-c", probe, str(Path(module.__file__)), str(bam)],
            capture_output=True, text=True, check=True,
        )
        peak, counted = out.stdout.split()
        results[label] = int(peak)
        assert int(counted) == (20 if label == "small" else 50000), "count must be right"

    # Streaming holds one 1 MB chunk plus a 3-byte overlap, so a header 2,500x
    # larger should cost a rounding error. Any materializing form allocates at
    # least 1x the header text, an order of magnitude above this bound at this
    # fixture size. Both directions verified against real implementations.
    excess = results["large"] - results["small"]
    budget = 5 * 1024 * 1024
    assert excess < budget, (
        "allocation grew by %.1f MB for a %.1f MB header (%.0fx size growth) -- "
        "the count is materializing the header instead of streaming it"
        % (excess / 1048576, header_bytes(large) / 1048576, growth)
    )


def test_collapse_helper_does_not_open_through_pysam():
    """The helper must not route the header through htslib.

    pysam/htslib materializes the whole header on open (measured 104.2 MB
    against an 81.3 MB header, 469.0 MB via to_dict), which would defeat the
    streaming reader. The tag scan is delegated to samtools deliberately: at
    1.11x the header it costs exactly what the reheader that follows it already
    costs, so it adds no peak the operation was not already paying.
    """
    source = (
        Path(__file__).resolve().parents[1] / "util" / "misc" / "collapse_bam_pg_header.py"
    ).read_text()
    code = "\n".join(
        line for line in source.splitlines() if not line.strip().startswith("#")
    )
    body = re.sub(r'"""(?:.|\n)*?"""', "", code)
    assert "import pysam" not in body, "helper must not import pysam"
    assert "AlignmentFile" not in body, "helper must not open BAMs through htslib"
    assert "view -H" not in body, "samtools view -H appends its own @PG records"
    assert '"-d"' in body and '"PG"' in body, (
        "tag scan should shell out to samtools -d PG, which at 1.11x the header "
        "is exactly what the reheader that follows already costs"
    )
