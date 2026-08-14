#!/usr/bin/env python3

"""Every in-repo tracking consumer must accept the gzipped tracking LRAA can now write.

LRAA writes quant.tracking gzipped, always. It is the only artifact that scales with library
size, so the compressed form is the only form that exists. A consumer that opens it with a bare
`open()` does not degrade gracefully: it reads the deflate header as text and dies, or worse, a
permissive parser yields garbage rows.

Two consumers were missed in the first pass (`write_isoform_unique_support_bams.py` and the
single-cell pseudobulk matrix builder) and were found by review, not by tests -- because no
test drove them with a compressed input at all. This pins the property for the whole set: the
same input, in both encodings, must produce the same output.

The R consumers are excluded deliberately: `read.csv` and `data.table::fread` both decompress
by magic-byte detection, verified directly rather than assumed.
"""

import gzip
import os
import subprocess
import sys

import pytest

REPO = os.path.dirname(os.path.dirname(os.path.realpath(__file__)))
UTIL = os.path.join(REPO, "util")

SC_CONSUMER = os.path.join(
    UTIL, "sc", "sc_expr_tracking_and_cell_clustering_to_cell_cluster_pseudobulk_count_matrix.py"
)
SPARSE_CONSUMER = os.path.join(UTIL, "sc", "singlecell_tracking_to_sparse_matrix.py")

HEADER = (
    "gene_id\ttranscript_id\ttranscript_splice_hash_code\tnum_exons\t"
    "mp_id\tread_name\tfrac_assigned\tread_weight\n"
)


def _rows():
    out = ["# LRAA version test\n", HEADER]
    for i in range(6):
        bc = f"BC{i % 3}"
        out.append(
            f"G{i % 2}\tT{i}\thash{i}\t2\tmp{i}\t{bc}^U{i}^r{i}\t1.000000\t1.000000\n"
        )
    return out


def _write_both(tmp_path):
    rows = _rows()
    plain = tmp_path / "sc.quant.tracking"
    gz = tmp_path / "sc.quant.tracking.gz"
    plain.write_text("".join(rows))
    with gzip.open(gz, "wt") as fh:
        fh.writelines(rows)
    return plain, gz


def _clusters(tmp_path):
    p = tmp_path / "clusters.tsv"
    p.write_text("cell\tcluster\n" + "".join(f"BC{i}\tclust{i % 2}\n" for i in range(3)))
    return p


def test_pseudobulk_matrix_builder_gives_identical_output_for_both_encodings(tmp_path):
    plain, gz = _write_both(tmp_path)
    clusters = _clusters(tmp_path)

    outs = {}
    for tag, tracking in (("plain", plain), ("gz", gz)):
        out = tmp_path / f"matrix.{tag}.tsv"
        r = subprocess.run(
            [sys.executable, SC_CONSUMER,
             "--cell_clusters_info", str(clusters),
             "--LRAA_tracking_file", str(tracking),
             "--output_matrix", str(out)],
            capture_output=True, text=True, cwd=str(tmp_path),
        )
        assert r.returncode == 0, (
            f"{tag} input failed:\n{r.stdout[-1500:]}\n{r.stderr[-1500:]}"
        )
        assert out.exists(), f"{tag} run produced no matrix"
        outs[tag] = out.read_text()

    assert outs["plain"], "fixture produced an empty matrix; the test proves nothing"
    assert outs["gz"] == outs["plain"], (
        "the pseudobulk matrix differs between plain and gzipped tracking input"
    )


def test_sparse_matrix_converter_gives_identical_output_for_both_encodings(tmp_path):
    """This one was already gzip-aware; pinned so it stays that way.

    Compares file CONTENTS, not just the set of names: a regression that emits correctly-named
    but empty or wrong matrices from the compressed path would satisfy a name-only check, and
    this converter's entire job is the contents.
    """
    import hashlib

    plain, gz = _write_both(tmp_path)

    digests = {}
    for tag, tracking in (("plain", plain), ("gz", gz)):
        outdir = tmp_path / f"out_{tag}"
        outdir.mkdir()
        r = subprocess.run(
            [sys.executable, SPARSE_CONSUMER,
             "--tracking", str(tracking),
             "--output_prefix", str(outdir / "sparse")],
            capture_output=True, text=True, cwd=str(tmp_path),
        )
        assert r.returncode == 0, (
            f"{tag} input failed:\n{r.stdout[-1500:]}\n{r.stderr[-1500:]}"
        )
        per_file = {}
        for root, _dirs, files in os.walk(outdir):
            for f in files:
                full = os.path.join(root, f)
                rel = os.path.relpath(full, outdir)
                # decompressed, because gzip embeds an mtime: two runs producing identical
                # content still differ byte-for-byte in the header
                if f.endswith(".gz"):
                    with gzip.open(full, "rb") as gh:
                        payload = gh.read()
                else:
                    payload = open(full, "rb").read()
                per_file[rel] = hashlib.md5(payload).hexdigest()
        digests[tag] = per_file

    assert digests["plain"], "fixture produced no sparse output; the test proves nothing"
    assert digests["gz"] == digests["plain"], (
        "sparse output differs between plain and gzipped tracking input:\n"
        f"  only in plain: {sorted(set(digests['plain']) - set(digests['gz']))}\n"
        f"  only in gz   : {sorted(set(digests['gz']) - set(digests['plain']))}\n"
        f"  differing    : {sorted(k for k in digests['plain'] if k in digests['gz'] and digests['plain'][k] != digests['gz'][k])}"
    )


def test_unique_support_bam_writer_gives_identical_output_for_both_encodings(tmp_path):
    """Driven with a real BAM, so it asserts success rather than the absence of error strings.

    The earlier version passed a nonexistent BAM and asserted that no decode error appeared.
    That was order-dependent: it only proved anything while the script happens to parse the
    tracking file before touching the BAM, so a refactor that opened the BAM first would make
    it green without the gzip path being exercised at all.
    """
    import pysam

    reads = ["r0", "r1", "r2", "r3", "r4", "r5"]
    bam = tmp_path / "reads.bam"
    header = {"HD": {"VN": "1.6", "SO": "coordinate"},
              "SQ": [{"SN": "chr1", "LN": 2000}]}
    with pysam.AlignmentFile(str(bam), "wb", header=header) as fh:
        for i, name in enumerate(reads):
            a = pysam.AlignedSegment()
            a.query_name = name
            a.query_sequence = "A" * 50
            a.flag = 0
            a.reference_id = 0
            a.reference_start = 100 + i
            a.mapping_quality = 60
            a.cigarstring = "50M"
            a.query_qualities = pysam.qualitystring_to_array("I" * 50)
            fh.write(a)
    pysam.index(str(bam))

    # Bulk-style read names, matching the bam's query names. The single-cell fixture encodes
    # barcode and UMI into read_name, which no bam record here carries, so every isoform would
    # come out with zero unique support and the comparison would hold vacuously.
    rows = ["# LRAA version test\n", HEADER] + [
        f"G{i % 2}\tT{i}\thash{i}\t2\tmp{i}\t{name}\t1.000000\t1.000000\n"
        for i, name in enumerate(reads)
    ]
    bulk_plain = tmp_path / "bulk.quant.tracking"
    bulk_gz = tmp_path / "bulk.quant.tracking.gz"
    bulk_plain.write_text("".join(rows))
    with gzip.open(bulk_gz, "wt") as fh:
        fh.writelines(rows)

    writer = os.path.join(UTIL, "misc", "write_isoform_unique_support_bams.py")
    listings = {}
    for tag, tracking in (("plain", bulk_plain), ("gz", bulk_gz)):
        workdir = tmp_path / f"w_{tag}"
        workdir.mkdir()
        r = subprocess.run(
            [sys.executable, writer, "--bam", str(bam), "--tracking", str(tracking)],
            capture_output=True, text=True, cwd=str(workdir),
        )
        assert r.returncode == 0, (
            f"{tag} input failed:\n{r.stdout[-1500:]}\n{r.stderr[-1500:]}"
        )
        listings[tag] = sorted(
            f for f in os.listdir(workdir) if f.endswith(".unique_reads.bam")
        )

    assert listings["plain"], (
        "no per-isoform bams were written; the fixture has no uniquely assigned reads and the "
        "test would pass no matter how the tracking file were read"
    )
    assert listings["gz"] == listings["plain"], (
        f"different bams from the two encodings: {listings}"
    )
