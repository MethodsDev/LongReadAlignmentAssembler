#!/usr/bin/env python3

"""A merged catalog must carry PAS and internal-priming calls for its own 3' ends.

`reconstruct_isoforms` builds new Transcript objects out of multipaths, so nothing an
input GTF carried in `_meta` reaches the merged model. Boundary READ COUNTS survive --
MultiPath copies them off the graph's boundary nodes -- which made the loss easy to miss:
the merged models looked annotated. MEASURED on the v0.37.0 single-cell cluster-guided
test run, per-cluster inputs carried PAS on 192 of 192 transcripts and the merged catalog
on 0 of 860. Single-cell catalogs are produced by this script, so every downstream
consumer of a merged GTF -- including the TSS/PolyA site beds -- had no polyadenylation
signal and no internal-priming call at all.

Recomputed from the genome rather than copied from the inputs, because both describe the
sequence around the model's OWN 3' terminus and a merged model's terminus need not
coincide with any single input's. Annotation only: a merge reconciles catalogs, and the
source runs already applied their own filtering policy.
"""

import importlib.machinery
import importlib.util
import os
import re
import subprocess
import sys

import pytest

REPO_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if os.path.join(REPO_ROOT, "pylib") not in sys.path:
    sys.path.insert(0, os.path.join(REPO_ROOT, "pylib"))

MERGE_SCRIPT = os.path.join(REPO_ROOT, "util", "merge_LRAA_GTFs.py")


# AATAAA occupies 1-based 1674..1679, and the + strand 3' end is 1700.  find_polyA_signal
# reports the signed distance to the hexamer's FIRST base, so the expected offset is
# 1674 - 1700 = -26, inside the default [-40, -10] window.  The A run after 1700 is what
# the internal-priming check reads.
PAS = "AATAAA"
PAS_OFFSET = "-26"
CONTIG_LEN = 2000


def _genome(tmp_path):
    seq = ["T"] * CONTIG_LEN
    seq[1673:1679] = list(PAS)  # 1-based 1674..1679, i.e. offset -21 from 1700
    seq[1700:1730] = list("A" * 30)  # oligo-dT template downstream of the 3' end
    fasta = tmp_path / "genome.fa"
    body = "".join(seq)
    with open(fasta, "wt") as fh:
        fh.write(">chr1\n")
        for i in range(0, len(body), 60):
            fh.write(body[i : i + 60] + "\n")
    return str(fasta)


def _input_gtf(path, transcript_id, exons):
    """A two-exon model with boundary attributes but NO PAS metadata at all.

    Deliberate: the merged model must be annotated from the genome, so an input that
    never carried the attribute still produces one that does.
    """
    attrs = 'gene_id "g1"; transcript_id "{}"; TSS "True"; PolyA "True"; ' 'TSS_read_count "40"; PolyA_read_count "25";'.format(
        transcript_id
    )
    exon_attrs = 'gene_id "g1"; transcript_id "{}";'.format(transcript_id)
    rows = [
        "\t".join(
            [
                "chr1",
                "LRAA",
                "transcript",
                str(exons[0][0]),
                str(exons[-1][1]),
                ".",
                "+",
                ".",
                attrs,
            ]
        )
    ]
    for lend, rend in exons:
        rows.append(
            "\t".join(
                ["chr1", "LRAA", "exon", str(lend), str(rend), ".", "+", ".", exon_attrs]
            )
        )
    with open(path, "wt") as fh:
        fh.write("\n".join(rows) + "\n")
    return str(path)


def _transcript_attrs(gtf_path):
    found = {}
    for line in open(gtf_path):
        if line.startswith("#"):
            continue
        cols = line.rstrip("\n").split("\t")
        if len(cols) < 9 or cols[2] != "transcript":
            continue
        attrs = dict(re.findall(r'(\w+) "([^"]*)"', cols[8]))
        found[attrs["transcript_id"]] = attrs
    return found


def test_merged_models_carry_PAS_and_internal_priming(tmp_path):
    genome = _genome(tmp_path)
    gtfs = [
        _input_gtf(tmp_path / "a.gtf", "t:a", [(1000, 1200), (1500, 1700)]),
        _input_gtf(tmp_path / "b.gtf", "t:b", [(1050, 1200), (1500, 1700)]),
    ]
    merged = str(tmp_path / "merged.gtf")

    result = subprocess.run(
        [sys.executable, MERGE_SCRIPT, "--genome", genome, "--gtf"]
        + gtfs
        + ["--output_gtf", merged],
        capture_output=True,
        text=True,
        cwd=str(tmp_path),
    )
    assert result.returncode == 0, result.stderr[-3000:]

    models = _transcript_attrs(merged)
    assert models, "merge produced no transcripts"

    for transcript_id, attrs in models.items():
        # The motif is real, found upstream of this model's own 3' end.
        assert attrs["PAS"] == PAS, (transcript_id, attrs)
        assert attrs["PAS_offset"] == PAS_OFFSET, (transcript_id, attrs)
        # 30 genomic A's follow the terminus, so the artifact call must fire.
        assert attrs["InternalPriming"] == "True", (transcript_id, attrs)


def test_annotation_does_not_delete_models(tmp_path):
    """The A-rich 3' end above is exactly what the read pass would filter.

    A merge must not apply that policy: the sources already decided, and dropping a
    model here would remove it from every cluster's catalog at once.
    """
    genome = _genome(tmp_path)
    gtfs = [
        _input_gtf(tmp_path / "a.gtf", "t:a", [(1000, 1200), (1500, 1700)]),
        _input_gtf(tmp_path / "b.gtf", "t:b", [(1000, 1200), (1600, 1700)]),
    ]
    merged = str(tmp_path / "merged.gtf")

    result = subprocess.run(
        [sys.executable, MERGE_SCRIPT, "--genome", genome, "--gtf"]
        + gtfs
        + ["--output_gtf", merged],
        capture_output=True,
        text=True,
        cwd=str(tmp_path),
    )
    assert result.returncode == 0, result.stderr[-3000:]

    models = _transcript_attrs(merged)
    # Two distinct intron chains in, both flagged internally primed, both still out.
    assert len(models) == 2, models
    assert all(attrs["InternalPriming"] == "True" for attrs in models.values())
