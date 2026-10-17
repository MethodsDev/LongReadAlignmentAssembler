#!/usr/bin/env python3

"""The merge must keep an alternative terminal exon that several clusters rebuilt.

A first (last) exon that starts (ends) inside another model's intron is not a partial
model: it carries intronic sequence the longer model splices out. The merge leaves
fracture_splice_graph_at_input_transcript_bounds off, so with no TSS/PolyA call at that
end the merge graph had no segment boundary there; the segment spanning the intron was
under-covered, _map_segment_to_graph_INITIAL dropped it (min_feature_frac_overlap), and the
truncated path was absorbed as contained. MEASURED on the LRAA v0.44.3 LP00714 (ALS)
ref-guided cluster-guided run: 701 basic FSM chains with a unique FSM read were in the
per-cluster GTFs and missing from the merged one, GPI's canonical ENST00000356487 chain
among them (13 of 14 clusters).

The fix fractures the merge graph at such a terminus when at least
merge_overhang_min_sources source GTFs carry the chain with an end in one alt-site window.
A model ending inside the other model's exon is still a partial model and still absorbed.

The pooled (all-cells) GTF, passed as --pooled_gtf, is a structure source in its own
right: its alternative terminal exons fracture the merge graph without needing a second
source, and its TSS/PolyA read counts are counted only beyond the per-cluster support
around a site, since the two were built from largely the same reads.
"""

import os
import random
import re
import subprocess
import sys

REPO_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
MERGE_SCRIPT = os.path.join(REPO_ROOT, "util", "merge_LRAA_GTFs.py")

CONTIG_LEN = 4000
# container: four exons; its first intron is 1201-1499
CONTAINER = [(1000, 1200), (1500, 1700), (2000, 2300), (2600, 2800)]
# alternative first exon: starts 100 bp inside the container's first intron (100/299 of
# that intron's segment covered, under the 0.5 mapping threshold)
ALT_FIRST = [(1400, 1700), (2000, 2300), (2600, 2800)]
# partial model: starts inside the container's second exon -- a truncation, not an alt exon
PARTIAL = [(1600, 1700), (2000, 2300), (2600, 2800)]
# starts at the container's first intronic base and ends an exon early, as GPI's 34,365,020
# model does: it covers the intronic stretch, so the merge graph has one exon segment
# spanning the whole intron, and with a shorter intron chain it contains no other model
EARLY = [(1201, 1700), (2000, 2300)]


def _genome(tmp_path):
    rng = random.Random(7)
    seq = [rng.choice("ACGT") for _ in range(CONTIG_LEN)]
    for donor, acceptor in ((1201, 1499), (1701, 1999), (2301, 2599)):  # GT...AG on the + strand
        seq[donor - 1 : donor + 1] = list("GT")
        seq[acceptor - 2 : acceptor] = list("AG")
    fasta = tmp_path / "genome.fa"
    body = "".join(seq)
    with open(fasta, "wt") as fh:
        fh.write(">chr1\n")
        for i in range(0, len(body), 60):
            fh.write(body[i : i + 60] + "\n")
    return str(fasta)


def _gtf(path, models):
    """models: list of (transcript_id, exons[, TSS_read_count]). No PolyA call on any
    end; a TSS call only where a TSS_read_count is given."""
    rows = []
    for tid, exons, *tss in models:
        if tss:
            tss_attrs = f'TSS "True"; TSS_read_count "{tss[0]}";'
        else:
            tss_attrs = 'TSS "False";'
        attrs = f'gene_id "g1"; transcript_id "{tid}"; {tss_attrs} PolyA "False";'
        rows.append("\t".join(["chr1", "LRAA", "transcript", str(exons[0][0]), str(exons[-1][1]), ".", "+", ".", attrs]))
        for lend, rend in exons:
            rows.append("\t".join(["chr1", "LRAA", "exon", str(lend), str(rend), ".", "+", ".", f'gene_id "g1"; transcript_id "{tid}";']))
    with open(path, "wt") as fh:
        fh.write("\n".join(rows) + "\n")
    return str(path)


def _merge_records(tmp_path, gtfs, pooled=None):
    """(exons, transcript attributes) per merged model."""
    out = str(tmp_path / "merged.gtf")
    cmd = [sys.executable, MERGE_SCRIPT, "--genome", _genome(tmp_path), "--CPU", "1", "--gtf"]
    cmd += gtfs + ["--output_gtf", out]
    if pooled is not None:
        cmd += ["--pooled_gtf", pooled]
    result = subprocess.run(cmd, capture_output=True, text=True, cwd=str(tmp_path))
    assert result.returncode == 0, result.stderr[-3000:]
    exons, attrs = {}, {}
    for line in open(out):
        cols = line.rstrip("\n").split("\t")
        if line.startswith("#") or len(cols) < 9:
            continue
        tid = re.search(r'transcript_id "([^"]+)"', cols[8]).group(1)
        if cols[2] == "transcript":
            attrs[tid] = cols[8]
        elif cols[2] == "exon":
            exons.setdefault(tid, []).append((int(cols[3]), int(cols[4])))
    return [(sorted(exons[t]), attrs.get(t, "")) for t in exons]


def _merge(tmp_path, gtfs, pooled=None):
    return [e for e, _ in _merge_records(tmp_path, gtfs, pooled)]


def _chain(exons):
    return tuple((exons[i][1] + 1, exons[i + 1][0] - 1) for i in range(len(exons) - 1))


ALT_CHAIN = _chain(ALT_FIRST)  # PARTIAL shares it; the two are told apart by start


def _starting_in(models, lo, hi):
    return [m for m in models if _chain(m) == ALT_CHAIN and lo <= m[0][0] <= hi]


def test_alt_first_exon_from_two_clusters_is_kept(tmp_path):
    gtfs = [
        _gtf(tmp_path / "c0.gtf", [("c0:container", CONTAINER), ("c0:early", EARLY)]),
        _gtf(tmp_path / "c1.gtf", [("c1:container", CONTAINER), ("c1:alt", ALT_FIRST)]),
        _gtf(tmp_path / "c2.gtf", [("c2:alt", ALT_FIRST)]),
    ]
    models = _merge(tmp_path, gtfs)
    alt = _starting_in(models, 1300, 1499)
    assert len(alt) == 1, models
    assert alt[0][0][0] == 1400, alt
    assert any(_chain(m) == _chain(CONTAINER) for m in models), models


def test_alt_first_exon_from_one_cluster_is_absorbed(tmp_path):
    # below merge_overhang_min_sources (2): the pre-fix behaviour stands
    gtfs = [
        _gtf(tmp_path / "c0.gtf", [("c0:container", CONTAINER), ("c0:early", EARLY)]),
        _gtf(tmp_path / "c1.gtf", [("c1:alt", ALT_FIRST)]),
    ]
    models = _merge(tmp_path, gtfs)
    assert not _starting_in(models, 1300, 1499), models


def test_partial_model_is_still_absorbed(tmp_path):
    # ends inside the container's exon, so no overhang: still a contained partial
    gtfs = [
        _gtf(tmp_path / "c0.gtf", [("c0:container", CONTAINER), ("c0:early", EARLY)]),
        _gtf(tmp_path / "c1.gtf", [("c1:partial", PARTIAL)]),
        _gtf(tmp_path / "c2.gtf", [("c2:partial", PARTIAL)]),
    ]
    models = _merge(tmp_path, gtfs)
    assert not _starting_in(models, 1500, 1700), models


def test_nearby_alt_ends_collapse_to_one_model(tmp_path):
    # three clusters, ends within a few bp of each other: one model, at the densest end
    gtfs = [_gtf(tmp_path / "c0.gtf", [("c0:container", CONTAINER), ("c0:early", EARLY)])]
    for i, start in enumerate((1396, 1400, 1404), start=1):
        gtfs.append(_gtf(tmp_path / f"c{i}.gtf", [(f"c{i}:alt", [(start, 1700), (2000, 2300), (2600, 2800)])]))
    models = _merge(tmp_path, gtfs)
    alt = _starting_in(models, 1300, 1499)
    assert len(alt) == 1, models
    assert 1396 <= alt[0][0][0] <= 1404, alt


def test_pooled_alt_first_exon_is_kept_without_a_second_source(tmp_path):
    # one source is enough when it is the pooled run: it resolved that end at full depth
    gtfs = [
        _gtf(tmp_path / "c0.gtf", [("c0:container", CONTAINER), ("c0:early", EARLY)]),
        _gtf(tmp_path / "c1.gtf", [("c1:container", CONTAINER)]),
    ]
    pooled = _gtf(tmp_path / "pooled.gtf", [("p:container", CONTAINER), ("p:alt", ALT_FIRST)])
    alt = _starting_in(_merge(tmp_path, gtfs, pooled), 1300, 1499)
    assert len(alt) == 1 and alt[0][0][0] == 1400, alt


def test_pooled_and_one_cluster_alt_end_is_not_two_votes_but_is_kept(tmp_path):
    # kept because the pooled model has it, with the end at an asserted coordinate
    gtfs = [
        _gtf(tmp_path / "c0.gtf", [("c0:container", CONTAINER), ("c0:early", EARLY)]),
        _gtf(tmp_path / "c1.gtf", [("c1:alt", [(1402, 1700), (2000, 2300), (2600, 2800)])]),
    ]
    pooled = _gtf(tmp_path / "pooled.gtf", [("p:alt", ALT_FIRST)])
    alt = _starting_in(_merge(tmp_path, gtfs, pooled), 1300, 1499)
    assert len(alt) == 1 and alt[0][0][0] in (1400, 1402), alt


def _tss_read_count(attrs):
    m = re.search(r'TSS_read_count "(\d+)"', attrs)
    return int(m.group(1)) if m else None


def _container_record(records):
    hits = [(e, a) for e, a in records if _chain(e) == _chain(CONTAINER)]
    assert len(hits) == 1, records
    return hits[0]


def test_pooled_tss_support_is_not_summed_with_the_clusters(tmp_path):
    # clusters: 10 + 10 reads at the TSS; pooled: 25 reads two bases away. Largely the
    # same reads, so the merged site carries 25 (the larger), not 45 (the sum).
    tss_model = [(1002, 1200)] + CONTAINER[1:]
    gtfs = [
        _gtf(tmp_path / "c1.gtf", [("c1:container", CONTAINER, 10)]),
        _gtf(tmp_path / "c2.gtf", [("c2:container", CONTAINER, 10)]),
    ]
    pooled = _gtf(tmp_path / "pooled.gtf", [("p:container", tss_model, 25)])
    exons, attrs = _container_record(_merge_records(tmp_path, gtfs, pooled))
    assert 'TSS "True"' in attrs, attrs
    assert _tss_read_count(attrs) == 25, attrs


def test_pooled_tss_below_cluster_support_adds_nothing(tmp_path):
    gtfs = [
        _gtf(tmp_path / "c1.gtf", [("c1:container", CONTAINER, 30)]),
        _gtf(tmp_path / "c2.gtf", [("c2:container", CONTAINER, 10)]),
    ]
    pooled = _gtf(tmp_path / "pooled.gtf", [("p:container", CONTAINER, 35)])
    exons, attrs = _container_record(_merge_records(tmp_path, gtfs, pooled))
    assert _tss_read_count(attrs) == 40, attrs
    assert exons[0][0] == 1000, exons


def test_pooled_tss_no_cluster_called_is_carried_in_full(tmp_path):
    gtfs = [
        _gtf(tmp_path / "c1.gtf", [("c1:container", CONTAINER)]),
        _gtf(tmp_path / "c2.gtf", [("c2:container", CONTAINER)]),
    ]
    pooled = _gtf(tmp_path / "pooled.gtf", [("p:container", CONTAINER, 12)])
    exons, attrs = _container_record(_merge_records(tmp_path, gtfs, pooled))
    assert 'TSS "True"' in attrs, attrs
    assert _tss_read_count(attrs) == 12, attrs


def test_tssless_end_within_full_window_of_a_tss_collapses_onto_it(tmp_path):
    # a TSS-anchored model at 1000 and the same chain with no call starting 40 bp
    # upstream (outside half the 50 bp alt-TSS window, inside the whole of it): one model
    upstream = [(960, 1200)] + CONTAINER[1:]
    gtfs = [
        _gtf(tmp_path / "c1.gtf", [("c1:container", CONTAINER, 10)]),
        _gtf(tmp_path / "c2.gtf", [("c2:upstream", upstream)]),
    ]
    hits = [(e, a) for e, a in _merge_records(tmp_path, gtfs) if _chain(e) == _chain(CONTAINER)]
    assert len(hits) == 1, hits
    assert hits[0][0][0][0] == 1000 and 'TSS "True"' in hits[0][1], hits


def test_tssless_end_beyond_the_window_stays_distinct(tmp_path):
    upstream = [(900, 1200)] + CONTAINER[1:]
    gtfs = [
        _gtf(tmp_path / "c1.gtf", [("c1:container", CONTAINER, 10)]),
        _gtf(tmp_path / "c2.gtf", [("c2:upstream", upstream)]),
    ]
    starts = sorted(e[0][0] for e, _ in _merge_records(tmp_path, gtfs) if _chain(e) == _chain(CONTAINER))
    assert starts == [900, 1000], starts
