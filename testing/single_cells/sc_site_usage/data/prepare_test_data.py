#!/usr/bin/env python3

"""How the inputs of this test were cut from the PBMC Kinnex single-cell data of the
LRAA paper (ref-guided, cluster-guided LRAA v0.44.0 run). Kept for provenance: it needs
the full data set and is not run by the test.

Seven genes chosen for clear, seed-stable site switches in the full analysis:
  EMP3, SELENOH, CIAO2A   alternative TSS, same splicing (tandem TSSs)
  AIF1                    alternative TSS on a different first exon
  POLR2K, CMPK1           alternative PolyA, same splicing (tandem 3' UTR)
  ELOVL5                  alternative PolyA, intronic site (different last exon)

Inputs written to this directory, all restricted to the genes' spans (+/- 1 kb) and to a
--cell_frac of the clustered cells (all, by default):
  cell_clusters.tsv                 cell barcode, cluster
  cluster_cell_types.tsv            cluster, majority cell type (for figure labels)
  reads.bam(.bai)                   primary alignments of those cells; base qualities and
                                    all tags but CB, XM, ts, NM dropped (not used here)
  sites.TSS.bed, sites.PolyA.bed    LRAA integrated site beds
  models.gtf                        LRAA models (SYMBOL^ transcript ids)
  gene_trans_map.tsv                transcript id -> symbol-prefixed ids
  cluster_quant.tar.gz              per-cluster quant.expr, these models only
  tracking.tsv.gz                   per-cluster read assignments, these models and cells
"""

import argparse
import csv
import gzip
import io
import os
import random
import re
import subprocess
import tarfile

import pysam

GENES = ["EMP3", "SELENOH", "CIAO2A", "AIF1", "POLR2K", "CMPK1", "ELOVL5"]
KEEP_TAGS = {"CB", "XM", "ts", "NM"}


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--eval_dir", default="/home/unix/bhaas/projects/LRAA_PAPER_Analyses/PBMCs_kinnex/LRAA_PBMCs_eval")
    p.add_argument("--cell_frac", type=float, default=1.0)
    p.add_argument("--flank", type=int, default=1000)
    p.add_argument("--seed", type=int, default=1)
    p.add_argument("--outdir", default=os.path.dirname(os.path.abspath(__file__)))
    a = p.parse_args()

    E = a.eval_dir
    cg = f"{E}/pbmcs_refGuided/pbmcs_refGuided_sc_cluster_guided/data"
    su = f"{E}/__Site_Usage_Analysis"
    out = lambda f: os.path.join(a.outdir, f)

    # cells
    rows = [l.rstrip("\n").split("\t") for l in open(f"{E}/pbmcs_refGuided/data/RefQuantOnly.cell_clusters.tsv")]
    rows = [r for r in rows if len(r) >= 2 and r[1].strip().isdigit()]
    rng = random.Random(a.seed)
    cells = {r[0]: r[1] for r in rows if rng.random() < a.cell_frac}
    with open(out("cell_clusters.tsv"), "wt") as ofh:
        print("cell_barcode\tcluster", file=ofh)
        for cb, cl in cells.items():
            print(f"{cb}\t{cl}", file=ofh)

    import pandas as pd
    u = pd.read_csv(f"{E}/pbmcs_refGuided/data/RefQuantOnly.cell_clusters_and_cell_types.wUMAP.tsv", sep="\t")
    t = u.groupby("seurat_clusters").cell_type_annot.agg(lambda x: x.value_counts().index[0]).reset_index()
    t.columns = ["cluster", "cell_type"]
    t.to_csv(out("cluster_cell_types.tsv"), sep="\t", index=False)

    # gene regions from the site analysis' gene spans
    regions = []
    for r in csv.DictReader(open(f"{su}/PBMCs.gene_spans.tsv"), delimiter="\t"):
        if r["gene_symbol"] in GENES:
            regions.append((r["chrom"], max(1, int(r["start"]) - a.flank), int(r["end"]) + a.flank))
    regions.sort()
    in_region = lambda c, s, e: any(c == rc and s <= re_ and e >= rs for rc, rs, re_ in regions)

    # reads
    src = pysam.AlignmentFile(f"{E}/../PBMCs_pbio.aligned.sorted.bam")
    tmp = out("reads.unsorted.bam")
    seen = set()
    with pysam.AlignmentFile(tmp, "wb", template=src) as dst:
        for c, s, e in regions:
            for r in src.fetch(c, s - 1, e):
                if r.is_secondary or r.is_supplementary or not r.has_tag("CB") or r.get_tag("CB") not in cells:
                    continue
                key = (r.query_name, r.reference_id, r.reference_start)
                if key in seen:
                    continue
                seen.add(key)
                r.set_tags([(t, v) for t, v in r.get_tags() if t in KEEP_TAGS])
                r.query_qualities = None
                dst.write(r)
    pysam.sort("-o", out("reads.bam"), tmp)
    os.remove(tmp)
    pysam.index(out("reads.bam"))

    # site beds
    for kind in ("TSS", "PolyA"):
        with open(out(f"sites.{kind}.bed"), "wt") as ofh:
            for line in open(f"{su}/site_beds/PBMCs_pbio_CG_refguided.integrated.{kind}.bed"):
                if line.startswith("#"):
                    ofh.write(line)
                    continue
                f = line.split("\t")
                if in_region(f[0], int(f[2]), int(f[2])):
                    ofh.write(line)

    # models overlapping the regions
    tids = set()
    gtf = f"{cg}/PBMCs_pbio_CG_refguided.withGeneSymbols.gtf"
    for line in open(gtf):
        f = line.split("\t")
        if len(f) > 8 and f[2] == "transcript" and in_region(f[0], int(f[3]), int(f[4])):
            tids.add(re.search(r'transcript_id "([^"]+)"', f[8]).group(1))
    bare = {t.split("^", 1)[-1] for t in tids}
    with open(out("models.gtf"), "wt") as ofh:
        for line in open(gtf):
            m = re.search(r'transcript_id "([^"]+)"', line)
            if m and m.group(1) in tids:
                ofh.write(line)

    with open(f"{cg}/PBMCs_pbio_CG_refguided.gene_transcript_splicehashcode.withGeneSymbols.tsv") as fh, \
            open(out("gene_trans_map.tsv"), "wt") as ofh:
        ofh.write(fh.readline())
        for line in fh:
            if line.split("\t")[1] in bare:
                ofh.write(line)

    # per-cluster quant.expr, these models only
    with tarfile.open(f"{cg}/PBMCs_pbio_CG_refguided.LRAA.final.cluster_quant.EXPRs.tar.gz") as src_tar, \
            tarfile.open(out("cluster_quant.tar.gz"), "w:gz") as dst_tar:
        for mem in src_tar.getmembers():
            if not mem.name.endswith("quant.expr"):
                continue
            lines = io.TextIOWrapper(src_tar.extractfile(mem)).read().splitlines(keepends=True)
            keep = [l for l in lines if l.startswith("#") or l.startswith("gene_id") or l.split("\t")[1] in bare]
            data = "".join(keep).encode()
            info = tarfile.TarInfo(os.path.basename(mem.name))
            info.size = len(data)
            dst_tar.addfile(info, io.BytesIO(data))

    # read assignments, these models and cells
    with gzip.open(f"{cg}/PBMCs_pbio.cluster_merged.quant.tracking.gz", "rt") as fh, \
            gzip.open(out("tracking.tsv.gz"), "wt") as ofh:
        for line in fh:
            if line.startswith("#") or line.startswith("gene_id"):
                ofh.write(line)
                continue
            f = line.split("\t", 6)
            if f[1].split("@")[-1] in bare and f[5].split("^", 1)[0] in cells:
                ofh.write(line)

    print(f"{len(cells)} cells, {len(seen)} alignments, {len(tids)} models")


if __name__ == "__main__":
    main()
