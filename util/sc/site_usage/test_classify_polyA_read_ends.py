"""pytest: classify_polyA_read_ends.py on a synthetic contig"""

import gzip
import os
import subprocess
import sys

import pandas as pd
import pysam

HERE = os.path.dirname(os.path.realpath(__file__))
SCRIPT = os.path.join(HERE, "classify_polyA_read_ends.py")


def _contig():
    seq = ["C"] * 2000
    # 3' end at 300 (+): PAS 22 nt upstream, A-rich downstream -> PAS rescues it (LRAA_site)
    seq[300 - 22 - 1:300 - 22 - 1 + 6] = list("AATAAA")
    seq[300:320] = ["A"] * 20
    # 3' end at 600 (+): no PAS, 20 A's downstream -> IP-like (no site)
    seq[600:620] = ["A"] * 20
    # 3' end at 900 (+): no PAS, A-rich -> IP-like at an LRAA site
    seq[900:920] = ["A"] * 20
    # 3' end at 1500 (-): 20 T's genomically upstream (= A's past the end in transcript sense)
    seq[1500 - 21:1500 - 1] = ["T"] * 20
    return "".join(seq)


def test_categories(tmp_path):
    fa = tmp_path / "g.fa"
    fa.write_text(">chrT\n" + _contig() + "\n")
    pysam.faidx(str(fa))
    gtf = tmp_path / "ref.gtf"
    gtf.write_text("chrT\tsrc\ttranscript\t1\t100\t.\t+\t.\ttranscript_id \"t1\";\n"
                   "chrT\tsrc\ttranscript\t1700\t1900\t.\t-\t.\ttranscript_id \"t2\";\n")
    atlas = tmp_path / "atlas.bed.gz"
    with gzip.open(atlas, "wt") as fh:
        fh.write("T\t1195\t1205\tc1\t1\t+\t0.5\n")  # 1-based 1196-1205; chr prefix added
    bed = tmp_path / "lraa.bed"
    bed.write_text("chrT\t299\t300\tPolyA:chrT:300:+\t0\t+\t10\n"
                   "chrT\t909\t910\tPolyA:chrT:910:+\t0\t+\t10\n")
    hist = tmp_path / "hist.tsv.gz"
    rows = [("+", 105, 0, 5),    # known: ref end 100 within 25
            ("+", 1220, 1, 3),   # known: PolyASite cluster ends 1205, within 25
            ("+", 300, 0, 7),    # LRAA site with PAS
            ("+", 600, 2, 11),   # IP-like, no site
            ("+", 900, 0, 13),   # LRAA site 910, IP-like
            ("-", 1500, 0, 17),  # IP-like on '-'
            ("+", 1800, 3, 19),  # other (soft clip > 2 kept as its own bin)
            ("+", 105, 0, 1)]    # repeated position: summed
    with gzip.open(hist, "wt") as fh:
        fh.write("chrom\tstrand\tpos\tresidual_soft_clip\treads\n")
        for s, p, c, n in rows:
            fh.write(f"chrT\t{s}\t{p}\t{c}\t{n}\n")
    prefix = tmp_path / "out"
    subprocess.run([sys.executable, SCRIPT, "--histogram", str(hist), "--ref_gtf", str(gtf),
                    "--polyasite_atlas", str(atlas), "--lraa_polyA_bed", str(bed), "--genome", str(fa),
                    "--tolerance", "25", "--output_prefix", str(prefix)], check=True)
    cat = pd.read_csv(f"{prefix}.categories.tsv", sep="\t")
    got = {(r.residual_soft_clip, r.category): r.reads for r in cat.itertuples()}
    assert got == {(0, "known"): 6, (1, "known"): 3, (0, "LRAA_site"): 7, (2, "IP_like"): 11,
                   (0, "LRAA_site_IP_like"): 13, (0, "IP_like"): 17, (3, "other"): 19}
    ev = pd.read_csv(f"{prefix}.evidence.tsv", sep="\t")
    assert ev.loc[ev.category == "LRAA_site", "pas"].all()
    assert set(ev.loc[ev.category == "IP_like", "downstream_A"]) == {20}
