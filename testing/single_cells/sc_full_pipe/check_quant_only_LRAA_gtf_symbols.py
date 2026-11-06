#!/usr/bin/env python3
# Checks a basic quant_only run quantified against an LRAA-built GTF (LRAA ids, not
# reference ids): gffcompare must have run so the LRAA features get gene names, and
# splice-pattern features must keep their hash codes whether named or not.

import glob, gzip, os, sys, tarfile

run_dir = sorted(glob.glob(os.path.join(sys.argv[1], "*_LRAA_singlecell_wf")))[-1]
out = os.path.join(run_dir, "out")

tracking = glob.glob(os.path.join(out, "incl_gene_symbols_gffcompare_tracking", "*"))
assert tracking and os.path.getsize(tracking[0]) > 0, "no gffcompare tracking: the input GTF was not compared to the reference"


def features(kind):
    tar_path = glob.glob(os.path.join(out, f"incl_gene_symbols_{kind}_sparse_tar_gz", "*"))[0]
    with tarfile.open(tar_path) as tar:
        read = lambda name: [l.decode().rstrip("\n") for l in gzip.open(tar.extractfile(
            next(m for m in tar.getmembers() if m.name.endswith("/" + name))))]
        return read("features.tsv.gz.orig"), read("features.tsv.gz")


orig, named = features("isoform")
lraa_named = sum(1 for o, n in zip(orig, named) if o.startswith("t:") and "^" in n)
assert lraa_named > 0, "no LRAA isoform got a gene name"
print(f"isoforms: {lraa_named} of {sum(o.startswith('t:') for o in orig)} LRAA ids named")

orig, named = features("splice_pattern")
lost = [(o, n) for o, n in zip(orig, named) if ":iso-" not in o and n.split("^", 1)[-1] != o]
assert not lost, f"{len(lost)} splice hash codes replaced, e.g. {lost[:3]}"
print(f"splice patterns: {sum('^' in n for n in named)} of {len(named)} named, every hash code kept")
