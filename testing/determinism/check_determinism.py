#!/usr/bin/env python3
"""Assert that repeated LRAA runs on identical input produce identical output.

Two distinct properties are checked, because a fix for one does not imply the
other and each has been broken independently in this codebase:

  content   the canonical, sorted content of every principal output is identical.
            This is the guard against order-dependent isoform selection --
            unhashable objects iterating in memory-address order, a
            tie-break-free sort resolving ties by that order, and float sums
            accumulated over an unordered container.

  bytes     the raw files are identical.  This is the guard against the merge
            emitting the same content in a different order.  `random.shuffle`
            over the job list once reordered the concatenation of per-job
            outputs, so content-identical runs produced different bytes and any
            md5-based regression check was measuring the RNG rather than LRAA.

A content failure is a correctness bug.  A bytes-only failure is an ordering bug
in the merge; both are reported, and both fail the test.

Runs must be fresh processes with independent output directories AND independent
caches: sharing a preprocessing cache hides differences that arise during
preprocessing, and reusing an output prefix makes a stale file look like
agreement.

Choosing a region matters more than it appears.  A region whose run-to-run
differences are pure permutation will pass a content check while genuinely
nondeterministic code sits underneath it -- observed directly on chr20:30-34Mb,
where sorted output was identical and raw bytes were not, against
chr20:36,000,000-37,000,000 where the differences are real content.  Validate
any determinism work on a region known to carry content churn; the bundled
dataset here is small and exercises fewer code paths than a real contig, so it
is a smoke test rather than a strong guard.  Point --genome/--bam/--region at a
churn-carrying region for the strong version.
"""
import argparse
import gzip
import hashlib
import os
import shutil
import subprocess
import sys

PRINCIPAL_SUFFIXES = (
    ".gtf",
    ".bed",
    ".quant.expr",
    # Tracking is always written gzipped. Naming the uncompressed form here would find no
    # file and silently shrink the comparison to the remaining outputs -- a determinism
    # check that quietly stops checking the output most sensitive to merge ordering.
    ".quant.tracking.gz",
    ".pre-cross-gene-EM.quant.expr",
    ".pre-cross-gene-EM.quant.tracking.gz",
    # Pure per-job concatenation, so it is the output most sensitive to merge
    # ordering -- it was the single file that --no_shuffle_parallel_jobs was
    # shown to stabilise while content elsewhere still varied.  Omitting it
    # would drop the most direct witness of the bug this test guards.
    ".genome_tx_arb.summary.tsv",
)

# LRAA stamps its own argv into a leading comment.  Two runs invoked from
# different paths or prefixes differ on that line while being otherwise
# identical, so it is excluded from the content comparison and reported
# separately rather than silently ignored.
COMMENT_PREFIXES = ("# LRAA version", "# LRAA CMD:", "# LRAA merge:")


def _read_all(path):
    """Raw bytes, decompressing when the name says gzip.

    Both digests below go through this. A gzip member embeds the compression timestamp in its
    header, so two runs producing byte-identical content still produce different FILE bytes --
    and a determinism check that compared those would report a difference on every single run,
    for every gzipped output, forever. Comparing the decompressed payload is the only reading
    of "identical" that means anything here.
    """
    if str(path).endswith(".gz"):
        with gzip.open(path, "rb") as fh:
            return fh.read()
    with open(path, "rb") as fh:
        return fh.read()


def digest_bytes(path):
    """sha256 of the payload -- decompressed for gzipped outputs, raw otherwise."""
    return hashlib.sha256(_read_all(path)).hexdigest()


def digest_content(path):
    """sha256 of the sorted, comment-stripped lines."""
    lines = [
        ln for ln in _read_all(path).splitlines()
        if not any(ln.startswith(p.encode()) for p in COMMENT_PREFIXES)
    ]
    lines.sort()
    h = hashlib.sha256()
    for ln in lines:
        h.update(ln)
        h.update(b"\n")
    return h.hexdigest()


def run_once(args, replicate, workroot):
    outdir = os.path.join(workroot, f"rep{replicate}")
    os.makedirs(outdir, exist_ok=True)
    cmd = [
        sys.executable, args.lraa,
        "--genome", os.path.abspath(args.genome),
        "--bam", os.path.abspath(args.bam),
        "--output_prefix", args.output_prefix,
        "--num_threads_per_worker", str(args.threads),
    ]
    if args.gtf:
        cmd += ["--gtf", os.path.abspath(args.gtf)]
    if args.region:
        cmd += ["--region", args.region]
    cmd += args.extra
    print(f"  replicate {replicate}: {' '.join(cmd)}", flush=True)
    proc = subprocess.run(cmd, cwd=outdir, capture_output=True, text=True)
    if proc.returncode != 0:
        sys.stderr.write(proc.stdout[-4000:] + proc.stderr[-4000:])
        raise SystemExit(f"replicate {replicate} failed with exit {proc.returncode}")
    return outdir


def collect(outdir, prefix):
    found = {}
    for name in os.listdir(outdir):
        if not name.startswith(prefix):
            continue
        if any(name.endswith(s) for s in PRINCIPAL_SUFFIXES):
            found[name] = os.path.join(outdir, name)
    return found


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--lraa", default=os.path.join(os.path.dirname(__file__), "..", "..", "LRAA"))
    ap.add_argument("--genome", required=True)
    ap.add_argument("--bam", required=True)
    ap.add_argument("--gtf")
    ap.add_argument("--region", help="e.g. chr20:36000000-37000000")
    ap.add_argument("--replicates", type=int, default=3)
    ap.add_argument("--threads", type=int, default=4)
    ap.add_argument("--output-prefix", dest="output_prefix", default="det")
    ap.add_argument("--workdir", default="det_work")
    ap.add_argument("--keep", action="store_true")
    ap.add_argument("extra", nargs=argparse.REMAINDER,
                    help="further LRAA arguments after --")
    args = ap.parse_args()
    args.extra = [a for a in args.extra if a != "--"]

    if args.replicates < 3:
        raise SystemExit("use at least 3 replicates; two runs agreeing is weak evidence")

    workroot = os.path.abspath(args.workdir)
    shutil.rmtree(workroot, ignore_errors=True)
    os.makedirs(workroot)

    print(f"running {args.replicates} independent replicates under {workroot}")
    outdirs = [run_once(args, i + 1, workroot) for i in range(args.replicates)]

    sets = [collect(d, args.output_prefix) for d in outdirs]
    names = set(sets[0])
    for i, s in enumerate(sets[1:], start=2):
        if set(s) != names:
            raise SystemExit(
                f"replicate {i} produced a different set of outputs: "
                f"{sorted(set(s) ^ names)}"
            )
    if not names:
        raise SystemExit("no principal outputs found; check --output-prefix")

    content_bad, bytes_bad = [], []
    for name in sorted(names):
        cdig = {digest_content(s[name]) for s in sets}
        bdig = {digest_bytes(s[name]) for s in sets}
        status = "ok"
        if len(cdig) != 1:
            content_bad.append(name)
            status = "CONTENT DIFFERS"
        elif len(bdig) != 1:
            bytes_bad.append(name)
            status = "bytes differ, content identical"
        print(f"  {name:60} {status}")

    print()
    if content_bad:
        print("FAIL: content differs between replicates -- order-dependent selection")
        for n in content_bad:
            print(f"  {n}")
    if bytes_bad:
        print("FAIL: byte order differs while content matches -- merge ordering")
        for n in bytes_bad:
            print(f"  {n}")

    if not args.keep and not (content_bad or bytes_bad):
        shutil.rmtree(workroot, ignore_errors=True)

    if content_bad or bytes_bad:
        raise SystemExit(1)
    print(f"PASS: {len(names)} principal outputs identical in content and bytes "
          f"across {args.replicates} replicates")


if __name__ == "__main__":
    main()
