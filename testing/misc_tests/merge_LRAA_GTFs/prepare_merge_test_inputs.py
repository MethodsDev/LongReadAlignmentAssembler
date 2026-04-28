#!/usr/bin/env python3

import argparse
import glob
import os
import re
from pathlib import Path

ANNOTATE_CHOICES = {
    "pbmc_sc_ref-guided.4.LRAA.ref-guided.gtf": "t:chr19:+:comp-62:iso-4",
    "pbmc_sc_ref-guided.11.LRAA.ref-guided.gtf": "t:chr19:+:comp-95:iso-1",
    "pbmc_sc_ref-guided.10.LRAA.ref-guided.gtf": "t:chr19:+:comp-28:iso-1",
    "pbmc_sc_ref-guided.12.LRAA.ref-guided.gtf": "t:chr19:+:comp-99:iso-2",
}


def _info_has_flag(info: str, key: str) -> bool:
    m = re.search(rf'{key} "(True|False)"', info)
    return bool(m and m.group(1) == "True")


def _ensure_attr(info: str, key: str, value: str) -> str:
    if re.search(rf'{key} "[^"]*"', info):
        return re.sub(rf'{key} "[^"]*"', f'{key} "{value}"', info)
    info = info.rstrip()
    if not info.endswith(";"):
        info += ";"
    return info + f' {key} "{value}";'


def _simulate_count(tid: str, coord: int, base: int) -> int:
    # Stable pseudo-support in a realistic range.
    return base + (abs(hash((tid, coord))) % 23)


def main():
    ap = argparse.ArgumentParser(description="Prepare merge test GTF inputs with boundary annotations/counts")
    ap.add_argument("--src", default="data", help="source dir containing input gtfs and genome fasta")
    ap.add_argument("--out", default="prepared_inputs", help="output dir for prepared inputs")
    ap.add_argument(
        "--num_genes",
        type=int,
        default=12,
        help="Number of transcripts to inject with mixed boundary configurations",
    )
    args = ap.parse_args()

    src = Path(args.src)
    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)

    gtf_files = sorted(glob.glob(str(src / "*.gtf")))
    if not gtf_files:
        raise SystemExit(f"No GTF files found in {src}")

    file_lines = {}
    transcript_records = []

    for gtf in gtf_files:
        bn = os.path.basename(gtf)
        out_lines = []
        for line_idx, line in enumerate(open(gtf)):
            if line.startswith("#") or not line.strip():
                out_lines.append(line)
                continue

            f = line.rstrip("\n").split("\t")
            if len(f) < 9 or f[2] != "transcript":
                out_lines.append(line)
                continue

            info = f[8]
            m_tid = re.search(r'transcript_id "([^"]+)"', info)
            tid = m_tid.group(1) if m_tid else None
            if tid is not None:
                transcript_records.append(
                    {
                        "file": bn,
                        "line_idx": line_idx,
                        "tid": tid,
                        "lend": int(f[3]),
                        "rend": int(f[4]),
                        "strand": f[6],
                        "is_multi_exon": False,
                    }
                )
            out_lines.append(line)
        file_lines[bn] = out_lines

    # Mark multi-exon transcripts when possible.
    for gtf in gtf_files:
        bn = os.path.basename(gtf)
        exon_counts = {}
        for line in open(gtf):
            if line.startswith("#") or not line.strip():
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 9 or f[2] != "exon":
                continue
            m_tid = re.search(r'transcript_id "([^"]+)"', f[8])
            if m_tid:
                tid = m_tid.group(1)
                exon_counts[tid] = exon_counts.get(tid, 0) + 1
        for rec in transcript_records:
            if rec["file"] == bn and exon_counts.get(rec["tid"], 0) > 1:
                rec["is_multi_exon"] = True

    # Auto-pick fallback transcripts per file if hardcoded IDs are absent.
    forced = set((k, v) for k, v in ANNOTATE_CHOICES.items())
    picked = []
    picked_set = set()
    for bn, tid in forced:
        if any(r["file"] == bn and r["tid"] == tid for r in transcript_records):
            picked.append((bn, tid))
            picked_set.add((bn, tid))

    if not picked:
        seen_files = set()
        for rec in transcript_records:
            if rec["file"] in seen_files:
                continue
            if not rec["is_multi_exon"]:
                continue
            picked.append((rec["file"], rec["tid"]))
            picked_set.add((rec["file"], rec["tid"]))
            seen_files.add(rec["file"])
            if len(picked) >= min(4, args.num_genes):
                break
        if len(picked) < min(4, args.num_genes):
            for rec in transcript_records:
                if rec["file"] in seen_files:
                    continue
                picked.append((rec["file"], rec["tid"]))
                picked_set.add((rec["file"], rec["tid"]))
                seen_files.add(rec["file"])
                if len(picked) >= min(4, args.num_genes):
                    break

    # Expand to a handful of transcripts with varied boundary configurations.
    for rec in transcript_records:
        key = (rec["file"], rec["tid"])
        if key in picked_set:
            continue
        if not rec["is_multi_exon"]:
            continue
        picked.append(key)
        picked_set.add(key)
        if len(picked) >= max(4, args.num_genes):
            break

    num_true_tss = 0
    num_true_polya = 0

    for rec in transcript_records:
        bn = rec["file"]
        line_idx = rec["line_idx"]
        line = file_lines[bn][line_idx]
        f = line.rstrip("\n").split("\t")
        info = f[8]

        # Ensure explicit boolean attrs exist.
        info = _ensure_attr(info, "TSS", "True" if _info_has_flag(info, "TSS") else "False")
        info = _ensure_attr(info, "PolyA", "True" if _info_has_flag(info, "PolyA") else "False")

        # Inject mixed boundary configurations across selected transcripts:
        # idx % 4:
        #   0 => TSS=True,  PolyA=True
        #   1 => TSS=True,  PolyA=False
        #   2 => TSS=False, PolyA=True
        #   3 => TSS=False, PolyA=False
        if (bn, rec["tid"]) in picked_set:
            idx = picked.index((bn, rec["tid"]))
            mode = idx % 4
            want_tss = mode in (0, 1)
            want_polya = mode in (0, 2)
            info = _ensure_attr(info, "TSS", "True" if want_tss else "False")
            info = _ensure_attr(info, "PolyA", "True" if want_polya else "False")

        tss_true = _info_has_flag(info, "TSS")
        polya_true = _info_has_flag(info, "PolyA")

        if tss_true:
            num_true_tss += 1
            tss_coord = rec["lend"] if rec["strand"] == "+" else rec["rend"]
            if not re.search(r'TSS_read_count "[0-9]+"', info):
                val = _simulate_count(rec["tid"], tss_coord, 8)
                info = _ensure_attr(info, "TSS_read_count", str(val))

        if polya_true:
            num_true_polya += 1
            polya_coord = rec["rend"] if rec["strand"] == "+" else rec["lend"]
            if not re.search(r'PolyA_read_count "[0-9]+"', info):
                val = _simulate_count(rec["tid"], polya_coord, 10)
                info = _ensure_attr(info, "PolyA_read_count", str(val))
        else:
            info = re.sub(r'\s*PolyA_read_count "[0-9]+";', "", info)

        if not tss_true:
            info = re.sub(r'\s*TSS_read_count "[0-9]+";', "", info)

        f[8] = info
        file_lines[bn][line_idx] = "\t".join(f) + "\n"

    for bn, out_lines in file_lines.items():
        with open(out / bn, "wt") as ofh:
            ofh.writelines(out_lines)

    # Copy genome fasta/index if present.
    for aux in ("chr19.genome.fa", "chr19.genome.fa.fai"):
        p = src / aux
        if p.exists():
            (out / aux).write_bytes(p.read_bytes())

    if num_true_tss == 0 or num_true_polya == 0:
        raise SystemExit(
            f"Prepared inputs invalid: TSS_true={num_true_tss}, PolyA_true={num_true_polya}."
        )

    print(f"Prepared {len(gtf_files)} GTF files in {out}")
    print(f"Boundary-annotated transcripts: TSS_true={num_true_tss}, PolyA_true={num_true_polya}")


if __name__ == "__main__":
    main()
