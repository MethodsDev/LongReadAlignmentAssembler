#!/usr/bin/env python3
# encoding: utf-8

"""Whole-genome alignment-mismapping filter (LRAA v0.40.0).

Removes isoforms that are alignment/strand-mismapping artifacts of a
much-higher-expressed transcript. This is a POST-MERGE, whole-genome stage: it
consumes the merged genome-wide gtf + quant.expr (never a per-chunk slice) and
rewrites both, plus a log naming every removed model and why.

Two detectors, unioned (see LRAA_Globals config block for the rationale):

  MIRROR (coordinate): a multi-exon model on the strand OPPOSITE a
    higher-expressed model, with high exonic base overlap and every internal
    splice site within `mismap_junction_tolerance` bp of that model's exon
    boundaries. Catches wrong-strand near-mirror artifacts.

  SEQUENCE (minimap2 cDNA all-vs-all): a model whose spliced cDNA is
    >= `mismap_min_seq_identity`% identical over >= `mismap_min_seq_coverage`
    of its length to a DIFFERENT-gene higher-expressed model. Catches run-ons,
    chimeras and same-strand mismappings.

Both gate on `mismap_max_expr_fraction`: a model is removed only when it carries
less than that fraction of the matched model's expression. That low-expression
gate is what makes the quant repairable by dropping the removed rows and
renormalizing TPM to 1e6 rather than requantifying (the removed reads are
artifacts owed to no surviving model).

Expression is read from the quant.expr `all_reads` column and joined to the GTF
on transcript_id ONLY (the quant `gene_id` column may carry a gene_symbol^
prefix); the same-gene test uses the GTF's own gene_id attribute.
"""

import os
import re
import csv
import logging
import subprocess
from collections import defaultdict

import pysam

import LRAA_Globals
from Transcript import GTF_contig_to_transcripts

logger = logging.getLogger(__name__)


# ---------------------------------------------------------------------------
# small helpers
# ---------------------------------------------------------------------------
def _reverse_complement(seq):
    return seq.translate(str.maketrans("ACGTNacgtn", "TGCANtgcan"))[::-1]


def _build_cdna(transcript, contig_seq_str):
    # mirrors IsoformReadRescue._build_transcript_sequence
    exon_segments = transcript.get_exon_segments()
    exon_seqs = [contig_seq_str[lend - 1 : rend] for lend, rend in exon_segments]
    if transcript.get_strand() == "-":
        exon_seqs = [_reverse_complement(s) for s in reversed(exon_seqs)]
    return "".join(exon_seqs)


def _internal_splice_sites(exon_segments):
    """Genomic coords of the internal exon boundaries (donor/acceptor sites)."""
    ex = sorted(exon_segments)
    sites = set()
    for i, (lend, rend) in enumerate(ex):
        if i > 0:
            sites.add(lend)
        if i < len(ex) - 1:
            sites.add(rend)
    return sites


def _exon_boundaries(exon_segments):
    b = set()
    for lend, rend in exon_segments:
        b.add(lend)
        b.add(rend)
    return b


def _exon_overlap_bp(a, b):
    """Overlapping exonic bp between two sorted [(lend,rend)] exon lists."""
    a = sorted(a)
    b = sorted(b)
    i = j = tot = 0
    while i < len(a) and j < len(b):
        s = max(a[i][0], b[j][0])
        e = min(a[i][1], b[j][1])
        if s <= e:
            tot += e - s + 1
        if a[i][1] < b[j][1]:
            i += 1
        else:
            j += 1
    return tot


def _exonic_length(exon_segments):
    return sum(rend - lend + 1 for lend, rend in exon_segments)


def _merged_query_cov(intervals, qlen):
    if qlen <= 0:
        return 0.0
    m = []
    for a, b in sorted((min(x, y), max(x, y)) for x, y in intervals):
        if m and a <= m[-1][1] + 1:
            m[-1][1] = max(m[-1][1], b)
        else:
            m.append([a, b])
    return sum(b - a + 1 for a, b in m) / qlen


# ---------------------------------------------------------------------------
# input loading
# ---------------------------------------------------------------------------
def _leading_comment_lines(path):
    out = []
    with open(path, "rt") as fh:
        for line in fh:
            if line.startswith("#"):
                out.append(line.rstrip("\n"))
            else:
                break
    return out


def _load_transcripts(gtf_filename):
    """transcript_id -> dict(obj, gene_id, strand, contig, exons, elen)."""
    contig_to_transcripts = GTF_contig_to_transcripts.parse_GTF_to_Transcripts(
        gtf_filename
    )
    info = {}
    for contig, transcripts in contig_to_transcripts.items():
        for t in transcripts:
            exons = sorted(t.get_exon_segments())
            info[t.get_transcript_id()] = dict(
                obj=t,
                gene_id=t.get_gene_id(),
                strand=t.get_strand(),
                contig=t.get_contig_acc(),
                exons=exons,
                elen=_exonic_length(exons),
            )
    return info


def _load_quant(quant_filename):
    """Returns (comments, fieldnames, rows_in_order, tid->all_reads)."""
    comments = _leading_comment_lines(quant_filename)

    def _non_comment(fh):
        for line in fh:
            if not line.startswith("#"):
                yield line

    rows = []
    expr = {}
    with open(quant_filename, "rt", newline="") as fh:
        reader = csv.DictReader(_non_comment(fh), delimiter="\t")
        fieldnames = reader.fieldnames
        for req in ("transcript_id", "all_reads", "TPM"):
            if req not in fieldnames:
                raise RuntimeError(
                    f"quant.expr {quant_filename} missing required column '{req}'"
                )
        for row in reader:
            rows.append(row)
            try:
                expr[row["transcript_id"]] = float(row["all_reads"])
            except (TypeError, ValueError):
                expr[row["transcript_id"]] = 0.0
    return comments, fieldnames, rows, expr


# ---------------------------------------------------------------------------
# detectors
# ---------------------------------------------------------------------------
def _detect_mirror(info, expr, tol, min_base_overlap, max_expr_fraction):
    """Return {tid: (partner_tid, base_overlap, expr_ratio)} flagged by the
    coordinate mirror rule."""
    flagged = {}
    by_contig = defaultdict(list)
    for tid, d in info.items():
        by_contig[d["contig"]].append(tid)

    for contig, tids in by_contig.items():
        recs = [(t, info[t]) for t in tids]
        for tid, d in recs:
            if len(d["exons"]) < 2:
                continue  # mirror needs splice sites
            m_lend, m_rend = d["exons"][0][0], d["exons"][-1][1]
            m_sites = _internal_splice_sites(d["exons"])
            best = None  # (partner_expr, partner_tid, base_overlap)
            for ntid, nd in recs:
                if nd["strand"] == d["strand"]:
                    continue
                # span overlap first (cheap)
                if not (m_lend <= nd["exons"][-1][1] and nd["exons"][0][0] <= m_rend):
                    continue
                ov_bp = _exon_overlap_bp(d["exons"], nd["exons"])
                if ov_bp <= 0:
                    continue
                base_ov = ov_bp / d["elen"]
                if base_ov < min_base_overlap:
                    continue
                nb = _exon_boundaries(nd["exons"])
                # every internal splice site within tol of an opp-strand boundary
                if not all(
                    any(abs(s - b) <= tol for b in nb) for s in m_sites
                ):
                    continue
                ne = expr.get(ntid, 0.0)
                if best is None or ne > best[0]:
                    best = (ne, ntid, base_ov)
            if best is None:
                continue
            partner_expr, partner_tid, base_ov = best
            me = expr.get(tid, 0.0)
            if partner_expr > 0 and me < max_expr_fraction * partner_expr:
                flagged[tid] = (partner_tid, base_ov, me / partner_expr)
    return flagged


def _run_minimap2_ava(cdna_fa, threads):
    """All-vs-all cDNA alignment; yields (qname, tname, identity, qlen, qstart, qend)."""
    cmd = [
        "minimap2",
        "-c",
        "-x",
        "map-ont",
        "-N",
        "20",
        "-p",
        "0.1",
        "-t",
        str(threads),
        cdna_fa,
        cdna_fa,
    ]
    logger.info("Running: %s", " ".join(cmd))
    proc = subprocess.Popen(
        cmd, stdout=subprocess.PIPE, stderr=subprocess.DEVNULL, text=True
    )
    for line in proc.stdout:
        f = line.rstrip("\n").split("\t")
        if len(f) < 12:
            continue
        qname, tname = f[0], f[5]
        qlen, qstart, qend = int(f[1]), int(f[2]), int(f[3])
        de = None
        for tag in f[12:]:
            if tag.startswith("de:f:"):
                de = float(tag[5:])
                break
        identity = (1.0 - de) * 100.0 if de is not None else 100.0
        yield qname, tname, identity, qlen, qstart, qend
    rc = proc.wait()
    if rc != 0:
        raise RuntimeError(f"minimap2 all-vs-all failed with exit code {rc}")


def _write_all_cdna(info, genome_fasta, cdna_fa):
    """Write spliced cDNA for every model; contigs fetched once each."""
    n = 0
    with pysam.FastaFile(genome_fasta) as fasta, open(cdna_fa, "wt") as ofh:
        # group by contig so each contig sequence is fetched exactly once
        by_contig = defaultdict(list)
        for tid, d in info.items():
            by_contig[d["contig"]].append(tid)
        for contig in by_contig:
            try:
                contig_seq = fasta.fetch(contig).upper()
            except (KeyError, ValueError):
                logger.warning("contig %s absent from genome fasta; skipping", contig)
                continue
            for tid in by_contig[contig]:
                seq = _build_cdna(info[tid]["obj"], contig_seq)
                if seq:
                    ofh.write(f">{tid}\n{seq}\n")
                    n += 1
    return n


def _detect_sequence(info, expr, genome_fasta, min_identity, min_cov,
                     max_expr_fraction, workdir, threads):
    """Return {tid: (partner_tid, identity, qcov, expr_ratio)} flagged by the
    cDNA sequence rule (cross-gene, near-identical, much-higher-expressed)."""
    cdna_fa = os.path.join(workdir, "lraa_cdna.fa")
    n = _write_all_cdna(info, genome_fasta, cdna_fa)
    logger.info("Wrote %d cDNA sequences for all-vs-all alignment", n)

    # single minimap2 pass: per cross-gene (q,t) pair, keep query intervals of
    # hits >= min_identity and the best identity seen for that pair.
    pair_intervals = defaultdict(list)
    pair_identity = {}
    pair_qlen = {}
    for q, t, identity, qlen, qs, qe in _run_minimap2_ava(cdna_fa, threads):
        if q == t:
            continue
        if info.get(q, {}).get("gene_id") == info.get(t, {}).get("gene_id"):
            continue  # same-gene excluded
        if identity < min_identity:
            continue
        key = (q, t)
        pair_intervals[key].append((qs, qe))
        pair_qlen[q] = qlen
        if key not in pair_identity or identity > pair_identity[key]:
            pair_identity[key] = identity

    # per query: the highest-expressed cross-gene target passing coverage
    maxhit = {}
    for (q, t), intervals in pair_intervals.items():
        if _merged_query_cov(intervals, pair_qlen.get(q, 0)) < min_cov:
            continue
        te = expr.get(t, 0.0)
        if q not in maxhit or te > maxhit[q][0]:
            qcov = _merged_query_cov(intervals, pair_qlen.get(q, 0))
            maxhit[q] = (te, t, pair_identity[(q, t)], qcov)

    flagged = {}
    for q, (te, t, identity, qcov) in maxhit.items():
        me = expr.get(q, 0.0)
        if te > 0 and me < max_expr_fraction * te:
            flagged[q] = (t, identity, qcov, me / te)
    return flagged


# ---------------------------------------------------------------------------
# output writing
# ---------------------------------------------------------------------------
_TRANSCRIPT_ID_RE = re.compile(r'transcript_id "([^"]+)"')


def _write_filtered_gtf(gtf_in, gtf_out, drop_set):
    """Stream the input GTF through unchanged, emitting every line except those
    of transcripts in drop_set.

    Byte-preserving for every surviving line: no round-trip through Transcript
    objects. That keeps exon coordinates and the attribute string identical to
    the input (a reserialize both duplicated attributes and inserted blank
    lines) and makes it structurally impossible for the filter to alter a model
    it did not remove -- only the dropped models disappear. Any line without a
    transcript_id (e.g. leading comments, blank lines) is passed through as-is.
    """
    with open(gtf_in, "rt") as ifh, open(gtf_out, "wt") as ofh:
        for line in ifh:
            m = _TRANSCRIPT_ID_RE.search(line)
            if m is not None and m.group(1) in drop_set:
                continue
            ofh.write(line)


def _write_filtered_quant(comments, fieldnames, rows, drop_set, info, quant_out):
    survivors = [r for r in rows if r["transcript_id"] not in drop_set]

    # renormalize TPM to 1e6 over survivors
    tpm_sum = 0.0
    for r in survivors:
        try:
            tpm_sum += float(r["TPM"])
        except (TypeError, ValueError):
            pass
    scale = (1e6 / tpm_sum) if tpm_sum > 0 else 1.0
    for r in survivors:
        try:
            r["TPM"] = f"{float(r['TPM']) * scale:.3f}"
        except (TypeError, ValueError):
            pass

    # recompute isoform_fraction / unique_gene_read_fraction ONLY for genes that
    # lost a member (unchanged genes keep their exact original values). Group on
    # the GTF bare gene_id (join by transcript_id), matching Quantify.
    affected_genes = set()
    for tid in drop_set:
        g = info.get(tid, {}).get("gene_id")
        if g is not None:
            affected_genes.add(g)
    has_iso = "isoform_fraction" in fieldnames
    has_ugf = "unique_gene_read_fraction" in fieldnames
    if affected_genes and (has_iso or has_ugf):
        gene_rows = defaultdict(list)
        for r in survivors:
            g = info.get(r["transcript_id"], {}).get("gene_id")
            if g in affected_genes:
                gene_rows[g].append(r)
        for g, grows in gene_rows.items():
            gene_all = 0.0
            for r in grows:
                try:
                    gene_all += float(r["all_reads"])
                except (TypeError, ValueError):
                    pass
            for r in grows:
                if gene_all > 0:
                    if has_iso:
                        r["isoform_fraction"] = f"{float(r['all_reads']) / gene_all:.3f}"
                    if has_ugf:
                        r["unique_gene_read_fraction"] = (
                            f"{float(r['uniq_reads']) / gene_all:.3f}"
                        )

    with open(quant_out, "wt", newline="") as ofh:
        for c in comments:
            ofh.write(c + "\n")
        writer = csv.DictWriter(
            ofh, fieldnames=fieldnames, delimiter="\t", lineterminator="\n"
        )
        writer.writeheader()
        for r in survivors:
            writer.writerow(r)


def _write_log(log_out, mirror, sequence, info, expr, total_models, n_survivors):
    with open(log_out, "wt") as ofh:
        ofh.write(
            "# LRAA alignment-mismapping filter: {} of {} models removed, {} retained\n".format(
                len(set(mirror) | set(sequence)), total_models, n_survivors
            )
        )
        ofh.write(
            "\t".join(
                [
                    "transcript_id",
                    "gene_id",
                    "detector",
                    "partner_transcript_id",
                    "detail",
                    "model_reads",
                    "partner_reads",
                    "expr_fraction_of_partner",
                ]
            )
            + "\n"
        )
        for tid in sorted(set(mirror) | set(sequence)):
            detectors = []
            partner = None
            detail = []
            ratio = None
            if tid in mirror:
                p, base_ov, r = mirror[tid]
                detectors.append("mirror")
                partner = p
                detail.append(f"base_overlap={base_ov:.2f}")
                ratio = r
            if tid in sequence:
                p, identity, qcov, r = sequence[tid]
                detectors.append("sequence")
                partner = partner or p
                detail.append(f"identity={identity:.1f}%,qcov={qcov:.2f}")
                if ratio is None:
                    ratio = r
                    partner = p
            ofh.write(
                "\t".join(
                    [
                        tid,
                        str(info.get(tid, {}).get("gene_id", "")),
                        "+".join(detectors),
                        str(partner),
                        ";".join(detail),
                        f"{expr.get(tid, 0.0):.1f}",
                        f"{expr.get(partner, 0.0):.1f}",
                        f"{ratio:.5f}" if ratio is not None else "",
                    ]
                )
                + "\n"
            )


# ---------------------------------------------------------------------------
# entry point
# ---------------------------------------------------------------------------
def run_mismapping_filter(
    gtf_in, quant_in, genome_fasta, gtf_out, quant_out, log_out, workdir, threads=1,
    exempt_contigs=None,
):
    """Filter a merged whole-genome GTF + quant.expr. Returns the drop-set.

    exempt_contigs: contigs whose models are carried forward untouched. Oversimplify
    contigs (e.g. chrM) hold reference models copied forward regardless of the reads,
    so a coverage-driven filter must not touch them: doing so removes DIFFERENT models
    in different cluster-guided inputs, and merge_LRAA_GTFs then refuses to merge an
    oversimplified contig whose per-input record sets disagree. Same reasoning as the
    read floor _run_oversimplify_best_overlap already declines to apply.
    """
    cfg = LRAA_Globals.config
    tol = int(cfg["mismap_junction_tolerance"])
    min_base_overlap = float(cfg["mismap_min_base_overlap"])
    min_identity = float(cfg["mismap_min_seq_identity"])
    min_cov = float(cfg["mismap_min_seq_coverage"])
    max_expr_fraction = float(cfg["mismap_max_expr_fraction"])

    logger.info("Alignment-mismapping filter: loading %s", gtf_in)
    info = _load_transcripts(gtf_in)
    comments, fieldnames, rows, expr = _load_quant(quant_in)

    logger.info(
        "Loaded %d transcripts; running mirror + sequence detectors", len(info)
    )
    mirror = _detect_mirror(info, expr, tol, min_base_overlap, max_expr_fraction)
    logger.info("Mirror detector flagged %d models", len(mirror))
    sequence = _detect_sequence(
        info, expr, genome_fasta, min_identity, min_cov,
        max_expr_fraction, workdir, threads,
    )
    logger.info("Sequence detector flagged %d models", len(sequence))

    drop_set = set(mirror) | set(sequence)

    if exempt_contigs:
        exempt_contigs = set(exempt_contigs)
        exempted = {
            tid for tid in drop_set
            if info.get(tid, {}).get("contig") in exempt_contigs
        }
        if exempted:
            logger.info(
                "Exempting %d flagged model(s) on carried-forward (oversimplify) "
                "contig(s) %s from removal",
                len(exempted), ",".join(sorted(exempt_contigs)),
            )
            drop_set -= exempted
            # keep the log's per-detector lists consistent with what was removed
            for tid in exempted:
                mirror.pop(tid, None)
                sequence.pop(tid, None)

    n_survivors = len(info) - len(drop_set)
    logger.info(
        "Removing %d of %d models (%d survive)", len(drop_set), len(info), n_survivors
    )

    _write_filtered_gtf(gtf_in, gtf_out, drop_set)
    _write_filtered_quant(comments, fieldnames, rows, drop_set, info, quant_out)
    _write_log(log_out, mirror, sequence, info, expr, len(info), n_survivors)
    return drop_set
