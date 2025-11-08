#!/usr/bin/env python3
"""
GeneCommunityCluster

Leiden-based community clustering of transcripts inside an initial overlap-derived
isoform component to refine gene groupings.

Design goals:
- Lightweight: binary edges only between sufficiently overlapping transcripts.
- Deterministic: seed-controlled Leiden run.
- Optional: can be disabled via config flag (use_community_clustering).

Public function:
    partition_with_leiden(transcripts, contig_acc, contig_strand, overlap_params, resolution, seed)
        Returns list[int] community_ids aligned to input transcript order.

If igraph or leidenalg are unavailable, raises ImportError (feature is gated by a config flag).

Edge inclusion criteria mimic existing reclustering thresholds:
    (overlap_len / shorter_len) >= min_recluster_overlap_shorter_iso_frac
AND (overlap_len / longer_len) >= min_recluster_overlap_longer_iso_frac
Containment (shorter fully overlapped) also forces an edge.

Future extensions (not implemented yet):
    - Singleton reassignment heuristics.
    - Edge sparsification for very large components.
"""
from typing import List
import LRAA_Globals

try:
    import igraph  # type: ignore
    import leidenalg  # type: ignore
except ImportError as e:
    # Defer error until function call so the rest of LRAA can run without clustering.
    IG_AVAILABLE = False
    IMPORT_ERROR = e
else:
    IG_AVAILABLE = True
    IMPORT_ERROR = None


def _merge_intervals(exons):
    merged = []
    for a, b in sorted(exons):
        if not merged or a > merged[-1][1] + 1:
            merged.append([a, b])
        else:
            if b > merged[-1][1]:
                merged[-1][1] = b
    return merged


def _transcript_overlap_len(t1, t2) -> int:
    m1 = _merge_intervals(t1.get_exon_segments())
    m2 = _merge_intervals(t2.get_exon_segments())
    i = j = 0
    ov = 0
    while i < len(m1) and j < len(m2):
        a1, b1 = m1[i]
        a2, b2 = m2[j]
        if b1 < a2:
            i += 1
        elif b2 < a1:
            j += 1
        else:
            ov += min(b1, b2) - max(a1, a2) + 1
            if b1 < b2:
                i += 1
            elif b2 < b1:
                j += 1
            else:
                i += 1
                j += 1
    return ov


def partition_with_leiden(transcripts: List, contig_acc: str, contig_strand: str,
                           resolution: float = None, seed: int = None) -> List[int]:
    """Run Leiden community detection on a binary overlap graph.

    Parameters
    ----------
    transcripts : List[Transcript]
        Transcripts belonging to a single initial overlap component.
    contig_acc : str
    contig_strand : str
    resolution : float
        Leiden resolution parameter; uses config if None.
    seed : int
        Random seed for reproducibility; uses config if None.

    Returns
    -------
    List[int]
        Community membership per transcript index (0..k-1).
    """
    if not LRAA_Globals.config.get("use_community_clustering", False):
        raise RuntimeError("partition_with_leiden called but use_community_clustering is False")

    if not IG_AVAILABLE:
        raise ImportError(f"Leiden clustering requested but igraph/leidenalg not available: {IMPORT_ERROR}")

    if resolution is None:
        resolution = float(LRAA_Globals.config.get("community_resolution", 1.0))
    if seed is None:
        seed = int(LRAA_Globals.config.get("community_random_seed", 42))

    # thresholds from config
    shorter_thr = float(LRAA_Globals.config.get("min_recluster_overlap_shorter_iso_frac", 0.50))
    longer_thr = float(LRAA_Globals.config.get("min_recluster_overlap_longer_iso_frac", 0.20))

    n = len(transcripts)
    if n == 0:
        return []
    if n == 1:
        return [0]

    # Build edges
    edges = []
    for i in range(n):
        ti = transcripts[i]
        len_i = ti.get_cdna_len()
        for j in range(i + 1, n):
            tj = transcripts[j]
            ov = _transcript_overlap_len(ti, tj)
            if ov <= 0:
                continue
            len_j = tj.get_cdna_len()
            shorter = min(len_i, len_j)
            longer = max(len_i, len_j)
            # containment or overlap criteria
            contained = (ov == shorter and shorter > 0)
            if contained or (
                shorter > 0
                and longer > 0
                and (ov / shorter) >= shorter_thr
                and (ov / longer) >= longer_thr
            ):
                edges.append((i, j))

    # If no edges, each transcript becomes its own community
    if not edges:
        return list(range(n))

    g = igraph.Graph(n=n, edges=edges, directed=False)

    try:
        part = leidenalg.find_partition(
            g,
            leidenalg.RBConfigurationVertexPartition,
            resolution_parameter=resolution,
            seed=seed,
        )
    except Exception as e:
        raise RuntimeError(f"Leiden clustering failed for {contig_acc}:{contig_strand}: {e}")

    membership = list(part.membership)
    # Normalize community ids to 0..k-1 in order of first occurrence for consistency
    remap = {}
    next_id = 0
    norm_membership = []
    for cid in membership:
        if cid not in remap:
            remap[cid] = next_id
            next_id += 1
        norm_membership.append(remap[cid])

    return norm_membership

__all__ = ["partition_with_leiden"]
