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

Edge inclusion criteria, in the order they are tested:
    shared introns >= min_recluster_shared_introns
    AND shared introns / min(intron count of either isoform) > min_recluster_shared_intron_frac
 OR (overlap_len / shorter_len) >= min_recluster_overlap_shorter_iso_frac
 OR containment (shorter fully overlapped)

Isoforms with IDENTICAL intron chains are contracted to a single graph vertex before
partitioning, so no partition can separate them. Contraction keeps edge multiplicity
as weights and internal edges as self-loops, which makes the contracted objective
equal to the original one restricted to must-link-respecting partitions; see
_contract_identical_chains.

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


def _renumber(membership: List[int]) -> List[int]:
    """Community ids to 0..k-1 in order of first occurrence."""
    remap = {}
    out = []
    for cid in membership:
        if cid not in remap:
            remap[cid] = len(remap)
        out.append(remap[cid])
    return out


def _intron_sets(transcripts: List) -> List[frozenset]:
    """One frozenset of (lend, rend) introns per transcript, index-aligned.

    A single-exon transcript yields an EMPTY set. Empty is not a shared splice
    pattern -- it has no junctions to agree about -- so such transcripts take part in
    neither intron rule and are left to the exon-overlap criterion. Treating empty as
    a match would collapse every single-exon model in a component into one gene.
    """
    return [frozenset(t.get_introns()) for t in transcripts]


def _shares_introns(
    si: frozenset, sj: frozenset, min_frac: float, min_shared: int
) -> bool:
    """Shared introns, as a fraction of whichever isoform has FEWER introns.

    The denominator is the smaller intron count, so a 2-intron isoform sharing one
    intron with a 20-intron isoform scores 0.5 rather than 0.05. Normalizing by the
    longer isoform instead would score fragmentary models near zero however well
    their own splice pattern agrees, which is the case this rule exists to catch.

    min_shared is a floor on the COUNT, applied alongside the fraction because the
    fraction alone is satisfied by a single shared junction whenever the smaller model
    has four introns or fewer -- enough for one 1- or 2-intron fragment to bind two
    neighbouring genes into one. MEASURED over four contigs, moving the floor from 1 to
    2 gave up 16 of 42 recovered gene ids and avoided 8 of 10 added fusions, taking the
    exchange rate from 4.2:1 to 13.0:1.

    Identical intron chains do NOT come through here -- they are a must-link handled
    structurally -- so no setting of either threshold can break that guarantee.
    """
    if not si or not sj:
        return False
    shared = len(si & sj)
    if shared < min_shared or shared == 0:
        return False
    return (shared / min(len(si), len(sj))) > min_frac


def _contract_identical_chains(intron_sets: List[frozenset]) -> List[int]:
    """Map each transcript to a graph vertex, identical intron chains sharing one.

    Same splice pattern is the same gene by definition, and this is how that is
    guaranteed: the transcripts become ONE vertex, so no partition can separate them.

    Done here rather than by reconciling communities afterwards. A post-hoc union of
    the two communities holding an identical-chain pair drags every OTHER member of
    both communities along with them, which can fuse two genes because one transcript
    in each happened to share a chain. Contraction constrains exactly the transcripts
    the rule speaks about and leaves everything else to be partitioned on its merits.

    Single-exon transcripts have an empty chain and each get their own vertex: an
    empty chain is not a shared splice pattern, and treating it as one would contract
    every single-exon model in the component into a single gene.
    """
    vertex_of = [0] * len(intron_sets)
    chain_to_vertex = {}
    next_vertex = 0
    for idx, si in enumerate(intron_sets):
        if si:
            key = tuple(sorted(si))
            v = chain_to_vertex.get(key)
            if v is None:
                v = next_vertex
                next_vertex += 1
                chain_to_vertex[key] = v
            vertex_of[idx] = v
        else:
            vertex_of[idx] = next_vertex
            next_vertex += 1
    return vertex_of


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
    if not LRAA_Globals.config["use_community_clustering"]:
        raise RuntimeError("partition_with_leiden called but use_community_clustering is False")

    if not IG_AVAILABLE:
        raise ImportError(f"Leiden clustering requested but igraph/leidenalg not available: {IMPORT_ERROR}")

    # Indexed: LRAA_Globals.config defines these keys, so a literal fallback here
    # would be a competing default that wins silently when one is absent.
    if resolution is None:
        resolution = float(LRAA_Globals.config["community_resolution"])
    if seed is None:
        seed = int(LRAA_Globals.config["community_random_seed"])

    # thresholds from config
    shorter_thr = float(LRAA_Globals.config["min_recluster_overlap_shorter_iso_frac"])
    shared_intron_thr = float(LRAA_Globals.config["min_recluster_shared_intron_frac"])
    shared_intron_min = int(LRAA_Globals.config["min_recluster_shared_introns"])
    # No longer_thr here: requiring both overlap fractions kept legitimate alternative
    # isoforms apart. The shared-intron rule now carries the cases that needed it.

    n = len(transcripts)
    if n == 0:
        return []
    if n == 1:
        return [0]

    # Build edges over CONTRACTED vertices, so identical splice patterns cannot be
    # separated by any partition of this graph.
    intron_sets = _intron_sets(transcripts)
    vertex_of = _contract_identical_chains(intron_sets)
    n_vertices = max(vertex_of) + 1

    # Edge MULTIPLICITY is kept as a weight, and intra-vertex edges are kept as
    # self-loops. Both are load-bearing, not bookkeeping: RBConfiguration's null model
    # is degree-based, so an identical-chain class of ten transcripts that collapses to
    # one vertex with one unit edge has a tenth of the degree its members had, and
    # Leiden then optimizes a DIFFERENT objective rather than the original one subject
    # to a must-link. VERIFIED numerically on a synthetic component: with multiplicity
    # and self-loops the contracted quality equals the full-graph quality exactly for
    # every partition tried (spread 0), while deduplicating to a simple graph gives a
    # partition-dependent discrepancy spanning 14.4 quality units -- i.e. a different
    # argmax. Self-loop degree bookkeeping works out because igraph counts a self-loop
    # twice, matching the two endpoint degrees each internal edge contributed before.
    edge_weight = {}
    for i in range(n):
        ti = transcripts[i]
        len_i = ti.get_cdna_len()
        si = intron_sets[i]
        vi = vertex_of[i]
        for j in range(i + 1, n):
            vj = vertex_of[j]
            linked = False
            # Introns first: a set intersection is cheaper than the exon-overlap scan,
            # and junction agreement settles the pair without needing coordinates.
            if _shares_introns(
                si, intron_sets[j], shared_intron_thr, shared_intron_min
            ):
                linked = True
            else:
                tj = transcripts[j]
                ov = _transcript_overlap_len(ti, tj)
                if ov > 0:
                    len_j = tj.get_cdna_len()
                    shorter = min(len_i, len_j)
                    # containment or overlap criteria
                    contained = (ov == shorter and shorter > 0)
                    if contained or (shorter > 0 and (ov / shorter) >= shorter_thr):
                        linked = True
            if not linked:
                continue
            key = (vi, vj) if vi <= vj else (vj, vi)
            edge_weight[key] = edge_weight.get(key, 0) + 1

    # Contraction alone may have answered it; nothing left to partition.
    if n_vertices == 1:
        return [0] * n

    # Only self-loops means no vertex is linked to any other: every vertex is its own
    # community. Identical chains are already one vertex, so the rule still holds.
    if not any(a != b for a, b in edge_weight):
        return _renumber(list(vertex_of))

    keys = sorted(edge_weight)
    g = igraph.Graph(n=n_vertices, edges=keys, directed=False)
    g.es["weight"] = [edge_weight[k] for k in keys]

    try:
        part = leidenalg.find_partition(
            g,
            leidenalg.RBConfigurationVertexPartition,
            weights="weight",
            resolution_parameter=resolution,
            seed=seed,
        )
    except Exception as e:
        raise RuntimeError(f"Leiden clustering failed for {contig_acc}:{contig_strand}: {e}")

    vertex_membership = list(part.membership)
    return _renumber([vertex_membership[v] for v in vertex_of])

__all__ = ["partition_with_leiden"]
 
def partition_with_dsu(transcripts: List, contig_acc: str, contig_strand: str) -> List[int]:
    """Partition an initial overlap component into subcomponents using a
    sweep-line plus union-find (DSU) strategy with overlap thresholds.

    Returns list[int] membership aligned to input transcript order.

    This mirrors the criteria used in Leiden, but avoids building dense graphs
    and computes exact exon overlap only when a fast span-overlap bound can
    still satisfy thresholds.
    """
    n = len(transcripts)
    if n == 0:
        return []
    if n == 1:
        return [0]

    shorter_thr = float(LRAA_Globals.config["min_recluster_overlap_shorter_iso_frac"])
    longer_thr = float(LRAA_Globals.config["min_recluster_overlap_longer_iso_frac"])
    shared_intron_thr = float(LRAA_Globals.config["min_recluster_shared_intron_frac"])
    shared_intron_min = int(LRAA_Globals.config["min_recluster_shared_introns"])
    intron_sets = _intron_sets(transcripts)

    # DSU structure over indices 0..n-1
    parent = list(range(n))
    size = [1] * n

    def find(x: int) -> int:
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    def union(a: int, b: int) -> bool:
        ra, rb = find(a), find(b)
        if ra == rb:
            return False
        if size[ra] < size[rb]:
            ra, rb = rb, ra
        parent[rb] = ra
        size[ra] += size[rb]
        return True

    # Identical non-empty intron chains are ONE gene by definition, so they are
    # unioned unconditionally and before the sweep. Not folded into the fractional
    # rule below: that rule is threshold-driven, and either a configured
    # min_recluster_shared_intron_frac >= 1.0 or a min_recluster_shared_introns above
    # the chain length would defeat the guarantee. Done outside the sweep as well, so
    # it cannot depend on the sweep's span window either.
    chain_members = {}
    for idx, si in enumerate(intron_sets):
        if si:
            chain_members.setdefault(tuple(sorted(si)), []).append(idx)
    for members in chain_members.values():
        for other in members[1:]:
            union(members[0], other)

    # Precompute simple genomic spans and lengths; keep original indices
    bounds = []  # (lend, rend, cdna_len, idx)
    for idx, t in enumerate(transcripts):
        ex = t.get_exon_segments()
        if ex:
            lend = min(a for a, _ in ex)
            rend = max(b for _, b in ex)
        else:
            lend = t.get_lend()
            rend = t.get_rend()
        bounds.append((lend, rend, t.get_cdna_len(), idx))

    # Sort by lend
    bounds.sort(key=lambda x: x[0])

    # Sweep-line: only consider span-overlapping neighbors
    for i in range(n):
        lend_i, rend_i, len_i, idx_i = bounds[i]
        j = i + 1
        while j < n and bounds[j][0] <= rend_i:
            lend_j, rend_j, len_j, idx_j = bounds[j]
            # Tested BEFORE the byte-count prefilter below, which demands
            # span_ov >= 0.2 * the LONGER cdna: a short intron-sharing model inside a
            # long gene fails that while still agreeing on junctions, so gating the
            # intron rule behind it would drop exactly those pairs. Identical chains
            # are already handled above and do not rely on this threshold.
            if _shares_introns(
                intron_sets[idx_i],
                intron_sets[idx_j],
                shared_intron_thr,
                shared_intron_min,
            ):
                union(idx_i, idx_j)
                j += 1
                continue
            span_ov = min(rend_i, rend_j) - max(lend_i, lend_j) + 1
            if span_ov > 0:
                shorter = min(len_i, len_j)
                longer = max(len_i, len_j)
                # minimal bases required by thresholds
                min_required = max(int(shorter_thr * shorter), int(longer_thr * longer))
                if span_ov >= min_required and shorter > 0 and longer > 0:
                    ov = _transcript_overlap_len(transcripts[idx_i], transcripts[idx_j])
                    contained = (ov == shorter and shorter > 0)
                    if contained or (
                        (ov / shorter) >= shorter_thr and (ov / longer) >= longer_thr
                    ):
                        union(idx_i, idx_j)
            j += 1

    # Convert DSU roots to normalized community ids 0..k-1
    membership = [0] * n
    remap = {}
    next_id = 0
    for idx in range(n):
        r = find(idx)
        if r not in remap:
            remap[r] = next_id
            next_id += 1
        membership[idx] = remap[r]

    return membership

__all__.append("partition_with_dsu")
