#!/usr/bin/env python3

import sys, re, os
from collections import defaultdict
from GenomeFeature import GenomeFeature
import Scored_path
import LRAA_Globals
import logging
import gzip

logger = logging.getLogger(__name__)


class Transcript(GenomeFeature):

    trans_id_counter = 0
    meta_ignore = ["gene_id", "transcript_id", "TSS", "PolyA"]

    def __init__(self, contig_acc, segment_coordinates_list, orient):

        segment_coordinates_list = sorted(segment_coordinates_list, key=lambda x: x[0])
        segment_coordinates_list = _merge_short_deletions(
            segment_coordinates_list, LRAA_Globals.config["read_aln_gap_merge_int"]
        )

        trans_lend = segment_coordinates_list[0][0]
        trans_rend = segment_coordinates_list[-1][1]

        super().__init__(contig_acc, trans_lend, trans_rend, orient)

        self._orient = orient
        self._exon_segments = segment_coordinates_list

        Transcript.trans_id_counter += 1

        self._id = ":".join(
            [
                contig_acc,
                str(trans_lend),
                str(trans_rend),
                "N:{}".format(Transcript.trans_id_counter),
                orient,
            ]
        )

        self._transcript_id = "t.{}".format(
            self._id
        )  # should reset to something more useful
        self._gene_id = "g.{}".format(self._id)  # ditto above

        self._meta = dict()

        # can import features from GTF feature attributes
        self._imported_TPM_val = None
        self._imported_has_TSS = None  # if parsed info from gtf, set True/False
        self._imported_has_POLYA = None

        self._likely_internal_primed = None

        self._scored_path_obj = (
            None  # optional - useful if transcript obj was built based on a scored path
        )

        self.multipaths_evidence_assigned = (
            set()
        )  # multipaths supporting the transcript structure.

        self._multipaths_evidence_weights = dict()

        self._multipath = (
            None  # multipath obj (when transcript is constructed based on a multipath)
        )
        self._simplepath = (
            None  # multipath obj (when transcript is constructed based on a multipath)
        )

        self._read_counts_assigned = None  # set during expression quantification

        self._isoform_fraction = None  # set during expression quantification

        self._is_novel_isoform_bool = True  # set to False if it's a known isoform.

        self._cdna_len = 0
        for exon_segment in segment_coordinates_list:
            self._cdna_len += exon_segment[1] - exon_segment[0] + 1

        return

    ## getters

    def get_exon_segments(self):
        return self._exon_segments.copy()

    def get_meta(self):
        """
        Returns a shallow copy of the transcript metadata dictionary.
        Useful for extracting provenance like source GTF or other annotations
        added during parsing.
        """
        return dict(self._meta) if self._meta is not None else {}

    def get_num_exon_segments(self):
        return len(self._exon_segments)

    def is_monoexonic(self):
        return self.get_num_exon_segments() == 1

    def has_introns(self):
        return len(self.get_exon_segments()) > 1

    def get_introns(self):
        intron_coordsets = list()
        exon_segments = self.get_exon_segments()
        if len(exon_segments) > 1:
            exon_segments = sorted(exon_segments, key=lambda x: x[0])
            for i in range(1, len(exon_segments)):
                intron_lend = exon_segments[i - 1][1] + 1
                intron_rend = exon_segments[i][0] - 1
                assert intron_lend < intron_rend
                intron_coordsets.append((intron_lend, intron_rend))
        return intron_coordsets

    def get_exons_string(self):
        return "{}:({}){}".format(
            self.get_contig_acc(), self.get_strand(), self.get_exon_segments()
        )

    def get_introns_string(self):
        return "{}:({}){}".format(
            self.get_contig_acc(), self.get_strand(), self.get_introns()
        )

    def get_strand(self):
        return self._orient

    def get_orient(self):
        return self.get_strand()

    def get_cdna_len(self):
        return self._cdna_len

    def get_transcript_id(self):
        if self._transcript_id is not None:
            return self._transcript_id

        elif "transcript_id" in self._meta:
            return self._meta["transcript_id"]

        else:
            return self._id

    def get_gene_name(self):
        if "gene_name" in self._meta:
            return self._meta["gene_name"]
        else:
            return None

    def get_gene_id(self):
        if self._gene_id is not None:
            return self._gene_id

        elif "gene_id" in self._meta:
            return self._meta["gene_id"]
        else:
            raise RuntimeError("gene_id not set for Transcript obj")

    def get_simple_path(self):
        assert self._simplepath is not None
        assert len(self._simplepath) > 0
        return self._simplepath

    def set_simple_path(self, simple_path):
        assert simple_path is not None
        assert len(simple_path) > 0
        self._simplepath = simple_path

    def get_left_boundary_sort_weight(self):
        assert self._simplepath is not None
        if re.match("TSS:|POLYA:", self._simplepath[0]):
            return 1
        else:
            return 0

    def get_right_boundary_sort_weight(self):
        assert self._simplepath is not None
        if re.match("TSS:|POLYA:", self._simplepath[-1]):
            return 1
        else:
            return 0

    def has_annotated_PolyA(self):
        if self._imported_has_POLYA is not None:
            return self._imported_has_POLYA
        else:
            return False

    def has_PolyA(self):

        if self._imported_has_POLYA is not None:
            return self._imported_has_POLYA

        assert self._simplepath is not None
        if re.match("POLYA:", self._simplepath[0]) or re.match(
            "POLYA:", self._simplepath[-1]
        ):
            return True
        else:
            return False

    def has_annotated_TSS(self):
        if self._imported_has_TSS is not None:
            return self._imported_has_TSS
        else:
            return False

    def has_TSS(self):

        if self._imported_has_TSS is not None:
            return self._imported_has_TSS

        assert self._simplepath is not None
        if re.match("TSS:", self._simplepath[0]) or re.match(
            "TSS:", self._simplepath[-1]
        ):
            return True
        else:
            return False


    def set_likely_internal_primed(self, TorF_boolean):
        # indicate if the transcript looks like it's internally primed.
        self._likely_internal_primed = TorF_boolean

    def __repr__(self):

        text = "Transcript: {} {}-{} [{}] {} {} segs: {}".format(
            self._contig_acc,
            self._lend,
            self._rend,
            self._orient,
            self.get_transcript_id(),
            self.get_gene_id(),
            self._exon_segments,
        )
        if self._simplepath is not None:
            text += " {} ".format(self._simplepath)

        if self._meta is not None and len(self._meta) > 0:
            text += "\t" + str(self._meta)

        return text

    def set_scored_path_obj(self, scored_path_obj):
        assert type(scored_path_obj) == Scored_path.Scored_path
        self._scored_path_obj = scored_path_obj

    def set_gene_id(self, gene_id):
        self._gene_id = gene_id

    def set_transcript_id(self, transcript_id):
        self._transcript_id = transcript_id

    def add_meta(self, meta_key, meta_val=None):

        if self._meta == None:
            self._meta = dict()

        if meta_val is None and type(meta_key) == dict:
            self._meta = meta_key.copy()

        elif meta_val is not None:
            self._meta[meta_key] = meta_val
        else:
            raise RuntimeError("Error: not sure how to handle input params")

        return

    def get_multipaths_evidence_assigned(self):
        return self.multipaths_evidence_assigned.copy()

    def includes_reference_transcript(self):
        return len(self.get_ref_trans_included()) > 0

    def get_ref_trans_included(self):
        ref_trans_included = set()
        for mp in self.get_multipaths_evidence_assigned():
            for read_name in mp.get_read_names():
                if "reftranscript:" in read_name:
                    ref_trans_included.add(read_name)
        return ref_trans_included

    def is_novel_isoform(self):
        return self._is_novel_isoform_bool

    def set_is_novel_isoform(self, boolean_val: bool):
        self._is_novel_isoform_bool = boolean_val

    def add_multipaths_evidence_assigned(self, multipaths):
        if self.multipaths_evidence_assigned == None:
            self.multipaths_evidence_assigned = set()

        if type(multipaths) in (list, set):
            self.multipaths_evidence_assigned.update(multipaths)

        else:
            self.multipaths_evidence_assigned.add(multipaths)

        self._ensure_multipaths_have_weights()

    def set_multipaths_evidence_weights(self, mp_weights: dict):
        for mp, weight in mp_weights.items():
            self._multipaths_evidence_weights[mp] = weight

    def _ensure_multipaths_have_weights(self):
        for mp in self.multipaths_evidence_assigned:
            if mp not in self._multipaths_evidence_weights:
                self._multipaths_evidence_weights[mp] = 1.0

    def get_multipath_weight(self, multipath):
        assert (
            multipath in self._multipaths_evidence_weights
        ), "Error, not finding multipath {} in mp weights for transcript {}".format(
            multipath, self.get_transcript_id()
        )
        return self._multipaths_evidence_weights[multipath]

    def has_multipath_as_evidence(self, multipath):
        return multipath in self.multipaths_evidence_assigned

    def absorb_other_transcript_multipaths(self, other_transcript):

        assert isinstance(other_transcript, Transcript)

        other_transcript_mps = other_transcript.get_multipaths_evidence_assigned()
        for other_mp in other_transcript_mps:
            other_mp_weight = other_transcript.get_multipath_weight(other_mp)

            if self.has_multipath_as_evidence(other_mp):
                this_mp_weight = self.get_multipath_weight(other_mp)
                if this_mp_weight < other_mp_weight:
                    # use higher weight
                    self.set_multipaths_evidence_weights({other_mp: other_mp_weight})
            else:
                # add as new evidence
                self.set_multipaths_evidence_weights({other_mp: other_mp_weight})
                self.add_multipaths_evidence_assigned(other_mp)

    def prune_reftranscript_as_evidence(self):

        mps = self.get_multipaths_evidence_assigned()
        for mp in mps:
            mp.prune_reftranscript_as_evidence()

    def set_read_counts_assigned(self, read_counts):
        self._read_counts_assigned = read_counts

    def get_read_counts_assigned(self):

        if self._imported_TPM_val is not None:
            return self._imported_TPM_val

        assert (
            self._read_counts_assigned is not None
        ), "Error, read counts assigned is None - maybe quant not run yet? " + str(self)
        return self._read_counts_assigned

    def has_annotated_TPM(self):
        if self._imported_TPM_val is not None:
            return True
        else:
            return False

    def get_TPM(self):
        if self._imported_TPM_val is not None:
            return self._imported_TPM_val
        else:
            return self.get_expr_fraction() * 1e6

    def get_expr_fraction(self):
        return self.get_read_counts_assigned() / LRAA_Globals.config["num_total_reads"]

    def set_isoform_fraction(self, frac):
        self._isoform_fraction = frac

    def get_isoform_fraction(self):
        assert (
            self._isoform_fraction is not None
        ), "Error, isoform fraction is None - maybe quant not run yet?"
        return self._isoform_fraction

    def lighten(self):
        # lighten transcript by removing nonessential memory allocs
        self._scored_path_obj = None
        self._multipath = None
        self.init_quant_info()

    def to_GTF_format(self, include_TPM=False):

        ## transcript line:

        gtf_text = ""

        gtf_text += "\t".join(
            [
                self._contig_acc,
                "LRAA",
                "transcript",
                str(self._lend),
                str(self._rend),
                ".",
                self._orient,
                ".",
                'gene_id "{}"; transcript_id "{}";'.format(
                    self.get_gene_id(), self.get_transcript_id()
                ),
            ]
        )

        if self._meta:
            for meta_key in sorted(self._meta.keys()):
                if meta_key not in Transcript.meta_ignore:
                    gtf_text += ' {} "{}";'.format(meta_key, self._meta[meta_key])

        # include other transcript features
        misc_transcript_features = dict()
        if include_TPM:
            misc_transcript_features["TPM"] = "{:.3f}".format(self.get_TPM())

        misc_transcript_features["PolyA"] = str(self.has_PolyA())
        misc_transcript_features["TSS"] = str(self.has_TSS())

        # Internal priming annotation: prefer internal flag, else fallback to imported meta if present
        if self._likely_internal_primed is not None:
            misc_transcript_features["InternalPriming"] = str(self._likely_internal_primed)
        elif self._meta and "InternalPriming" in self._meta:
            # ensure we still propagate an imported value even if the internal flag wasn't explicitly set
            misc_transcript_features["InternalPriming"] = str(self._meta["InternalPriming"])

        for misc_feature, misc_val in misc_transcript_features.items():
            gtf_text += ' {} "{}";'.format(misc_feature, misc_val)

        gtf_text += "\n"

        for segment in self._exon_segments:
            gtf_text += (
                "\t".join(
                    [
                        self._contig_acc,
                        "LRAA",
                        "exon",
                        str(segment[0]),
                        str(segment[1]),
                        ".",
                        self._orient,
                        ".",
                        'gene_id "{}"; transcript_id "{}";'.format(
                            self.get_gene_id(), self.get_transcript_id()
                        ),
                    ]
                )
                + "\n"
            )

        if self._scored_path_obj:
            # include construction info as comments.
            gtf_text += "# derived from scored path obj:\n"

            mpgn_list = self._scored_path_obj.get_path_mpgn_list()
            for mpgn in mpgn_list:
                gtf_text += "# " + str(mpgn) + "\n"

            gtf_text += "\n"

        return gtf_text

    def init_quant_info(self):
        self.multipaths_evidence_assigned = set()
        self._multipaths_evidence_weights = dict()
        self._read_counts_assigned = None
        self._isoform_fraction = None

        return


    @classmethod
    def recluster_transcripts_to_genes(cls, transcripts, contig_acc, contig_strand):
        import networkx as nx
        from typing import Hashable, List, Dict

        # Build first-pass graph linking transcripts sharing any simple path node
        G = nx.Graph()
        G.add_nodes_from(range(len(transcripts)))
        seen: Dict[Hashable, int] = {}
        for i, t in enumerate(transcripts):
            sp = t.get_simple_path()
            for nid in set(sp):
                if nid in seen:
                    G.add_edge(seen[nid], i)
                else:
                    seen[nid] = i

        initial_clusters = [
            [transcripts[i] for i in comp] for comp in nx.connected_components(G)
        ]

        # Optional debug export of initial cluster memberships
        try:
            debug_init_clusters = LRAA_Globals.config.get("debug_write_init_clusters", True)
        except Exception:
            debug_init_clusters = True
        if debug_init_clusters:
            # Write TSV: cluster_index (1-based) \t size \t transcript_ids (comma-separated)
            with open("__init_clusters.tsv", "wt") as fh_dbg:
                print("cluster_index\tsize\ttranscript_ids", file=fh_dbg)
                for idx, clust in enumerate(initial_clusters, start=1):
                    tids = [t.get_transcript_id() for t in clust]
                    print(f"{idx}\t{len(clust)}\t{','.join(tids)}", file=fh_dbg)



        # Overlap thresholds
        min_overlap_shorter_frac = LRAA_Globals.config.get(
            "min_recluster_overlap_shorter_iso_frac", 0.50
        )
        min_overlap_longer_frac = LRAA_Globals.config.get(
            "min_recluster_overlap_longer_iso_frac", 0.10
        )

        def transcript_overlap_len(t1: 'Transcript', t2: 'Transcript') -> int:
            def merge_intervals(exons):
                merged = []
                for a, b in sorted(exons):
                    if not merged or a > merged[-1][1] + 1:
                        merged.append([a, b])
                    else:
                        if b > merged[-1][1]:
                            merged[-1][1] = b
                return merged

            m1 = merge_intervals(t1.get_exon_segments())
            m2 = merge_intervals(t2.get_exon_segments())
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

        def transcript_len(t: 'Transcript') -> int:
            return t.get_cdna_len()

        refined_clusters: List[List['Transcript']] = []
        # Progress bar (tqdm) expected to be available
        from tqdm import tqdm  # type: ignore
        cluster_iter = tqdm(
            initial_clusters,
            desc=f"Reclustering transcripts {contig_acc}({contig_strand})",
            unit="cluster",
        )

        for cluster in cluster_iter:
            if logger.isEnabledFor(logging.INFO):
                logger.info(
                    "recluster: processing initial cluster of %d transcripts (contig=%s strand=%s)",
                    len(cluster),
                    contig_acc,
                    contig_strand,
                )
            use_comm = bool(LRAA_Globals.config.get("use_community_clustering", False))
            max_for_leiden = int(
                LRAA_Globals.config.get("max_transcripts_for_community_clustering", 1500)
            )
            if len(cluster) == 1:
                refined_clusters.append(cluster)
            elif use_comm and len(cluster) <= max_for_leiden:
                if logger.isEnabledFor(logging.INFO):
                    logger.info(
                        "recluster: using Leiden for cluster of %d transcripts (threshold=%d)",
                        len(cluster),
                        max_for_leiden,
                    )
                # Leiden community clustering branch (fatal on failure)
                from GeneCommunityCluster import partition_with_leiden
                membership = partition_with_leiden(
                    cluster,
                    contig_acc,
                    contig_strand,
                    resolution=LRAA_Globals.config.get("community_resolution", 1.0),
                    seed=LRAA_Globals.config.get("community_random_seed", 42),
                )
                comm_to_indices = defaultdict(list)
                for idx, cid in enumerate(membership):
                    comm_to_indices[cid].append(idx)
                for cid in sorted(comm_to_indices.keys()):
                    comp_transcripts = [cluster[i] for i in comm_to_indices[cid]]
                    refined_clusters.append(comp_transcripts)
            else:
                if use_comm and len(cluster) > max_for_leiden:
                    logger.info(
                        "recluster: skipping Leiden due to size (%d > %d); using DSU fallback",
                        len(cluster),
                        max_for_leiden,
                    )
                elif not use_comm:
                    logger.info(
                        "recluster: community clustering disabled; using DSU fallback for cluster of %d",
                        len(cluster),
                    )
                # Use modular DSU-based partitioning from GeneCommunityCluster
                from GeneCommunityCluster import partition_with_dsu
                membership = partition_with_dsu(cluster, contig_acc, contig_strand)
                comp_to_indices = defaultdict(list)
                for idx, cid in enumerate(membership):
                    comp_to_indices[cid].append(idx)
                if logger.isEnabledFor(logging.INFO):
                    logger.info(
                        "recluster: DSU fallback produced %d subcomponents",
                        len(comp_to_indices),
                    )
                for cid in sorted(comp_to_indices.keys()):
                    comp_transcripts = [cluster[i] for i in comp_to_indices[cid]]
                    refined_clusters.append(comp_transcripts)

        revised_transcripts = []
        for i, cluster in enumerate(refined_clusters):
            new_gene_id = f"g:{contig_acc}:{contig_strand}:comp-{i+1}"
            for j, t in enumerate(cluster):
                new_transcript_id = f"t:{contig_acc}:{contig_strand}:comp-{i+1}:iso-{j+1}"
                t.set_gene_id(new_gene_id)
                t.set_transcript_id(new_transcript_id)
                revised_transcripts.append(t)


        return revised_transcripts

class GTF_contig_to_transcripts:

    @classmethod
    def parse_GTF_to_Transcripts(
        cls,
        gtf_filename,
        chr_restrict=None,
        strand_restrict=None,
        lend_restrict=None,
        rend_restrict=None,
    ):

        gene_id_to_meta = defaultdict(dict)
        transcript_id_to_meta = defaultdict(dict)
        transcript_id_to_genome_info = defaultdict(dict)

        local_debug = False

        opener = gzip.open if re.search("\\.gz$", gtf_filename) is not None else open

        with opener(gtf_filename, "rt") as fh:
            for line in fh:
                if line[0] == "#":
                    continue
                line = line.rstrip()

                if not re.match("\\w", line):
                    continue
                vals = line.split("\t")
                if len(vals) < 9:
                    logger.warn(
                        "GTF line has fewer fields than expected. Skipping. {}".format(
                            line
                        )
                    )
                    continue

                contig = vals[0]
                feature_type = vals[2]
                lend = int(vals[3])
                rend = int(vals[4])
                strand = vals[6]

                if chr_restrict is not None:
                    if chr_restrict != contig:
                        if local_debug:
                            logger.debug(
                                "skipping {} since contig {} != restricted to {}".format(
                                    line, contig, chr_restrict
                                )
                            )
                        continue

                if strand_restrict is not None:
                    if strand != strand_restrict:
                        if local_debug:
                            logger.debug(
                                "skipping {} since strand {} != restricted to {}".format(
                                    line, strand, strand_restrict
                                )
                            )
                        continue

                if lend_restrict is not None and rend_restrict is not None:
                    if lend < lend_restrict or rend > rend_restrict:
                        if local_debug:
                            logger.debug(
                                "skipping {} as lend {} and rend {} are not within range {}-{}".format(
                                    line, lend, rend, lend_restrict, rend_restrict
                                )
                            )
                        continue

                info = vals[8]
                info_dict = cls._parse_info(info)

                if feature_type == "gene":
                    gene_id = info_dict["gene_id"]
                    gene_id_to_meta[gene_id] = info_dict

                if "transcript_id" not in info_dict:
                    if local_debug:
                        logger.debug(
                            "skipping {} as transcript_id not in info_dict".format(line)
                        )
                    continue

                transcript_id = info_dict["transcript_id"]

                if transcript_id not in transcript_id_to_meta:
                    transcript_id_to_meta[transcript_id] = info_dict
                else:
                    transcript_id_to_meta[transcript_id].update(info_dict)

                if feature_type != "exon":
                    continue

                transcript_id_to_genome_info[transcript_id]["contig"] = contig
                transcript_id_to_genome_info[transcript_id]["strand"] = strand
                if "coords" in transcript_id_to_genome_info[transcript_id].keys():
                    transcript_id_to_genome_info[transcript_id]["coords"].append(
                        [lend, rend]
                    )
                else:
                    transcript_id_to_genome_info[transcript_id]["coords"] = [
                        [lend, rend]
                    ]

        # convert to transcript objects

        contig_to_transcripts = defaultdict(list)

        if LRAA_Globals.DEBUG:
            debug_fh = open("__transcript_gtf_features_parsed.tsv", "wt")

        for transcript_id in transcript_id_to_genome_info:
            transcript_info_dict = transcript_id_to_genome_info[transcript_id]
            contig = transcript_info_dict["contig"]
            strand = transcript_info_dict["strand"]
            coords_list = transcript_info_dict["coords"]

            transcript_meta = transcript_id_to_meta[transcript_id]
            gene_id = transcript_meta["gene_id"]
            gene_meta = gene_id_to_meta[gene_id]
            transcript_meta.update(gene_meta)

            # print("Transcript meta: {}".format(str(transcript_meta)))

            transcript_obj = Transcript(contig, coords_list, strand)
            transcript_obj.add_meta(transcript_meta)

            transcript_obj.set_gene_id(gene_id)
            transcript_obj.set_transcript_id(transcript_id)

            if "PolyA" in transcript_meta:
                transcript_obj._imported_has_POLYA = (
                    True if transcript_meta["PolyA"] == "True" else False
                )
            else:
                transcript_obj._imported_has_POLYA = False

            if "TSS" in transcript_meta:
                transcript_obj._imported_has_TSS = (
                    True if transcript_meta["TSS"] == "True" else False
                )
            else:
                transcript_obj._imported_has_TSS = False

            if "TPM" in transcript_meta:
                transcript_obj._imported_TPM_val = float(transcript_meta["TPM"])

            # import InternalPriming flag if present so it will be re-exported via to_GTF_format
            if "InternalPriming" in transcript_meta:
                transcript_obj._likely_internal_primed = (
                    True if transcript_meta["InternalPriming"] == "True" else False
                )

            contig_to_transcripts[contig].append(transcript_obj)

            if LRAA_Globals.DEBUG:
                print(
                    "\t".join(
                        [contig, gene_id, transcript_id, strand, str(coords_list)]
                    ),
                    file=debug_fh,
                )

        if LRAA_Globals.DEBUG:
            debug_fh.close()

        return contig_to_transcripts

    # private
    @classmethod
    def _parse_info(cls, info):

        info_dict = dict()

        parts = info.split(";")
        for part in parts:
            part = part.strip()
            m = re.match('^(\\S+) \\"([^\\"]+)\\"', part)
            if m is not None:
                token = m.group(1)
                val = m.group(2)

                info_dict[token] = val

        return info_dict


# utility functions


def _merge_short_deletions(segments_list, merge_dist):

    if len(segments_list) == 1:
        return segments_list  # nothing to do here.

    merged_segments_list = list()
    merged_segments_list.append(segments_list[0].copy())

    for segment in segments_list[1:]:
        if segment[0] - merged_segments_list[-1][1] <= merge_dist:
            # merge
            merged_segments_list[-1][1] = segment[1]
        else:
            merged_segments_list.append(segment.copy())

    return merged_segments_list


if __name__ == "__main__":

    # testing gtf parser
    usage = "usage: {} gtf_filename\n\n".format(sys.argv[0])

    if len(sys.argv) < 2:
        exit(usage)

    gtf_filename = sys.argv[1]

    contig_to_transcripts = GTF_contig_to_transcripts.parse_GTF_to_Transcripts(
        gtf_filename
    )

    for contig, transcript_list in contig_to_transcripts.items():
        for transcript_obj in transcript_list:
            print("\t".join([contig, str(transcript_obj)]))

    sys.exit(0)
