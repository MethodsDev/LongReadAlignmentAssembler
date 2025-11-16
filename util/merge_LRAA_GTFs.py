#!/usr/bin/env python3

import sys, os, re

sys.path.insert(
    0, os.path.sep.join([os.path.dirname(os.path.realpath(__file__)), "../pylib"])
)

from Splice_graph import Splice_graph  # type: ignore
from Transcript import Transcript, GTF_contig_to_transcripts  # type: ignore
from LRAA import LRAA  # type: ignore
import LRAA_Globals  # type: ignore
import logging
import traceback
import argparse
from collections import defaultdict
import Util_funcs  # type: ignore

FORMAT = (
    "%(asctime)-15s %(levelname)s %(module)s.%(name)s.%(funcName)s:\n\t%(message)s\n"
)

logger = logging.getLogger()

# Reconfigure logging to ensure visible output on terminal even if other modules already set handlers.
def _configure_logging(debug=False):
    # Remove any pre-existing handlers to avoid duplicate or suppressed output
    try:
        for h in list(logger.handlers):
            logger.removeHandler(h)
    except Exception:
        pass

    logger.setLevel(logging.DEBUG if debug else logging.INFO)

    # Stream handler for stdout (INFO and below)
    class StdoutFilter(logging.Filter):
        def filter(self, record):
            return record.levelno <= logging.INFO

    sh_out = logging.StreamHandler(stream=sys.stdout)
    sh_out.setLevel(logging.DEBUG)
    sh_out.addFilter(StdoutFilter())
    sh_out.setFormatter(logging.Formatter(FORMAT))
    logger.addHandler(sh_out)

    # Stream handler for stderr (WARNING and above)
    sh_err = logging.StreamHandler(stream=sys.stderr)
    sh_err.setLevel(logging.WARNING)
    sh_err.setFormatter(logging.Formatter(FORMAT))
    logger.addHandler(sh_err)

    # Basic exception hook to ensure uncaught exceptions print traceback to stderr
    def _excepthook(exc_type, exc, tb):
        try:
            logger.error("Uncaught exception: %s", exc)
            traceback.print_exception(exc_type, exc, tb, file=sys.stderr)
        except Exception:
            pass
    sys.excepthook = _excepthook


_configure_logging(debug=False)


def main():

    parser = argparse.ArgumentParser(
        description="Merge LRAA gtfs",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )

    parser.add_argument(
        "--gtf",
        type=str,
        required=True,
        nargs="+",
        help="LRAA gtfs to merge",
    )

    parser.add_argument("--genome", type=str, required=True, help="target genome file")

    parser.add_argument(
        "--output_gtf",
        type=str,
        default="LRAA.merged.gtf",
        help="prefix for output filenames",
    )

    parser.add_argument(
        "--exclude_SE_transcripts",
        action="store_true",
        default=False,
        help="exclude single-exon isoforms from the merge",
    )

    parser.add_argument(
        "--debug",
        "-d",
        action="store_true",
        default=False,
        help="debug mode, more verbose",
    )

    # Community clustering controls (mirrors main LRAA CLI)
    parser.add_argument(
        "--no_use_community_clustering",
        action="store_true",
        default=False,
        help=(
            "Disable Leiden community clustering for transcript→gene assignment (enabled by default)."
        ),
    )

    parser.add_argument(
        "--community_resolution",
        type=float,
        default=1.0,
        help=(
            "Leiden resolution parameter (higher → more communities)."
        ),
    )

    parser.add_argument(
        "--community_random_seed",
        type=int,
        default=42,
        help=(
            "Random seed for Leiden community clustering."
        ),
    )

    parser.add_argument(
        "--store_backend",
        default=LRAA_Globals.config.get("store_backend", "memory"),
        choices=["auto", "lmdb", "sqlite", "memory"],
        help=(
            "Storage backend for read tracking stores (ReadNameStore/MpReadIdStore). "
            "Choices: auto (prefer lmdb, fallback sqlite), lmdb, sqlite, memory."
        ),
    )

    args = parser.parse_args()

    if args.debug:
        LRAA_Globals.DEBUG = True
        _configure_logging(debug=True)
        logger.debug("Debug logging enabled for merge script.")

    LRAA_Globals.LRAA_MODE = "MERGE"

    # Set store backend behavior and propagate to env so helpers honor it
    resolved_backend = (args.store_backend or "memory").lower()
    LRAA_Globals.config["store_backend"] = resolved_backend
    try:
        if resolved_backend == "memory":
            os.environ["LRAA_READSTORE_BACKEND"] = "memory"
        elif resolved_backend == "sqlite":
            os.environ["LRAA_READSTORE_BACKEND"] = "sqlite"
        elif resolved_backend == "lmdb":
            os.environ["LRAA_READSTORE_BACKEND"] = "lmdb"
            try:
                import lmdb  # type: ignore  # noqa: F401
            except Exception:
                sys.exit("Error: --store_backend lmdb requested but python-lmdb is not installed.")
        else:  # auto
            os.environ.pop("LRAA_READSTORE_BACKEND", None)
    except Exception:
        pass

    logger.info("merge store backend: %s", resolved_backend)

    # Map community clustering flags into global config
    try:
        LRAA_Globals.config["use_community_clustering"] = not args.no_use_community_clustering
        LRAA_Globals.config["community_resolution"] = float(args.community_resolution)
        LRAA_Globals.config["community_random_seed"] = int(args.community_random_seed)
        logger.debug(
            "Community clustering settings (merge): use=%s, resolution=%.3f, seed=%d",
            LRAA_Globals.config["use_community_clustering"],
            LRAA_Globals.config["community_resolution"],
            LRAA_Globals.config["community_random_seed"],
        )
    except Exception as _e:
        logger.warning(f"Failed to apply community clustering config overrides: {_e}")

    gtf_list = args.gtf
    genome_fasta_file = args.genome
    output_gtf = args.output_gtf
    exclude_SE = args.exclude_SE_transcripts

    if len(gtf_list) < 2 and not LRAA_Globals.DEBUG:
        exit("Error, need at least two gtf files to merge")

    ofh = open(output_gtf, "wt")
    tracking_records = []  # accumulate provenance rows across contigs/strands

    contig_strand_to_input_transcripts = defaultdict(list)

    for gtf_file in gtf_list:
        logger.info(f"-capturing input transcripts from gtf {gtf_file}")
        contig_to_input_transcripts = (
            GTF_contig_to_transcripts.parse_GTF_to_Transcripts(gtf_file)
        )
        for contig, transcript_obj_list in contig_to_input_transcripts.items():
            for transcript in transcript_obj_list:
                if exclude_SE and not transcript.has_introns():
                    logger.debug("-excluding SE isoform: {}".format(transcript))
                    continue

                # annotate source GTF basename for later provenance encoding
                try:
                    src_basename = os.path.basename(gtf_file)
                    transcript.add_meta("source_gtf", src_basename)
                except Exception:
                    pass

                transcript_strand = transcript.get_strand()
                contig_strand_token = "{}^{}".format(contig, transcript_strand)
                contig_strand_to_input_transcripts[contig_strand_token].append(
                    transcript
                )

    for stranded_contig, transcript_list in contig_strand_to_input_transcripts.items():

        contig_acc, contig_strand = stranded_contig.split("^")

        contig_seq_str = Util_funcs.retrieve_contig_seq_from_fasta_file(
            contig_acc, genome_fasta_file
        )

        ## Build Splice Graph
        logger.info(f"\n// -building splice graph for {contig_acc}")
        sg = Splice_graph()
        sg.build_splice_graph_for_contig(
            contig_acc,
            contig_strand,
            contig_seq_str,
            alignments_bam_file=None,
            region_lend=None,
            region_rend=None,
            input_transcripts=transcript_list,
            quant_mode=False,
        )

        lraa_obj = LRAA(sg)

        lraa_obj.assign_transcripts_paths_in_graph(transcript_list)

        lraa_obj.build_multipath_graph(
            contig_acc,
            contig_strand,
            contig_seq_str,
            bam_file=None,
            input_transcripts=transcript_list,
        )

        # Define merged isoforms
        logger.info(f"\n// -begin merge of isoforms for {contig_acc}")
        transcripts = lraa_obj.reconstruct_isoforms()

        # Optional final reclustering/refinement + reporting
        # This will also trigger neighbor-Jaccard pair reporting if configured
        try:
            transcripts = Transcript.recluster_transcripts_to_genes(
                transcripts, contig_acc, contig_strand
            )
        except Exception as e:
            # fallback to original transcripts if reclustering fails for any reason
            logger.warning(
                f"Reclustering/refinement skipped due to error: {e}. Proceeding with original transcripts."
            )

        ## report transcripts in GTF format
        logger.info(
            "writing gtf output for {} [{}] containing {} transcripts".format(
                contig_acc, contig_strand, len(transcripts)
            )
        )

        for transcript in transcripts:
            ofh.write(transcript.to_GTF_format(include_TPM=False) + "\n")

            # Build tracking provenance rows from synthetic read names supporting this transcript
            try:
                # gather all names across assigned multipath evidences (usually one)
                name_to_weight = dict()  # key: (src, tid) -> { 'TSS':bool, 'PolyA':bool, 'n':int }
                for mp in transcript.get_multipaths_evidence_assigned():
                    for rn in mp.get_read_names():
                        if not isinstance(rn, str):
                            continue
                        if not rn.startswith("fake_for_merge|"):
                            continue
                        # parse tokens like key=value separated by |
                        try:
                            parts = rn.split("|")
                            kv = {}
                            for p in parts[1:]:  # skip prefix
                                if "=" in p:
                                    k, v = p.split("=", 1)
                                    kv[k] = v
                            src = kv.get("src", "unknown")
                            tid = kv.get("tid", "unknown")
                            flags_str = kv.get("flags", "-")
                            n_str = kv.get("n", "1")
                            has_TSS = 1 if ("TSS" in flags_str.split(",")) else 0
                            has_PolyA = 1 if ("PolyA" in flags_str.split(",")) else 0
                            n = int(n_str) if n_str.isdigit() else 1
                            key = (src, tid)
                            prev = name_to_weight.get(key)
                            if prev is None:
                                name_to_weight[key] = {
                                    "TSS": has_TSS,
                                    "PolyA": has_PolyA,
                                    "n": n,
                                }
                            else:
                                # ensure flags reflect any True seen; keep n as the maximum encountered
                                prev["TSS"] = 1 if (prev["TSS"] or has_TSS) else 0
                                prev["PolyA"] = 1 if (prev["PolyA"] or has_PolyA) else 0
                                prev["n"] = max(prev["n"], n)
                        except Exception:
                            continue

                for (src, tid), info in name_to_weight.items():
                    tracking_records.append(
                        {
                            "merged_transcript_id": transcript.get_transcript_id(),
                            "merged_gene_id": transcript.get_gene_id(),
                            "contig": contig_acc,
                            "strand": contig_strand,
                            "source_gtf": src,
                            "source_transcript_id": tid,
                            "source_has_TSS": info["TSS"],
                            "source_has_PolyA": info["PolyA"],
                            # contribution_count intentionally omitted from output (was info["n"]) as it's configuration-dependent
                        }
                    )
            except Exception:
                pass

    logger.info("Done.")

    # Write tracking output next to the GTF
    try:
        tracking_path = output_gtf + ".tracking.tsv"
        with open(tracking_path, "wt") as tfh:
            header = [
                "merged_transcript_id",
                "merged_gene_id",
                "contig",
                "strand",
                "source_gtf",
                "source_transcript_id",
                "source_has_TSS",
                "source_has_PolyA",
            ]
            tfh.write("\t".join(header) + "\n")
            for rec in tracking_records:
                row = [
                    rec["merged_transcript_id"],
                    rec["merged_gene_id"],
                    rec["contig"],
                    rec["strand"],
                    rec["source_gtf"],
                    rec["source_transcript_id"],
                    str(rec["source_has_TSS"]),
                    str(rec["source_has_PolyA"]),
            # contribution_count removed
                ]
                tfh.write("\t".join(row) + "\n")
        logger.info(f"Wrote tracking file: {tracking_path}")
    except Exception as _e:
        logger.warning(f"Failed to write tracking TSV due to: {_e}")

    sys.exit(0)


if __name__ == "__main__":
    main()
