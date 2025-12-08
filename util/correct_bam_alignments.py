#!/usr/bin/env python3
# encoding: utf-8

"""
Correct BAM alignments using splice graph-guided soft-clip realignment.

This utility leverages LRAA's alignment correction machinery to fix soft-clipped
alignments based on a reference GTF and genome. It builds a splice graph from
normalized reads and reference transcripts, then uses it to correct alignments
where soft-clipped sequences can be realigned across known splice junctions.

Outputs:
- Full BAM with corrections applied
- BAM containing only original (faulty) alignments that were corrected
- BAM containing only the corrected versions of those alignments
- BAM containing all original alignments, with tags marking correctable ones
- Statistics report on correction rates per contig/strand
"""

import sys
import os
import re
import argparse
import logging
import subprocess
import tempfile
import shutil
import json
import hashlib
from pathlib import Path
from collections import defaultdict
import multiprocessing as mp
import pysam

# Add pylib to path
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
PYLIB_DIR = os.path.join(os.path.dirname(SCRIPT_DIR), "pylib")
sys.path.insert(0, PYLIB_DIR)

from Splice_graph import Splice_graph
from Pretty_alignment import Pretty_alignment
from Bam_alignment_extractor import Bam_alignment_extractor
from Transcript import GTF_contig_to_transcripts
import LRAA_Globals
import Util_funcs

# Setup logging
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s - %(levelname)s - %(name)s - %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S",
)
logger = logging.getLogger(__name__)


def parse_args():
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(
        description="Correct BAM alignments using splice graph-guided realignment",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Basic usage
  python correct_bam_alignments.py --bam input.bam --gtf ref.gtf --genome genome.fa --output_prefix corrected
  
  # Skip normalization for small/targeted data
  python correct_bam_alignments.py --bam input.bam --gtf ref.gtf --genome genome.fa --output_prefix corrected --skip_normalization
  
  # Process specific contigs with multiple CPUs
  python correct_bam_alignments.py --bam input.bam --gtf ref.gtf --genome genome.fa --output_prefix corrected --contig chr1,chr2 --CPU 8
        """,
    )
    
    # Required arguments
    parser.add_argument(
        "--bam",
        type=str,
        required=True,
        help="Input BAM file (will be indexed if needed)",
    )
    parser.add_argument(
        "--gtf",
        type=str,
        required=True,
        help="Reference GTF file with transcript annotations",
    )
    parser.add_argument(
        "--genome",
        type=str,
        required=True,
        help="Genome FASTA file",
    )
    parser.add_argument(
        "--output_prefix",
        type=str,
        required=True,
        help="Prefix for output files",
    )
    
    # Optional arguments
    parser.add_argument(
        "--normalize_max_cov_level",
        type=int,
        default=1000,
        help="Coverage normalization level for splice graph construction (default: 1000)",
    )
    parser.add_argument(
        "--skip_normalization",
        action="store_true",
        help="Skip BAM normalization step (use for small/targeted datasets)",
    )
    parser.add_argument(
        "--CPU",
        type=int,
        default=4,
        help="Number of parallel processes (default: 4)",
    )
    parser.add_argument(
        "--contig",
        type=str,
        default=None,
        help="Comma-separated list of contigs to process (default: all)",
    )
    parser.add_argument(
        "--config_update",
        type=str,
        default=None,
        help="JSON file with config updates (LRAA config format)",
    )
    parser.add_argument(
        "--debug",
        action="store_true",
        help="Enable debug mode with verbose logging",
    )
    parser.add_argument(
        "--keep_temp",
        action="store_true",
        help="Keep temporary per-contig files (for debugging)",
    )
    
    return parser.parse_args()


def initialize_config(args):
    """Initialize LRAA configuration."""
    if args.debug:
        logging.getLogger().setLevel(logging.DEBUG)
        LRAA_Globals.DEBUG = True
    
    # Apply config updates if provided
    if args.config_update and os.path.exists(args.config_update):
        with open(args.config_update, "r") as f:
            config_updates = json.load(f)
            for key, value in config_updates.items():
                if key in LRAA_Globals.config:
                    LRAA_Globals.config[key] = value
                    logger.info(f"Config update: {key} = {value}")
    
    # Initialize splice graph parameters
    Splice_graph.init_sg_params({
        "read_aln_gap_merge_int": LRAA_Globals.config["read_aln_gap_merge_int"],
        "inter_exon_segment_merge_dist": 50,
        "max_genomic_contig_length": 1e10,
        "min_alt_splice_freq": LRAA_Globals.config["min_alt_splice_freq"],
        "min_alt_unspliced_freq": LRAA_Globals.config["min_alt_unspliced_freq"],
        "max_intron_length_for_exon_segment_filtering": 10000,
        "min_intron_support": 2,
        "min_terminal_splice_exon_anchor_length": LRAA_Globals.config[
            "min_terminal_splice_exon_anchor_length"
        ],
        "remove_unspliced_introns": False,
    })


def ensure_bam_indexed(bam_file):
    """Ensure BAM file is indexed, create index if missing."""
    index_file = bam_file + ".bai"
    if not os.path.exists(index_file):
        logger.info(f"Creating BAM index for {bam_file}")
        pysam.index(bam_file)
    return index_file


def get_fasta_reader(genome_file):
    """Get a pysam FastaFile reader."""
    if not os.path.exists(genome_file + ".fai"):
        logger.info(f"Creating FASTA index for {genome_file}")
        pysam.faidx(genome_file)
    return pysam.FastaFile(genome_file)


def get_contigs_from_bam(bam_file, restrict_contigs=None):
    """Get list of contigs to process from BAM file."""
    with pysam.AlignmentFile(bam_file, "rb") as bam:
        all_contigs = bam.references
    
    if restrict_contigs:
        # Filter to requested contigs
        requested = set(restrict_contigs.split(","))
        contigs = [c for c in all_contigs if c in requested]
        logger.info(f"Restricting to {len(contigs)} contigs: {contigs}")
    else:
        contigs = list(all_contigs)
        logger.info(f"Processing all {len(contigs)} contigs from BAM")
    
    return contigs


def normalize_bam_for_splice_graph(input_bam, normalize_max_cov_level, work_dir):
    """
    Normalize BAM file by strand for splice graph construction.
    
    Returns path to normalized BAM file.
    """
    if normalize_max_cov_level <= 0:
        return input_bam
    
    logger.info(f"Normalizing BAM to max coverage level: {normalize_max_cov_level}")
    
    work_dir = Path(work_dir)
    work_dir.mkdir(parents=True, exist_ok=True)
    
    input_path = Path(input_bam).resolve()
    base_name = input_path.stem
    
    normalized_bam = work_dir / f"{base_name}.norm_{normalize_max_cov_level}.bam"
    checkpoint = normalized_bam.parent / (normalized_bam.name + ".ok")
    
    if checkpoint.exists():
        logger.info(f"Using existing normalized BAM: {normalized_bam}")
        return str(normalized_bam)
    
    # Run normalization
    normalize_script = os.path.join(
        os.path.dirname(SCRIPT_DIR), "util", "normalize_bam_by_strand.py"
    )
    
    cmd = [
        normalize_script,
        "--input_bam", str(input_bam),
        "--output_bam", str(normalized_bam),
        "--max_cov", str(normalize_max_cov_level),
        "--workdir", str(work_dir / "norm_work"),
    ]
    
    logger.info(f"Running: {' '.join(cmd)}")
    subprocess.check_call(cmd)
    
    # Create checkpoint
    checkpoint.touch()
    
    # Index the normalized BAM
    ensure_bam_indexed(str(normalized_bam))
    
    logger.info(f"Normalized BAM created: {normalized_bam}")
    return str(normalized_bam)


def parse_gtf_for_contig_strand(gtf_file, contig_acc, contig_strand):
    """
    Parse GTF file for transcripts on a specific contig and strand.
    
    Returns list of Transcript objects.
    """
    logger.debug(f"Parsing GTF for {contig_acc}{contig_strand}")
    
    # Parse transcripts from GTF with contig/strand restriction
    contig_to_transcripts = GTF_contig_to_transcripts.parse_GTF_to_Transcripts(
        gtf_file,
        chr_restrict=contig_acc,
        strand_restrict=contig_strand,
        lend_restrict=None,
        rend_restrict=None,
    )
    
    # Get transcripts for this contig
    if contig_acc not in contig_to_transcripts:
        logger.debug(f"No transcripts found for {contig_acc}")
        return []
    
    transcripts = contig_to_transcripts[contig_acc]
    
    logger.debug(
        f"Loaded {len(transcripts)} transcripts for {contig_acc}{contig_strand}"
    )
    return transcripts


def build_splice_graph(
    contig_acc,
    contig_strand,
    contig_seq,
    bam_file,
    gtf_transcripts,
):
    """
    Build splice graph for a contig/strand.
    
    Returns Splice_graph object.
    """
    logger.info(f"Building splice graph for {contig_acc}{contig_strand}")
    
    sg = Splice_graph()
    sg.build_splice_graph_for_contig(
        contig_acc=contig_acc,
        contig_strand=contig_strand,
        contig_seq_str=contig_seq,
        alignments_bam_file=bam_file,
        region_lend=None,
        region_rend=None,
        input_transcripts=gtf_transcripts if gtf_transcripts else None,
        quant_mode=False,
        restrict_splice_type=None,
        SE_read_encapsulation_mask=None,
    )
    
    if sg.is_empty():
        logger.warning(f"Splice graph is empty for {contig_acc}{contig_strand}")
        return None
    
    logger.info(f"Splice graph built successfully for {contig_acc}{contig_strand}")
    return sg


def process_contig_strand(args_tuple):
    """
    Process a single contig/strand combination.
    
    This is the main worker function for parallel processing.
    Returns tuple: (contig, strand, stats_dict, temp_files_dict, error_msg)
    """
    (
        contig_acc,
        contig_strand,
        bam_file,
        normalized_bam,
        gtf_file,
        genome_file,
        output_prefix,
        work_dir,
    ) = args_tuple
    
    logger.info(f"Processing {contig_acc}{contig_strand}")
    
    stats = {
        "contig": contig_acc,
        "strand": contig_strand,
        "total_reads": 0,
        "corrected_reads": 0,
        "left_corrections": 0,
        "right_corrections": 0,
        "both_corrections": 0,
    }
    
    try:
        # Get contig sequence
        fasta = get_fasta_reader(genome_file)
        contig_seq = fasta.fetch(contig_acc)
        
        # Parse GTF for this contig/strand
        gtf_transcripts = parse_gtf_for_contig_strand(
            gtf_file, contig_acc, contig_strand
        )
        
        # Build splice graph from normalized BAM
        splice_graph = build_splice_graph(
            contig_acc,
            contig_strand,
            contig_seq,
            normalized_bam,
            gtf_transcripts,
        )
        
        if splice_graph is None:
            # This is acceptable - no reads or no introns on this contig/strand
            logger.info(
                f"Skipping {contig_acc}{contig_strand}: no reads or no multi-exon reads (no splice graph needed)"
            )
            return (contig_acc, contig_strand, stats, None, None)
        
        # Process alignments from original BAM
        result = process_alignments_for_contig_strand(
            contig_acc,
            contig_strand,
            contig_seq,
            bam_file,
            splice_graph,
            output_prefix,
            work_dir,
            stats,
        )
        
        return (contig_acc, contig_strand, stats, result, None)
        
    except Exception as e:
        error_msg = f"{type(e).__name__}: {str(e)}"
        logger.error(
            f"Error processing {contig_acc}{contig_strand}: {error_msg}",
            exc_info=True
        )
        return (contig_acc, contig_strand, stats, None, error_msg)


def process_alignments_for_contig_strand(
    contig_acc,
    contig_strand,
    contig_seq,
    bam_file,
    splice_graph,
    output_prefix,
    work_dir,
    stats,
):
    """
    Process alignments for a contig/strand, correct them, and write temp BAM files.
    
    Returns dict with paths to temporary BAM files.
    """
    work_dir = Path(work_dir)
    work_dir.mkdir(parents=True, exist_ok=True)
    
    # Temporary output files for this contig/strand
    strand_token = contig_strand if contig_strand else "both"
    temp_full = work_dir / f"{output_prefix}.{contig_acc}.{strand_token}.full.bam"
    temp_faulty = work_dir / f"{output_prefix}.{contig_acc}.{strand_token}.faulty.bam"
    temp_corrected = work_dir / f"{output_prefix}.{contig_acc}.{strand_token}.corrected.bam"
    temp_annotated = work_dir / f"{output_prefix}.{contig_acc}.{strand_token}.annotated.bam"
    
    # Open input BAM
    input_bam = pysam.AlignmentFile(bam_file, "rb")
    
    # Open output BAMs
    full_bam = pysam.AlignmentFile(str(temp_full), "wb", template=input_bam)
    faulty_bam = pysam.AlignmentFile(str(temp_faulty), "wb", template=input_bam)
    corrected_bam = pysam.AlignmentFile(str(temp_corrected), "wb", template=input_bam)
    annotated_bam = pysam.AlignmentFile(str(temp_annotated), "wb", template=input_bam)
    
    # Fetch alignments for this contig
    try:
        alignments = input_bam.fetch(contig_acc)
    except Exception as e:
        logger.error(f"Error fetching alignments for {contig_acc}: {e}")
        input_bam.close()
        full_bam.close()
        faulty_bam.close()
        corrected_bam.close()
        return None
    
    # Process alignments
    alignments_to_correct = []
    alignment_map = {}  # read_name -> original pysam record
    
    for aln in alignments:
        # Skip unmapped, secondary, and supplementary alignments
        if aln.is_unmapped:
            full_bam.write(aln)
            annotated_bam.write(aln)
            continue
        if aln.is_secondary or aln.is_supplementary:
            full_bam.write(aln)
            annotated_bam.write(aln)
            continue
        
        # Check strand
        aln_strand = "-" if aln.is_reverse else "+"
        if contig_strand is not None and aln_strand != contig_strand:
            full_bam.write(aln)
            annotated_bam.write(aln)
            continue
        
        stats["total_reads"] += 1
        
        # Convert to Pretty_alignment
        try:
            pretty_aln = Pretty_alignment.get_pretty_alignment(aln)
            
            # Check if candidate for correction
            if pretty_aln.is_softclip_realign_candidate():
                alignments_to_correct.append(pretty_aln)
                alignment_map[pretty_aln.get_read_name()] = aln
        except Exception as e:
            logger.warning(f"Failed to create Pretty_alignment for {aln.query_name}: {e}")
            full_bam.write(aln)
            annotated_bam.write(aln)
            continue
    
    # Perform correction on candidates
    if alignments_to_correct:
        logger.info(
            f"Correcting {len(alignments_to_correct)} candidate alignments "
            f"for {contig_acc}{contig_strand}"
        )
        
        # Store original segments for comparison
        original_segments = {}
        for pa in alignments_to_correct:
            original_segments[pa.get_read_name()] = [
                list(seg) for seg in pa.get_pretty_alignment_segments()
            ]
        
        # Run correction
        Pretty_alignment.try_correct_alignments(
            alignments_to_correct, splice_graph, contig_seq
        )
        
        # Process corrected alignments
        for pretty_aln in alignments_to_correct:
            read_name = pretty_aln.get_read_name()
            orig_aln = alignment_map[read_name]
            orig_segs = original_segments[read_name]
            curr_segs = pretty_aln.get_pretty_alignment_segments()
            
            # Check if alignment was actually corrected
            was_corrected = orig_segs != curr_segs
            
            if was_corrected:
                stats["corrected_reads"] += 1
                
                # Determine correction type
                correction_type = determine_correction_type(orig_segs, curr_segs)
                if "left" in correction_type:
                    stats["left_corrections"] += 1
                if "right" in correction_type:
                    stats["right_corrections"] += 1
                if correction_type == "both":
                    stats["both_corrections"] += 1
                
                # Create corrected pysam record
                try:
                    corrected_aln = pretty_aln.to_corrected_pysam_alignment(orig_aln)
                    
                    # Validate CIGAR before writing
                    read_len = len(corrected_aln.query_sequence)
                    cigar_query_len = sum(
                        length for op, length in corrected_aln.cigartuples
                        if op in (0, 1, 4, 7, 8)  # M, I, S, =, X consume query
                    )
                    
                    if cigar_query_len != read_len:
                        logger.warning(
                            f"CIGAR validation failed for {read_name}: "
                            f"CIGAR query length {cigar_query_len} != read length {read_len}. "
                            f"Using original alignment instead."
                        )
                        # Write original instead of corrupted correction
                        full_bam.write(orig_aln)
                        
                        # Add XC tag with ~ delimiters for faulty read
                        import copy
                        faulty_annotated = copy.deepcopy(orig_aln)
                        faulty_annotated.set_tag("XC", f"~{correction_type}~", value_type="Z")
                        annotated_bam.write(faulty_annotated)
                        
                        stats["corrected_reads"] -= 1  # Don't count as corrected
                        if "left" in correction_type:
                            stats["left_corrections"] -= 1
                        if "right" in correction_type:
                            stats["right_corrections"] -= 1
                        if correction_type == "both":
                            stats["both_corrections"] -= 1
                        continue
                    
                    # Add tags to mark correction
                    corrected_aln.set_tag("XC", correction_type, value_type="Z")
                    corrected_aln.set_tag("OC", orig_aln.cigarstring, value_type="Z")
                    corrected_aln.set_tag("OA", orig_aln.reference_start, value_type="i")
                    
                    # Create annotated original with XC tag using ~ delimiters
                    import copy
                    annotated_orig = copy.deepcopy(orig_aln)
                    annotated_orig.set_tag("XC", f"~{correction_type}~", value_type="Z")
                    
                    # Write to output files
                    full_bam.write(corrected_aln)
                    
                    # Write faulty read with same XC tag (~ delimiters)
                    faulty_with_tag = copy.deepcopy(orig_aln)
                    faulty_with_tag.set_tag("XC", f"~{correction_type}~", value_type="Z")
                    faulty_bam.write(faulty_with_tag)
                    
                    corrected_bam.write(corrected_aln)
                    annotated_bam.write(annotated_orig)
                    
                except Exception as e:
                    logger.warning(
                        f"Failed to convert corrected alignment for {read_name}: {e}"
                    )
                    # Fall back to original
                    full_bam.write(orig_aln)
                    
                    # Add XC tag with ~ delimiters for faulty read
                    import copy
                    faulty_annotated = copy.deepcopy(orig_aln)
                    faulty_annotated.set_tag("XC", f"~{correction_type}~", value_type="Z")
                    annotated_bam.write(faulty_annotated)
                    
                    stats["corrected_reads"] -= 1  # Don't count as corrected
                    if "left" in correction_type:
                        stats["left_corrections"] -= 1
                    if "right" in correction_type:
                        stats["right_corrections"] -= 1
                    if correction_type == "both":
                        stats["both_corrections"] -= 1
            else:
                # Was a candidate but correction didn't change anything (faulty)
                # Tag it with ~ delimiters to indicate it's a failed correction
                potential_correction_type = determine_potential_correction_type(orig_aln)
                
                full_bam.write(orig_aln)
                
                if potential_correction_type != "none":
                    # Add XC tag with ~ delimiters for faulty candidate
                    import copy
                    faulty_annotated = copy.deepcopy(orig_aln)
                    faulty_annotated.set_tag("XC", f"~{potential_correction_type}~", value_type="Z")
                    annotated_bam.write(faulty_annotated)
                else:
                    # No soft clips, shouldn't have been a candidate
                    annotated_bam.write(orig_aln)
    
    # Close files
    input_bam.close()
    full_bam.close()
    faulty_bam.close()
    corrected_bam.close()
    annotated_bam.close()
    
    # Sort and index temp BAMs (corrections can change positions, causing unsorted output)
    logger.debug(f"Sorting and indexing temp BAMs for {contig_acc}{contig_strand}")
    
    for bam_path in [temp_full, temp_faulty, temp_corrected, temp_annotated]:
        if os.path.exists(bam_path) and os.path.getsize(bam_path) > 0:
            sorted_bam = str(bam_path) + ".sorted.bam"
            try:
                pysam.sort("-o", sorted_bam, str(bam_path))
                # Replace unsorted with sorted
                os.rename(sorted_bam, str(bam_path))
                # Index the sorted BAM
                pysam.index(str(bam_path))
            except Exception as e:
                logger.warning(f"Failed to sort/index {bam_path}: {e}")
    
    logger.info(
        f"Completed {contig_acc}{contig_strand}: "
        f"{stats['corrected_reads']}/{stats['total_reads']} corrected"
    )
    
    return {
        "full": str(temp_full),
        "faulty": str(temp_faulty),
        "corrected": str(temp_corrected),
        "annotated": str(temp_annotated),
    }


def determine_correction_type(orig_segments, corrected_segments):
    """
    Determine what type of correction was applied.
    
    Returns: "left", "right", or "both"
    """
    if not orig_segments or not corrected_segments:
        return "unknown"
    
    left_changed = orig_segments[0] != corrected_segments[0]
    right_changed = orig_segments[-1] != corrected_segments[-1]
    
    if left_changed and right_changed:
        return "both"
    elif left_changed:
        return "left"
    elif right_changed:
        return "right"
    else:
        return "none"


def determine_potential_correction_type(pysam_aln):
    """
    Determine which ends have soft clips that could potentially be corrected.
    Used for tagging faulty candidates.
    
    Returns: "left", "right", "both", or "none"
    """
    if not pysam_aln.cigartuples:
        return "none"
    
    # Check for soft clips at ends (operation 4 is soft clip)
    left_softclip = pysam_aln.cigartuples[0][0] == 4
    right_softclip = pysam_aln.cigartuples[-1][0] == 4
    
    if left_softclip and right_softclip:
        return "both"
    elif left_softclip:
        return "left"
    elif right_softclip:
        return "right"
    else:
        return "none"


def merge_bam_files(temp_bam_list, output_bam, header_template):
    """
    Merge temporary BAM files into final output.
    
    Uses samtools merge for efficiency.
    """
    if not temp_bam_list:
        logger.warning(f"No temp files to merge for {output_bam}")
        return
    
    # Filter to existing files
    existing_files = [f for f in temp_bam_list if os.path.exists(f) and os.path.getsize(f) > 0]
    
    if not existing_files:
        logger.warning(f"No valid temp files to merge for {output_bam}")
        # Create empty BAM
        with pysam.AlignmentFile(header_template, "rb") as template:
            pysam.AlignmentFile(output_bam, "wb", template=template).close()
        return
    
    if len(existing_files) == 1:
        # Just copy the single file
        shutil.copy(existing_files[0], output_bam)
    else:
        # Use samtools merge
        cmd = ["samtools", "merge", "-f", output_bam] + existing_files
        logger.info(f"Merging {len(existing_files)} BAM files into {output_bam}")
        subprocess.check_call(cmd)
    
    # Index output
    pysam.index(output_bam)
    logger.info(f"Created and indexed: {output_bam}")


def write_stats_report(all_stats, output_file):
    """Write correction statistics to TSV file."""
    with open(output_file, "w") as f:
        # Header
        f.write(
            "contig\tstrand\ttotal_reads\tcorrected_reads\tpct_corrected\t"
            "left_corrections\tright_corrections\tboth_corrections\n"
        )
        
        # Data rows
        for stats in all_stats:
            pct = (
                100.0 * stats["corrected_reads"] / stats["total_reads"]
                if stats["total_reads"] > 0
                else 0.0
            )
            f.write(
                f"{stats['contig']}\t{stats['strand']}\t{stats['total_reads']}\t"
                f"{stats['corrected_reads']}\t{pct:.2f}\t"
                f"{stats['left_corrections']}\t{stats['right_corrections']}\t"
                f"{stats['both_corrections']}\n"
            )
        
        # Summary
        total_reads = sum(s["total_reads"] for s in all_stats)
        total_corrected = sum(s["corrected_reads"] for s in all_stats)
        total_left = sum(s["left_corrections"] for s in all_stats)
        total_right = sum(s["right_corrections"] for s in all_stats)
        total_both = sum(s["both_corrections"] for s in all_stats)
        total_pct = (
            100.0 * total_corrected / total_reads if total_reads > 0 else 0.0
        )
        
        f.write(
            f"TOTAL\tALL\t{total_reads}\t{total_corrected}\t{total_pct:.2f}\t"
            f"{total_left}\t{total_right}\t{total_both}\n"
        )
    
    logger.info(f"Wrote statistics to {output_file}")


def write_failure_report(failed_contigs, output_file):
    """Write failure details to TSV file."""
    with open(output_file, "w") as f:
        # Header
        f.write("contig\tstrand\treason\n")
        
        # Data rows
        for contig, strand, error_msg in failed_contigs:
            # Clean up error message for TSV (replace tabs/newlines)
            clean_msg = error_msg.replace("\t", " ").replace("\n", " ") if error_msg else "Unknown error"
            f.write(f"{contig}\t{strand}\t{clean_msg}\n")
    
    logger.info(f"Wrote failure report to {output_file}")


def main():
    """Main entry point."""
    args = parse_args()
    
    # Track if we had any errors
    had_errors = False
    
    logger.info("=" * 80)
    logger.info("BAM Alignment Correction Utility")
    logger.info("=" * 80)
    
    # Initialize configuration
    initialize_config(args)
    
    # Validate inputs
    if not os.path.exists(args.bam):
        sys.exit(f"Error: BAM file not found: {args.bam}")
    if not os.path.exists(args.gtf):
        sys.exit(f"Error: GTF file not found: {args.gtf}")
    if not os.path.exists(args.genome):
        sys.exit(f"Error: Genome file not found: {args.genome}")
    
    # Ensure BAM is indexed
    ensure_bam_indexed(args.bam)
    
    # Create work directory
    work_dir = Path(f"__{args.output_prefix}_correction_work")
    work_dir.mkdir(parents=True, exist_ok=True)
    logger.info(f"Working directory: {work_dir}")
    
    # Normalize BAM if requested
    if args.skip_normalization:
        logger.info("Skipping normalization (--skip_normalization)")
        normalized_bam = args.bam
    else:
        normalized_bam = normalize_bam_for_splice_graph(
            args.bam,
            args.normalize_max_cov_level,
            work_dir / "normalization",
        )
    
    # Get contigs to process
    contigs = get_contigs_from_bam(args.bam, args.contig)
    
    # Prepare work items for parallel processing
    work_items = []
    for contig in contigs:
        for strand in ["+", "-"]:
            work_items.append(
                (
                    contig,
                    strand,
                    args.bam,
                    normalized_bam,
                    args.gtf,
                    args.genome,
                    args.output_prefix,
                    str(work_dir / "temp_bams"),
                )
            )
    
    logger.info(f"Processing {len(work_items)} contig/strand combinations")
    
    # Process in parallel
    all_stats = []
    temp_files = {"full": [], "faulty": [], "corrected": [], "annotated": []}
    
    if args.CPU > 1 and len(work_items) > 1:
        logger.info(f"Using {args.CPU} parallel processes")
        with mp.Pool(processes=args.CPU) as pool:
            results = pool.map(process_contig_strand, work_items)
    else:
        logger.info("Processing sequentially (single CPU)")
        results = [process_contig_strand(item) for item in work_items]
    
    # Collect results and track actual errors (not just empty contigs)
    failed_contigs = []
    for contig, strand, stats, temp_bam_dict, error_msg in results:
        all_stats.append(stats)
        if temp_bam_dict:
            temp_files["full"].append(temp_bam_dict["full"])
            temp_files["faulty"].append(temp_bam_dict["faulty"])
            temp_files["corrected"].append(temp_bam_dict["corrected"])
            temp_files["annotated"].append(temp_bam_dict["annotated"])
        elif error_msg is not None:
            # Only track actual errors, not "no data" situations
            failed_contigs.append((contig, strand, error_msg))
            had_errors = True
    
    # Merge temporary files into final outputs
    logger.info("Merging temporary BAM files into final outputs")
    
    final_full = f"{args.output_prefix}.corrected.bam"
    final_faulty = f"{args.output_prefix}.faulty_only.bam"
    final_corrected = f"{args.output_prefix}.corrected_only.bam"
    final_annotated = f"{args.output_prefix}.annotated_originals.bam"
    
    merge_bam_files(temp_files["full"], final_full, args.bam)
    merge_bam_files(temp_files["faulty"], final_faulty, args.bam)
    merge_bam_files(temp_files["corrected"], final_corrected, args.bam)
    merge_bam_files(temp_files["annotated"], final_annotated, args.bam)
    
    # Write statistics
    stats_file = f"{args.output_prefix}.correction_stats.tsv"
    write_stats_report(all_stats, stats_file)
    
    # Write failure report if there were errors
    if failed_contigs:
        failure_file = f"{args.output_prefix}.failed_contigs.tsv"
        write_failure_report(failed_contigs, failure_file)
        logger.warning(f"Failure details written to: {failure_file}")
    
    # Cleanup
    if not args.keep_temp:
        logger.info("Cleaning up temporary files")
        shutil.rmtree(work_dir)
    else:
        logger.info(f"Temporary files retained in: {work_dir}")
    
    # Summary
    logger.info("=" * 80)
    if had_errors:
        logger.error("Correction completed with ERRORS!")
        logger.error(f"{len(failed_contigs)} contig/strand combinations failed")
        logger.error(f"See {args.output_prefix}.failed_contigs.tsv for details")
    else:
        logger.info("Correction complete - SUCCESS!")
        total_processed = sum(s["total_reads"] for s in all_stats)
        total_corrected = sum(s["corrected_reads"] for s in all_stats)
        if total_processed > 0:
            logger.info(f"Processed {total_processed} reads, corrected {total_corrected} ({100*total_corrected/total_processed:.2f}%)")
        else:
            logger.info("No reads processed")
    logger.info(f"Full corrected BAM: {final_full}")
    logger.info(f"Faulty alignments only: {final_faulty}")
    logger.info(f"Corrected alignments only: {final_corrected}")
    logger.info(f"Annotated originals BAM: {final_annotated}")
    logger.info(f"Statistics report: {stats_file}")
    if failed_contigs:
        logger.info(f"Failed contigs report: {failure_file}")
    logger.info("=" * 80)
    
    # Exit with proper status code
    if had_errors:
        sys.exit(1)


if __name__ == "__main__":
    main()
