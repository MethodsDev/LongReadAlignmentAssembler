#!/usr/bin/env python3

import sys, os, re
import subprocess
import argparse
import logging
import pysam
import random
sys.path.insert(0, os.path.sep.join([os.path.dirname(os.path.realpath(__file__)), "../pylib"]))
from Pipeliner import Pipeliner, Command

logging.basicConfig(level=logging.INFO, 
                    format='%(asctime)s : %(levelname)s : %(message)s',
                    datefmt='%H:%M:%S')
logger = logging.getLogger(__name__)




def main():

    parser = argparse.ArgumentParser(description="normalize bam in a strand-specific manner", formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument("--input_bam", type=str, required=True, help="input bam filename")
    parser.add_argument("--output_bam", type=str, required=True, help="output for normalied bam file")
    parser.add_argument("--normalize_max_cov_level", type=int, default=1000, help="normalize to max read coverage level before assembly (default: 1000)")
    parser.add_argument("--read_start_bin_size", type=int, default=100, help="group alignments by start positions every --read_start_bin_size number of bases along the genome")
    parser.add_argument("--random_seed", type=int, default=42, help="random seed for reproducible sampling (default: 42)")

    args = parser.parse_args()
    
    input_bam_filename = args.input_bam
    output_bam_filename = args.output_bam
    normalize_max_cov_level = args.normalize_max_cov_level
    read_start_bin_size = args.read_start_bin_size
    random_seed = args.random_seed

    # Set random seed for reproducible sampling
    random.seed(random_seed)
    logger.info(f"Using random seed: {random_seed}")

    
    pipeliner = Pipeliner("__chckpts")
    
    # first separate input bam into strand-specific files
    SS_output_prefix = os.path.basename(input_bam_filename) + ".SS"

    scriptdir = os.path.abspath(os.path.dirname(__file__))
    cmd = " ".join([os.path.join(scriptdir, "separate_bam_by_strand.py"),
                    "--bam {}".format(input_bam_filename),
                    "--output_prefix {}".format(SS_output_prefix)])
    
    
    pipeliner.add_commands([Command(cmd, "sep_by_strand.ok")])

    pipeliner.run()
    
    ## run normalizations
    bamsifter_prog = os.path.join(os.path.abspath(os.path.dirname(__file__)), "../plugins/bamsifter/bamsifter")
    SS_bam_files = [SS_output_prefix + x for x in (".+.bam", ".-.bam") ]

    SS_norm_bam_files = list()
    
    for SS_bam_file in SS_bam_files:

        norm_bam_filename = f"{SS_bam_file}.{normalize_max_cov_level}.bam"

        norm_bam_checkpoint = norm_bam_filename + ".ok"

        if not os.path.exists(norm_bam_checkpoint):
            sift_bam(SS_bam_file, norm_bam_filename, normalize_max_cov_level, read_start_bin_size)
            subprocess.check_call("touch {}".format(norm_bam_checkpoint), shell=True)
        
        SS_norm_bam_files.append(norm_bam_filename)
        

    # merge the norm SS bam filenames into the final output file

    cmd = f"samtools merge -f {output_bam_filename} " + " ".join(SS_norm_bam_files)
    pipeliner.add_commands([Command(cmd, "SS_merge.ok")])


    cmd = f"samtools index {output_bam_filename}"
    pipeliner.add_commands([Command(cmd, "index_merged.ok")])

    pipeliner.run()

    logger.info("Done.  See SS-normalized bam: {}".format(output_bam_filename))

    sys.exit(0)


def sift_bam(SS_bam_file, norm_bam_filename, normalize_max_cov_level, read_start_bin_size):

    bamfile_reader = pysam.AlignmentFile(SS_bam_file, "rb")

    bam_writer = pysam.AlignmentFile(norm_bam_filename, "wb", template=bamfile_reader)

    prev_chrom = -1
    prev_read_bin = -1

    read_queue = list()


    def write_reads(read_list):

        ## do sampling if needed
        if len(read_list) > normalize_max_cov_level:
            rand_indices = sorted(random.sample(range(len(read_list)), normalize_max_cov_level))
            read_list = [read_list[i] for i in rand_indices]
        
        for read in read_list:
            bam_writer.write(read)
        
        return

    
    for read in bamfile_reader.fetch():

        start_pos = read.reference_start + 1
        read_bin = int(start_pos / read_start_bin_size)
        chrom = read.reference_id

                    
        if prev_chrom != chrom or read_bin != prev_read_bin:

            if len(read_queue) > 0:
                write_reads(read_queue)
                read_queue.clear()

        read_queue.append(read)
        
        prev_chrom = chrom
        prev_read_bin = read_bin

    if len(read_queue) > 0:
        write_reads(read_queue)

    bam_writer.close()
    bamfile_reader.close()
    
    return
    
        
        

    
if __name__=='__main__':
    main()
