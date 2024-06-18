#!/usr/bin/env python3

import sys, os, re
import Transcript
import MultiPath
import MultiPathCounter
import Simple_path_utils as SPU
from collections import defaultdict
import LRAA_Globals
from LRAA_Globals import SPACER, DEBUG
import logging
from math import log

logger = logging.getLogger(__name__)



class Quantify:

    def __init__(self):

        self._path_node_id_to_gene_ids = defaultdict(set)

        self._gene_id_to_transcript_objs = defaultdict(set)

        self._read_name_to_multipath = dict()
        
        return


    def quantify(self, splice_graph, transcripts, mp_counter, run_EM=True):

        assert type(transcripts) == list
        assert type(transcripts[0]) == Transcript.Transcript
        assert type(mp_counter) == MultiPathCounter.MultiPathCounter

        
        # assign path nodes to gene
        # also assign gene_id to transcript objs
        self._assign_path_nodes_to_gene(transcripts)

        self._assign_reads_to_transcripts(splice_graph, mp_counter)

        transcript_to_fractional_read_assignment = self._estimate_isoform_read_support(transcripts, run_EM)

        # see documentation for _estimate_isoform_read_support() below
        
        return transcript_to_fractional_read_assignment

        
    def _assign_path_nodes_to_gene(self, transcripts):

        for transcript in transcripts:
            
            simplepath = transcript._simplepath

            if simplepath is None:
                logger.warn("simplepath is not avaialble for transcript: {}".format(transcript))
                continue
            
            #assert simplepath is not None, "Error, simplepath not set for transcript obj: {}".format(transcript)

            transcript_id = transcript.get_transcript_id()
            gene_id = transcript.get_gene_id()
            self._gene_id_to_transcript_objs[gene_id].add(transcript)
            
            for node_id in simplepath:
                if node_id != SPACER:
                    self._path_node_id_to_gene_ids[node_id].add(gene_id)

        return



    def _assign_reads_to_transcripts(self, splice_graph, mp_counter, fraction_read_align_overlap=LRAA_Globals.config['fraction_read_align_overlap']):

        logger.info("# Assigning reads to transcripts")
        
        # assign to gene based on majority voting of nodes.
        # TODO:// might want or need this to involve length and/or feature type weighted shared node voting

        mp_count_pairs = mp_counter.get_all_MultiPathCountPairs()

        gene_unanchored_mp_count_pairs = list()


        num_paths_total = 0
        num_read_counts_total = 0

        num_paths_anchored_to_gene = 0
        num_read_counts_anchored_to_gene = 0
        
        num_paths_assigned = 0
        num_read_counts_assigned = 0
        
        
        for mp_count_pair in mp_count_pairs:
            mp, count = mp_count_pair.get_multipath_and_count()

            num_paths_total += 1
            num_read_counts_total += count
            
            sp = mp.get_simple_path()

            top_gene = self._get_gene_with_best_node_matches_to_simplepath(sp)

            if top_gene is None:
                gene_unanchored_mp_count_pairs.append(mp_count_pair)
                logger.debug("mp_count_pair unanchored: " + str(mp_count_pair))
                
                continue

            logger.debug("mp_count_pair {} anchored to gene: {}".format(mp_count_pair, top_gene))
            num_paths_anchored_to_gene += 1
            num_read_counts_anchored_to_gene += count
            
            
            ## assign reads to transcripts
            gene_isoforms = self._gene_id_to_transcript_objs[top_gene]
            transcripts_assigned = self._assign_path_to_transcript(splice_graph, mp, gene_isoforms, fraction_read_align_overlap,
                                                                   trim_TSS_polyA = False, test_exact=True)

            if transcripts_assigned is None:
                # try again with TSS and polyA trimming
                transcripts_assigned = self._assign_path_to_transcript(splice_graph, mp, gene_isoforms, fraction_read_align_overlap,
                                                                       trim_TSS_polyA = True, test_exact=True)          
            

            if transcripts_assigned is None:
                # try again with TSS and PolyA trimmed off
                transcripts_assigned = self._assign_path_to_transcript(splice_graph, mp, gene_isoforms, fraction_read_align_overlap,
                                                                       trim_TSS_polyA = True, test_exact=False)

            if transcripts_assigned is None:
                # last resort, do majority voting
                transcripts_assigned = self._assign_path_to_transcript_by_majority_voting(splice_graph, mp, gene_isoforms)
                

                
            if transcripts_assigned is None:
                logger.debug("mp_count_pair {} maps to gene but no isoform(transcript)".format(mp_count_pair))
            else:
                logger.debug("mp_count_pair {} maps to transcripts: {}".format(mp_count_pair, transcripts_assigned))
                for transcript in transcripts_assigned:
                    transcript.add_read_names(mp.get_read_names())
                    
                num_paths_assigned += 1
                num_read_counts_assigned += count

        if num_paths_total == 0:
            num_paths_total = 1e-5 # make nonzero to avoid div-by-zero below

        if num_read_counts_total == 0:
            num_read_counts_total = 1e-5 # ditto above
            
        ## audit summary
        audit_txt = "\n".join(
            [
             "num_paths_total: {}, num_read_counts_total: {}".format(num_paths_total, num_read_counts_total),
             "\tnum_paths_anchored_to_gene: {} = {:.2f}%, num_read_counts_anchored_to_gene: {} = {:.2f}%\n".format(num_paths_anchored_to_gene,
                                                                                                                 num_paths_anchored_to_gene/num_paths_total*100,
                                                                                                                 num_read_counts_anchored_to_gene,
                                                                                                                 num_read_counts_anchored_to_gene/num_read_counts_total*100),
            "\tnum_paths_assigned_to_trans: {} = {:.2f}%, num_read_counts_assigned_to_trans: {} = {:.2f}%\n".format(num_paths_assigned,
                                                                                                                  num_paths_assigned/num_paths_total*100,
                                                                                                                  num_read_counts_assigned,
                                                                                                                  num_read_counts_assigned/num_read_counts_total*100)])
        
                
        logger.info(audit_txt)
                
        return


    def _get_gene_with_best_node_matches_to_simplepath(self, simplepath):

        gene_ranker = defaultdict(int)

        for node in simplepath:
            if node != SPACER:
                if node in self._path_node_id_to_gene_ids:
                    gene_set = self._path_node_id_to_gene_ids[node]
                    for gene_id in gene_set:
                        gene_ranker[gene_id] += 1

        if len(gene_ranker) == 0:
            return None
        else:
            genes_ranked = sorted(gene_ranker.keys(), key=lambda x: gene_ranker[x], reverse=True)
            return genes_ranked[0]
        

        
    def _assign_path_to_transcript(self, splice_graph, mp, transcripts, fraction_read_align_overlap,
                                   trim_TSS_polyA = False, test_exact=True):
        
        assert type(mp) == MultiPath.MultiPath
        assert type(transcripts) == set, "Error, type(transcripts) is {} not set ".format(type(transcripts))
        assert type(list(transcripts)[0]) == Transcript.Transcript
        assert fraction_read_align_overlap >= 0 and fraction_read_align_overlap <= 1.0, "Error, fraction_read_align_overlap must be between 0 and 1.0"

        contig_strand = splice_graph.get_contig_strand()
        
        read_sp = mp.get_simple_path()
        if trim_TSS_polyA:
            read_sp, read_TSS_id, read_polyA_id = SPU.trim_TSS_and_PolyA(read_sp, contig_strand)

        # store read name to mp for later debugging.
        for read_name in mp.get_read_names():
            self._read_name_to_multipath[read_name] = mp

        
        transcripts_compatible_with_read = list()
        
        for transcript in transcripts:
            transcript_sp = transcript._simplepath
                        
            assert transcript_sp is not None

            if trim_TSS_polyA:
                transcript_sp, transcript_TSS_id, transcript_polyA_id = SPU.trim_TSS_and_PolyA(transcript_sp, contig_strand)
                

            if test_exact:

                if transcript_sp == read_sp:
                    logger.debug("[trim_TSS_polyA={} test_exact={}]  Read {} IDENTICAL with transcript {}".format(trim_TSS_polyA, test_exact, read_sp, transcript_sp))
                    transcripts_compatible_with_read.append(transcript)
                    

            else:

                if (SPU.are_overlapping_and_compatible_NO_gaps_in_overlap(transcript_sp, read_sp)
                    and
                    SPU.fraction_read_overlap(splice_graph, read_sp, transcript_sp) >= fraction_read_align_overlap):

                    logger.debug("[trim_TSS_polyA={} test_exact={}]  Read {} COMPATIBLE with transcript {}".format(trim_TSS_polyA, test_exact, read_sp, transcript_sp))
                    #print("Read {} compatible with transcript {}".format(read_sp, transcript_sp))
                    transcripts_compatible_with_read.append(transcript)

                else:
                    logger.debug("[trim_TSS_polyA={} test_exact={}]  Read {} NOT_compatible with transcript {}".format(trim_TSS_polyA, test_exact,read_sp, transcript_sp))
                        

        if len(transcripts_compatible_with_read) == 0:
            return None
        else:
            return transcripts_compatible_with_read



        
    def _assign_path_to_transcript_by_majority_voting(self, splice_graph, mp, transcripts):
        
        assert type(mp) == MultiPath.MultiPath
        assert type(transcripts) == set, "Error, type(transcripts) is {} not set ".format(type(transcripts))
        assert type(list(transcripts)[0]) == Transcript.Transcript

        contig_strand = splice_graph.get_contig_strand()
        
        read_sp = mp.get_simple_path()
 
        # store read name to mp for later debugging.
        for read_name in mp.get_read_names():
            self._read_name_to_multipath[read_name] = mp

        
        scored_transcripts = list()
        
        for transcript in transcripts:
            transcript_sp = transcript._simplepath

            shared_simple_nodes = [simple_node for simple_node in read_sp if simple_node in transcript_sp]

            if len(shared_simple_nodes) > 0:
                scored_transcripts.append([len(shared_simple_nodes), transcript])

        if len(scored_transcripts) > 0:
            scored_transcripts = sorted(scored_transcripts, key=lambda x: x[0], reverse=True)
            top_transcript_score_pair = scored_transcripts.pop(0)
            top_transcript_score, top_transcript = top_transcript_score_pair
            top_transcripts = [top_transcript]
            # capture ties
            for alt_top_transcript in scored_transcripts:
                if alt_top_transcript[0] == top_transcript_score:
                    top_transcripts.append(alt_top_transcript[1])
            
            return top_transcripts

        else:
            return None
                

    def _estimate_isoform_read_support(self, transcripts, run_EM):

        """

        Given the reads assigned to the transcript (accessed with transcript.get_read_names() ) 
            Read counts are assigned to transcripts taking into account multiple read mappings
            Without EM (run_EM is False), read counts are equally divided among the isoforms they're assigned.
            With EM, they're assigned fractionally according to inferred isoform expression levels in EM cycles.

            Final read counts and isoform fraction values are stored in the transcript objects themselves, and accessed as:
                transcript.get_read_counts_assigned() and transcript.get_isoform_fraction()
                
        returns transcript_to_fractional_read_assignment
             with structure [transcript_id][read_name] = frac_read_assigned

        """

        
        read_name_to_transcripts = defaultdict(set)

        transcript_to_fractional_read_assignment = defaultdict(dict)
        
        for transcript in transcripts:
            read_names = transcript.get_read_names()
            for read_name in read_names:
                read_name_to_transcripts[read_name].add(transcript)
                
        num_mapped_reads = len(read_name_to_transcripts)

        transcript_to_read_count = defaultdict(float)

        transcript_to_expr_val = defaultdict(float)
        
        ## first round of EM for now - split evenly across mapped transcripts.
        for transcript in transcripts:

            transcript_id = transcript.get_transcript_id()
            
            transcript_read_count_total = 0
            read_names = transcript.get_read_names()
            for read_name in read_names:
                num_transcripts_with_assigned_read = len(read_name_to_transcripts[read_name])
                frac_read_assignment = 1 / num_transcripts_with_assigned_read
                transcript_read_count_total += frac_read_assignment
                transcript_to_fractional_read_assignment[transcript_id][read_name] = frac_read_assignment
                
            transcript_to_read_count[transcript_id] = transcript_read_count_total
            transcript_to_expr_val[transcript_id] = transcript_read_count_total / num_mapped_reads if num_mapped_reads > 0 else 0 # * 1e6
            logger.debug(f"-assigning transcript {transcript_id} read count: {transcript_read_count_total} and expr val {transcript_read_count_total}/{num_mapped_reads} = {transcript_to_expr_val[transcript_id]}")

            
        ## DEBUGGING
        if False:
            logger.debug("# Isoform read assignments:\n")
            for read_name in read_name_to_transcripts:
                transcripts_read_assigned = read_name_to_transcripts[read_name]
                logger.debug("read_name {} assigned to {}".format(read_name, transcripts_read_assigned))
                if len(transcripts_read_assigned) > 1:
                    logger.debug("*** Splitting read: {} across {} transcripts: {}".format(read_name, len(transcripts_read_assigned), transcripts_read_assigned))




                
        if run_EM:

            def compute_log_likelihood():
                # compute log likelihood
                log_likelihood = 0
                for transcript_id, reads_to_fracs  in transcript_to_fractional_read_assignment.items():
                    transcript_expr = transcript_to_expr_val[transcript_id]
                    if transcript_expr > 0:
                        for read, frac in reads_to_fracs.items():
                            log_likelihood += frac * log(transcript_expr)

                return log_likelihood
            

            ## compute likelihood.

            ## go through multiple rounds of EM

            prev_log_likelihood = compute_log_likelihood()
            logger.debug("Log likelihood before starting EM: {:.5E}".format(prev_log_likelihood))
            
            for EM_round in range(1, 1000):
                
                logger.debug("EM round {}".format(EM_round))

                ## fractionally assign reads based on expr values
                transcript_to_read_count.clear()

                for transcript in transcripts:
                    transcript_id = transcript.get_transcript_id()
                    transcript_read_count_total = 0
                    read_names = transcript.get_read_names()
                    transcript_expr = transcript_to_expr_val[transcript_id] 
                    for read_name in read_names:
                        transcripts_with_read = read_name_to_transcripts[read_name]
                        sum_expr = 0
                        for tran_with_read in transcripts_with_read:
                            tran_with_read_id = tran_with_read.get_transcript_id()
                            sum_expr += transcript_to_expr_val[tran_with_read_id]
                        frac_read_assignment = transcript_expr / sum_expr
                        transcript_to_read_count[transcript_id] += frac_read_assignment
                        transcript_to_fractional_read_assignment[transcript_id][read_name] = frac_read_assignment

                ## recompute expr_vals
                transcript_to_expr_val.clear()

                for transcript in transcripts:
                    transcript_id = transcript.get_transcript_id()
                    transcript_read_count = transcript_to_read_count[transcript_id]
                    transcript_to_expr_val[transcript_id] = transcript_read_count/num_mapped_reads if num_mapped_reads > 0 else 0 # * 1e6
                    logger.debug(f"-assigning transcript {transcript_id} read count: {transcript_read_count_total} and expr val {transcript_read_count}/{num_mapped_reads} = {transcript_to_expr_val[transcript_id]}")
                    

                log_likelihood = compute_log_likelihood()
                    
                logger.debug("\tlog_likelihood: {:.5E}".format(log_likelihood))
                if prev_log_likelihood is not None:
                    delta = log_likelihood - prev_log_likelihood
                    logger.debug(f"\t\tEM_round[{EM_round}] delta = {delta:.5E}")
                    if delta < 1e-5:
                        logger.debug(f"\t\tEM_round[{EM_round}] delta = {delta:.5E} ...  max likelihood reached. Stopping EM.")
                        break
                    
                prev_log_likelihood = log_likelihood
                

        ## assign final read counts to each transcript object.
        gene_to_transcripts = defaultdict(list)
        for transcript in transcripts:
            transcript_id = transcript.get_transcript_id()
            transcript_read_count = transcript_to_read_count[transcript_id]
            transcript.set_read_counts_assigned(transcript_read_count)
            gene_id = transcript.get_gene_id()
            gene_to_transcripts[gene_id].append(transcript)

        
        ## isoform isoform fraction
        for gene_id in gene_to_transcripts:
            transcripts_of_gene = gene_to_transcripts[gene_id]

            # evaluate isoform fraction
            sum_gene_reads = 0
            for transcript_of_gene in transcripts_of_gene:
                transcript_read_count = transcript_of_gene.get_read_counts_assigned()
                sum_gene_reads += transcript_read_count

            logger.debug("gene_id {} has total reads: {}".format(gene_id, sum_gene_reads))

            for transcript_of_gene in transcripts_of_gene:
                transcript_id =  transcript_of_gene.get_transcript_id()
                transcript_read_count = transcript_of_gene.get_read_counts_assigned()
                isoform_frac = transcript_read_count / sum_gene_reads if sum_gene_reads > 0 else 0
                logger.debug("\ttranscript_id {} has {} reads = {} isoform fraction of {}".format(
                    transcript_id,
                    transcript_read_count,
                    isoform_frac,
                    gene_id))
                transcript_of_gene.set_isoform_fraction(isoform_frac)
                    
        return transcript_to_fractional_read_assignment
        
    
                
    def report_quant_results(self, transcripts, transcript_to_fractional_read_assignment, ofh_quant_vals, ofh_read_tracking):
                
        ## generate final report.

        ## sort descendingly by read support
        transcripts = sorted(transcripts, key=lambda x: (x.get_read_counts_assigned(), x.get_transcript_id()), reverse=True)
        
        
        # first, get sum of reads per gene
        gene_to_read_count = defaultdict(int)
        for transcript in transcripts:
            gene_id = transcript.get_gene_id()
            counts = transcript.get_read_counts_assigned()
            gene_to_read_count[gene_id] += counts

        
            
        for transcript in transcripts:
            transcript_id = transcript.get_transcript_id()
            gene_id = transcript.get_gene_id()
            counts = transcript.get_read_counts_assigned()
            isoform_frac = transcript.get_isoform_fraction()
            
            readnames = transcript.get_read_names()
            readnames = sorted(readnames)

            num_uniquely_assigned_reads = 0
            
            for readname in readnames:
                frac_read_assigned = transcript_to_fractional_read_assignment[transcript_id][readname]
                print("\t".join([gene_id, transcript_id, readname, "{:.3f}".format(frac_read_assigned)]), file=ofh_read_tracking)
                if frac_read_assigned == 1:
                    num_uniquely_assigned_reads += 1

            gene_read_count = gene_to_read_count[gene_id]
            unique_gene_read_fraction = num_uniquely_assigned_reads / gene_read_count if gene_read_count > 0 else 0

            report_txt = "\t".join([gene_id,
                                    transcript_id,
                                    f"{num_uniquely_assigned_reads}",
                                    f"{counts:.1f}",
                                    f"{isoform_frac:.3f}",
                                    f"{unique_gene_read_fraction:0.3f}"])
            
            logger.info(report_txt)
            print(report_txt, file=ofh_quant_vals)
            
                
            """
            if (DEBUG):
                print("transcript_id\t{}\n{}".format(transcript_id, transcript._simplepath), file=ofh_read_tracking)
                for readname in readnames:
                    print("read:\t{}\n{}".format(readname,
                                                 self._read_name_to_multipath[readname].get_simple_path()),
                          file=ofh_read_tracking)
                print("\n", file=ofh_read_tracking)
            else:
                print("\t".join([gene_id, transcript_id, ",".join(readnames)]), file=ofh_read_tracking)
            """

        return 
    



    @staticmethod
    def get_gene_read_counts(frac_read_assignments, transcript_id_to_transcript_obj):
        gene_id_to_read_count = defaultdict(int)
        for transcript_id, transcript_read_frac_assignments in frac_read_assignments.items():
            gene_id = transcript_id_to_transcript_obj[transcript_id].get_gene_id()
            for read, frac_assigned in transcript_read_frac_assignments.items():
                gene_id_to_read_count[gene_id] += frac_assigned

        return gene_id_to_read_count
        


    

    @staticmethod
    def filter_isoforms_by_min_isoform_fraction(transcripts, min_isoform_fraction, run_EM):

        min_frac_gene_unique_reads = LRAA_Globals.config['min_frac_gene_unique_reads']
        
        logger.info("Filtering transcripts according to min isoform fraction: {}".format(min_isoform_fraction))

        transcript_id_to_transcript_obj = dict([ (x.get_transcript_id(), x) for x in transcripts] )
        

        def get_idoform_unique_assigned_read_count(transcript_id, frac_read_assignments):
            num_unique_reads = 0
            for read in frac_read_assignments[transcript_id]:
                if frac_read_assignments[transcript_id][read] >= 0.9999: # close enough to 1.0
                    num_unique_reads += 1

            return num_unique_reads
                                
        
        isoforms_were_filtered = True # init for loop

        q = Quantify()


        filtering_round = 0

        frac_read_assignments = None


        
        while isoforms_were_filtered:

            filtering_round += 1
            num_total_isoforms = 0
            num_filtered_isoforms = 0
            transcripts_retained = list()
            isoforms_were_filtered = False # update to True if we do filter an isoform out.

            frac_read_assignments = q._estimate_isoform_read_support(transcripts, run_EM)
            gene_id_to_read_count = Quantify.get_gene_read_counts(frac_read_assignments, transcript_id_to_transcript_obj)
            
            genes_represented = set()
            for transcript in transcripts:
                num_total_isoforms += 1
                transcript_id = transcript.get_transcript_id()
                gene_id = transcript.get_gene_id()
                genes_represented.add(gene_id)
                gene_read_count = gene_id_to_read_count[gene_id]
                transcript_unique_read_count = get_idoform_unique_assigned_read_count(transcript_id, frac_read_assignments)
                
                frac_gene_unique_reads = transcript_unique_read_count / gene_read_count
                
                logger.debug("Transcript_id: {} has unique read frac of gene total reads: {}".format(transcript_id, frac_gene_unique_reads))

                if (not isoforms_were_filtered and
                    (frac_gene_unique_reads < min_frac_gene_unique_reads
                      or
                       transcript.get_isoform_fraction() < min_isoform_fraction)
                    ):
                    
                    isoforms_were_filtered = True
                    num_filtered_isoforms += 1

                    logger.debug("Filtering out transcript_id {} as low fraction of unique reads: {}".format(transcript_id, frac_gene_unique_reads))
                    
                else:
                    transcripts_retained.append(transcript)

            
            logger.debug("isoform filtering round {} involved filtering of {} isoforms / {} total isoforms of {} genes".format(
                filtering_round,
                num_filtered_isoforms,
                num_total_isoforms,
                len(genes_represented)) )
            
            # reset list of transcripts
            transcripts = transcripts_retained
            
        return (transcripts, frac_read_assignments)


    @staticmethod
    def prune_likely_degradation_products(transcripts, splice_graph, run_EM):

        logger.info("Pruning likely degradation products")
        
        sg = splice_graph
        
        # run an initial quant.
        q = Quantify()
        frac_read_assignments = q._estimate_isoform_read_support(transcripts, run_EM)

        transcript_id_to_transcript_obj = dict([ (x.get_transcript_id(), x) for x in transcripts] )
        gene_read_counts = Quantify.get_gene_read_counts(frac_read_assignments, transcript_id_to_transcript_obj)
        
        transcripts_ret = list() # transcripts not pruned and returned.
                
        # first organize by gene_id
        gene_id_to_transcripts = defaultdict(set)
        for transcript in transcripts:
            gene_id = transcript.get_gene_id()
            gene_id_to_transcripts[gene_id].add(transcript)


        for gene_id, transcript_set in gene_id_to_transcripts.items():

            if len(transcript_set) == 1:
                transcripts_ret.extend(list(transcript_set))
                continue
            
            # compare isoforms for compatibility ignoring TSS and PolyA sites
            # if an isoform is fully contained by another and has substantially less expression, prune it.

            transcript_list = list(transcript_set)
            contig_strand = transcript_list[0].get_strand()


            """  #TODO: figure out sorting order to speed this up.
            transcript_list = sorted(transcript_list, key=lambda x: (x.get_coords()[0], x.get_coords()[1], x.get_left_boundary_sort_weight(), x.get_right_boundary_sort_weight()))
            
            if contig_strand == '-':
                transcript_list = list(sorted(transcript_list, key=lambda x: (-1 * x.get_coords()[1],
                                                                              -1 * x.get_coords()[0],
                                                                              -1 * x.get_right_boundary_sort_weight(),
                                                                              -1 * x.get_left_boundary_sort_weight()
                                                                           ) ) ) 
            """

            # sort by desc cdna len
            transcript_list = list(reversed(sorted(transcript_list, key=lambda x: x.get_cdna_len())))
            
            transcript_prune_as_degradation = set()
            for i in range(len(transcript_list)):
                transcript_i = transcript_list[i]
                transcript_i_id = transcript_i.get_transcript_id()

                gene_i_id = transcript_i.get_gene_id()

                gene_read_count = gene_read_counts[gene_i_id]
                
                if transcript_i in transcript_prune_as_degradation:
                    continue

                
                transcript_i_simple_path = transcript_i.get_simple_path()
                i_path_trimmed, i_TSS_id, i_polyA_id = SPU.trim_TSS_and_PolyA(transcript_i_simple_path, contig_strand)
                transcript_i_read_counts_assigned = transcript_i.get_read_counts_assigned()
                
                #for j in range(i+1, len(transcript_list)):

                for j in range(len(transcript_list)):

                    if i==j:
                        continue
                    
                    transcript_j = transcript_list[j]

                    if transcript_j in transcript_prune_as_degradation:
                        continue

                    transcript_j_id = transcript_j.get_transcript_id()
                    gene_j_id = transcript_j.get_gene_id()

                    assert gene_i_id == gene_j_id
                    
                    transcript_j_simple_path = transcript_j.get_simple_path()
                    j_path_trimmed, j_TSS_id, j_polyA_id = SPU.trim_TSS_and_PolyA(transcript_j_simple_path, contig_strand)
                    transcript_j_read_counts_assigned = transcript_j.get_read_counts_assigned()

                    if len(j_path_trimmed) > len(i_path_trimmed):
                        # no way can i subsume j
                        continue
                    
                    # Simple_path_utils.simple_paths_overlap_and_compatible_spacefree_region_path_A(self.get_splice_graph(), my_path, other_path)
                    frac_expression_i = transcript_j_read_counts_assigned / transcript_i_read_counts_assigned

                    frac_gene_expression_j = transcript_j_read_counts_assigned / gene_read_count

                    
                    logger.debug("Exploring path: {} as subsuming {}".format(transcript_i_simple_path, transcript_j_simple_path))
                    
                    if SPU.path_A_contains_path_B(i_path_trimmed, j_path_trimmed):

                        logger.debug("splice compatible & contained transcript_j_id {} has frac gene expression: {}".format(transcript_j_id, frac_gene_expression_j))
                                     
                        subsume_J = False
                        

                        if LRAA_Globals.config['collapse_alt_TSS_and_PolyA']:
                            logger.debug("Collapsing compatible path: {} into {}".format(transcript_j, transcript_i))
                            subsume_J = True


                        elif j_TSS_id is not None:
                            

                            j_TSS_read_count = sg.get_node_obj_via_id(j_TSS_id).get_read_support()

                            frac_gene_express_j_TSS = j_TSS_read_count / gene_read_count
                            logger.debug("frac_gene_express_j_TSS {} {} : {:.3f}".format(transcript_j_id, transcript_j_simple_path, frac_gene_express_j_TSS))

                            
                            if frac_gene_express_j_TSS < LRAA_Globals.config['min_frac_gene_alignments_define_TSS_site']:
                                logger.debug(
                                    "based on j_TSS count frac_gene_expression: {:.3f}, path_i: {} is subsuming path_j: {}".format(
                                        frac_gene_express_j_TSS, transcript_i_simple_path, transcript_j_simple_path)
                                        )
                                subsume_J = True
                                
                                                            
                            elif i_TSS_id is not None:

                                i_TSS_read_count = sg.get_node_obj_via_id(i_TSS_id).get_read_support()

                                frac_i_TSS =  j_TSS_read_count/i_TSS_read_count
                                logger.debug("frac_i_TSS: {:.3f} of path_j: {} to path_i{}".format(frac_i_TSS, transcript_j_simple_path, transcript_i_simple_path))


                                if frac_i_TSS < LRAA_Globals.config['max_frac_alt_TSS_from_degradation']:
                                    logger.debug("based on frac_i_TSS: {:.3f}, path_i: {} is subsuming path_j: {}".format(frac_i_TSS, transcript_i_simple_path, transcript_j_simple_path))
                                    subsume_J = True

                            elif i_TSS_id is None:
                                subsume_J = False


                        elif i_TSS_id is not None and j_TSS_id is None:
                            # no j_TSS but have i_TSS
                            subsume_J = True


                        else:
                            # neither has a TSS assigned.
                            subsume_J = True


                            
                        #####  PolyA check ######    
                        ## But dont subsume if they have polyA and they differ
                        if subsume_J:
                            # here we might resurrect it based on polyA status
                            if frac_gene_expression_j >= LRAA_Globals.config['min_frac_alignments_define_polyA_site']: 
                            
                                if i_polyA_id is not None and j_polyA_id is not None:
                                    if i_polyA_id != j_polyA_id:
                                        subsume_J = False
                                        logger.debug("resurrecting {} based on alt polyA".format(transcript_j_id))
                                elif j_polyA_id is not None:
                                    logger.debug("resurrecting {} based on defined polyA".format(transcript_j_id))
                                    subsume_J = False

                        
                        if subsume_J:
                            logger.debug("Pruning {} as likely degradation product of {}".format(transcript_j, transcript_i))
                            transcript_prune_as_degradation.add(transcript_j)
                            transcript_i.add_read_names(transcript_j.get_read_names())


                            
            # retain the ones not pruned
            for transcript in transcript_set:
                if transcript not in transcript_prune_as_degradation:
                    transcripts_ret.append(transcript)


        # after pruning transcripts, rerun quant
        frac_read_assignments = q._estimate_isoform_read_support(transcripts_ret, run_EM)
        
        return (transcripts_ret, frac_read_assignments)

