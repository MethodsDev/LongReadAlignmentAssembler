import sys, os, re
from collections import defaultdict
from GenomeFeature import GenomeFeature
import PASA_scored_path


class Transcript (GenomeFeature):

    trans_id_counter = 0
    

    def __init__(self, contig_acc, segment_coordinates_list, orient):

        segment_coordinates_list = sorted(segment_coordinates_list, key=lambda x: x[0])

        trans_lend = segment_coordinates_list[0][0]
        trans_rend = segment_coordinates_list[-1][1]

        super().__init__(contig_acc, trans_lend, trans_rend, orient)

        self._orient = orient
        self._exon_segments = segment_coordinates_list

        Transcript.trans_id_counter += 1
        
        self._id = ":".join([contig_acc, str(trans_lend), str(trans_rend), "N:{}".format(Transcript.trans_id_counter), orient])

        self._transcript_id = "t.{}".format(self._id) # should reset to something more useful
        self._gene_id = "g.{}".format(self._id) # ditto above

        self._meta = dict()

        self._scored_path_obj = None #optional - useful if transcript obj was built based on a scored path

        self.read_names = list() # list of read names supporting the transcript structure.

        self._multipath = None # multipath obj

        self._simplepath = None

        self._read_counts_assigned = None # set during expression quantification

        self._isoform_fraction = None # set during expression quantification

        self._cdna_len = 0
        for exon_segment in segment_coordinates_list:
            self._cdna_len += exon_segment[1] - exon_segment[0] + 1
        
        return

    ## getters

    def get_exon_segments(self):
        return self._exon_segments.copy()

    def get_introns(self):
        intron_coordsets = list()
        exon_segments = self.get_exon_segments()
        if len(exon_segments) > 1:
            exon_segments = sorted(exon_segments, key=lambda x: x[0])
            for i in range(1, len(exon_segments)):
                intron_lend = exon_segments[i-1][1]+1
                intron_rend = exon_segments[i][0] - 1
                assert intron_lend < intron_rend
                intron_coordsets.append( (intron_lend, intron_rend) )
        return intron_coordsets
    
    
    def get_strand(self):
        return(self._orient)

    def get_orient(self):
        return(self.get_strand())

    def get_cdna_len(self):
        return self._cdna_len

    def get_transcript_id(self):
        if self._transcript_id is not None:
            return self._transcript_id
        
        elif 'transcript_id' in self._meta:
            return(self._meta['transcript_id'])

        else:
            return(self._id)

    
    def get_gene_id(self):
        if self._gene_id is not None:
            return self._gene_id
        
        elif 'gene_id' in self._meta:
            return(self._meta['gene_id'])
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
            
    
    def __repr__(self):

        text = "Transcript: {} {}-{} [{}] {} {} segs: {}".format(self._contig_acc,
                                                                 self._lend,
                                                                 self._rend,
                                                                 self._orient,
                                                                 self.get_transcript_id(),
                                                                 self.get_gene_id(),
                                                                 self._exon_segments)
        if self._simplepath is not None:
            text += " {} ".format(self._simplepath)
        
        if self._meta is not None and len(self._meta) > 0:
            text += "\t" + str(self._meta)

        return text
        
        
    def set_scored_path_obj(self, scored_path_obj):
        assert(type(scored_path_obj) == PASA_scored_path.PASA_scored_path)
        self._scored_path_obj = scored_path_obj
        
        
    def set_gene_id(self, gene_id):
        self._gene_id = gene_id

    def set_transcript_id(self, transcript_id):
        self._transcript_id = transcript_id
        

    def add_meta(self, meta_key, meta_val=None):
        
        if self._meta == None:
            self._meta = dict()


        if meta_val is None and type(meta_key)==dict:
            self._meta = meta_key.copy()

        elif meta_val is not None:
            self._meta[meta_key] = meta_val
        else:
            raise RuntimeError("Error: not sure how to handle input params")
            
        return

    def get_read_names(self):
        return self.read_names.copy()
    

    def add_read_names(self, read_names):
        if self.read_names == None:
            self.read_names = list()

        if type(read_names) in (list, set):
            self.read_names.extend(list(read_names))
        else:
            self.read_names.append(read_names)
        

    def set_read_counts_assigned(self, read_counts):
        self._read_counts_assigned = read_counts

    def get_read_counts_assigned(self):
        assert self._read_counts_assigned is not None, "Error, read counts assigned is None - maybe quant not run yet?"
        return self._read_counts_assigned

    def set_isoform_fraction(self, frac):
        self._isoform_fraction = frac

    def get_isoform_fraction(self):
        assert self._isoform_fraction is not None, "Error, isoform fraction is None - maybe quant not run yet?"
        return self._isoform_fraction
    

    def lighten(self):
        # lighten transcript by removing nonessential memory allocs
        self._scored_path_obj = None
        self._multipath = None
        

    def to_GTF_format(self):

        ## transcript line:

        gtf_text = ""
        
        if self.read_names:
            gtf_text = f"#{self._id}\t" + ",".join(self.read_names) + "\n"
        
        gtf_text += "\t".join([self._contig_acc,
                              "PASA-SALRAA",
                              "transcript",
                              str(self._lend),
                              str(self._rend),
                              ".",
                              self._orient,
                              ".",
                              "gene_id \"{}\"; transcript_id \"{}\";".format(self.get_gene_id(), self.get_transcript_id())])

        if self._meta:
            for meta_key in sorted(self._meta.keys()):
                gtf_text += " {} \"{}\";".format(meta_key, self._meta[meta_key])

        gtf_text += "\n"

        for segment in self._exon_segments:
            gtf_text += "\t".join([self._contig_acc,
                                   "PASA-SALRAA",
                                   "exon",
                                   str(segment[0]),
                                   str(segment[1]),
                                   ".",
                                   self._orient,
                                   ".",
                                   "gene_id \"{}\"; transcript_id \"{}\";".format(self.get_gene_id(), self.get_transcript_id())]) + "\n"


        if self._scored_path_obj:
            # include construction info as comments.
            gtf_text += "# derived from scored path obj:\n"

            mpgn_list = self._scored_path_obj.get_path_mpgn_list()
            for mpgn in mpgn_list:
                gtf_text += "# " + str(mpgn) + "\n"

            gtf_text += "\n"

        
        return gtf_text
        
        
        
class GTF_contig_to_transcripts:

    @classmethod
    def parse_GTF_to_Transcripts(cls, gtf_filename):

        gene_id_to_meta = defaultdict(dict)
        transcript_id_to_meta = defaultdict(dict)
        transcript_id_to_genome_info = defaultdict(dict)
        
        with open(gtf_filename, "rt") as fh:
            for line in fh:
                if line[0] == "#":
                    continue
                line = line.rstrip()
                
                if not re.match("\\w", line):
                    continue
                vals = line.split("\t")
                if len(vals) < 9:
                    logger.warn("GTF line has fewer fields than expected. Skipping. {}".format(line))
                    continue
                                
                contig = vals[0]
                feature_type = vals[2]
                lend = int(vals[3])
                rend = int(vals[4])
                strand = vals[6]
                info = vals[8]
                
                info_dict = cls._parse_info(info)

                if feature_type == 'gene':
                    gene_id = info_dict['gene_id']
                    gene_id_to_meta[gene_id] = info_dict
                    
                if 'transcript_id' not in info_dict:
                    continue

                if feature_type != 'exon':
                    continue

                
                transcript_id = info_dict['transcript_id']
                transcript_id_to_meta[transcript_id] = info_dict

                transcript_id_to_genome_info[transcript_id]['contig'] = contig
                transcript_id_to_genome_info[transcript_id]['strand'] = strand
                if ('coords' in transcript_id_to_genome_info[transcript_id].keys()):
                    transcript_id_to_genome_info[transcript_id]['coords'].append([lend, rend])
                else:
                    transcript_id_to_genome_info[transcript_id]['coords'] = [ [lend,rend] ]

        # convert to transcript objects

        contig_to_transcripts = defaultdict(list)

        for transcript_id in transcript_id_to_genome_info:
            transcript_info_dict = transcript_id_to_genome_info[transcript_id]
            contig = transcript_info_dict['contig']
            strand = transcript_info_dict['strand']
            coords_list = transcript_info_dict['coords']

            transcript_meta = transcript_id_to_meta[transcript_id]
            gene_id = transcript_meta['gene_id']
            gene_meta = gene_id_to_meta[gene_id]
            transcript_meta.update(gene_meta)
            
            transcript_obj = Transcript(contig, coords_list, strand)
            transcript_obj.add_meta(transcript_meta)
            
            transcript_obj.set_gene_id(gene_id)
            transcript_obj.set_transcript_id(transcript_id)
            
            contig_to_transcripts[contig].append(transcript_obj)

            
        return contig_to_transcripts


    
    # private
    @classmethod
    def _parse_info(cls, info):

        info_dict = dict()
        
        parts = info.split(";")
        for part in parts:
            part = part.strip()
            m = re.match('^(\S+) \\"([^\\"]+)\\"', part)
            if m:
                token = m.group(1)
                val = m.group(2)

                info_dict[token] = val

        return info_dict

        
if __name__=='__main__':

    # testing gtf parser
    usage = "usage: {} gtf_filename\n\n".format(sys.argv[0])

    if len(sys.argv) < 2:
        exit(usage)

    gtf_filename = sys.argv[1]
    
    contig_to_transcripts = GTF_contig_to_transcripts.parse_GTF_to_Transcripts(gtf_filename)

    for contig, transcript_list in contig_to_transcripts.items():
        for transcript_obj in transcript_list:
            print("\t".join([contig, str(transcript_obj)]))


    sys.exit(0)


