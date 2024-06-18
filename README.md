# LRAA - Long Read Alignment Assembler

** Under Development **

Utility for resolving and/or quantification of isoforms from aligned PacBio Iso-seq reads.


## Isoform Discovery

Given a bam file from aligned reads (using minimap2), perform isoform discovery like so:

    LRAA --genome genome.fasta --bam aligned_IsoSeq.mm2.bam



## Isoform Quantification

    LRAA --genome genome.fasta --bam aligned_IsoSeq.mm2.bam --gtf target_isoforms.gtf --quant_only


    
