# LRAA - Long Read Alignment Assembler

Isoform discovery and quantification based on long read isoform sequence alignments (PacBio or ONT). 

LRAA has three modes, described below:

- De novo (reference annotation-free) isoform identification and quantification
- Reference-guided isoform detection and quantification
- Isoform expression quantification only



## Isoform Discovery

### De novo Reference Annotation-free Isoform Discovery (and quantification)

Given a bam file from aligned reads (using minimap2), perform isoform discovery like so:

    LRAA --genome genome.fasta --bam aligned_IsoSeq.mm2.bam

>If using with reads that are not PacBio HiFi and have error rates that are > 2%, use the --LowFi parameter.  By default, any read alignments with <98% identity are ignored.

### Reference Annotation-guided Isoform Discovery (and quantification)

    LRAA --genome genome.fasta --bam aligned_IsoSeq.mm2.bam --gtf reference_annotation.gtf

>Note that input refernece annotations not found with evidence of expression are excluded from the output. Also, reference annotation structures may be extended if evidence supports alternative TSS or PolyA sites. Transcript identfiers will all be reassigned. Use GFFcompare to the reference to determine relationships to the input ref annotations.

>If using with reads that are not PacBio HiFi and have error rates that are > 2%, use the --LowFi parameter.  By default, any read alignments with <98% identity are ignored.    

## Isoform Quantification Only

    LRAA --genome genome.fasta --bam aligned_IsoSeq.mm2.bam --gtf target_isoforms.gtf --quant_only

>No novel isoform detection is performed. All original gene_id and transcript_ids are retained in the final output, including those with no evidence of expression (0 reads and 0 TPM).
    
>If using with reads that are not PacBio HiFi and have error rates that are > 2%, use the --LowFi parameter.  By default, any read alignments with <98% identity are ignored.


## Isoform ID and/or Quant for single cell analyses

When using LRAA with single cell data, it's required that the cell barcodes and UMIs are encoded as 'CB' and 'XM' annotations, respectively, in the minimap2 aligned bam file.

Run LRAA in (pseudo)bulk mode (only current execution mode). Once LRAA completes, one of the output files will be a '*.tracking' file that includes read assignments and incorporates the cell barcode and UMI tags.

To construct a cell-by-transcript expression matrix, run the included script:

    util/sc/singlecell_tracking_to_sparse_matrix.R --tracking LRAA.tracking --output_prefix sample_name

The resulting sparse matrix can be used with tools such as [Seurat](https://satijalab.org/seurat/) for single cell analyses.


## LRAA is available on Dockstore

Find LRAA on Dockstore [here](https://dockstore.org/workflows/github.com/MethodsDev/LongReadAlignmentAssembler/LRAA:main?tab=info).
