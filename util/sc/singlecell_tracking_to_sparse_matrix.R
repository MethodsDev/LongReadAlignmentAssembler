#!/usr/bin/env Rscript


suppressPackageStartupMessages(library("argparse"))

parser = ArgumentParser()
parser$add_argument("--tracking", help="input data file quant.tracking", required=TRUE, nargs=1)
parser$add_argument("--output_prefix", help="output prefix", required=TRUE, nargs=1)

args = parser$parse_args()    
dat_filename = args$tracking
output_prefix = args$output_prefix

library(tidyverse)
library(Matrix)

message("-reading in ", dat_filename)
data = read.csv(dat_filename, sep="\t", header=T)

# extract cell barcode and umi info
data = data %>% separate_wider_delim(read_name, "^", names=c("cell_barcode", "UMI", "core_read_name"))


# get gene/cell counts
gene_cell_counts = data %>% group_by(gene_id, cell_barcode) %>% summarize(sum_UMIs = sum(frac_assigned) ) %>% ungroup()

message("-writing gene cell counts table")
gene_cell_counts_tsv = paste0(output_prefix, ".gene_cell_counts.tsv")
write.table(gene_cell_counts, file=gene_cell_counts_tsv, sep="\t", row.names=F, quote=F)

# get isoform/cell counts
isoform_cell_counts = data	%>% group_by(transcript_id, cell_barcode) %>% summarize(sum_UMIs = sum(frac_assigned) ) %>% ungroup()
message("-writing isoform cell counts table")
isoform_cell_counts_tsv = paste0(output_prefix, ".isoform_cell_counts.tsv")
write.table(isoform_cell_counts, file=isoform_cell_counts_tsv, sep="\t", row.names=F, quote=F)


make_sparse_matrix_outputs = function(counts_data, outdirname) {
   message("-making sparse matrix outputs for: ", outdirname)

   if (! dir.exists(outdirname)) {
      dir.create(outdirname, recursive = TRUE)
   }

   colnames(counts_data) = c('feature_id', 'cell_barcode', 'UMI_counts')

   feature_ids = counts_data %>% select(feature_id) %>% unique() %>% pull(feature_id)
   cell_barcodes = counts_data %>% select(cell_barcode) %>% unique() %>% pull(cell_barcode)

   counts_data$feature_id = factor(counts_data$feature_id, levels=feature_ids)
   counts_data$cell_barcode = factor(counts_data$cell_barcode, levels=cell_barcodes)

   counts_data$feature_id = as.numeric(counts_data$feature_id)
   counts_data$cell_barcode = as.numeric(counts_data$cell_barcode)

   sparseM = sparseMatrix(j=counts_data$cell_barcode, i=counts_data$feature_id, x=counts_data$UMI_counts)

   Matrix::writeMM(obj =  sparseM, file =  paste0(outdirname, "/matrix.mtx"))
   write(feature_ids, file = paste0(outdirname, "/features.tsv"))
   write(cell_barcodes, file = paste0(outdirname, "/barcodes.tsv"))

   system(paste0("gzip ", outdirname, "/*")) 
   
   message("done with ", outdirname)


}

make_sparse_matrix_outputs(gene_cell_counts, paste0(output_prefix, "^gene-sparseM") )

make_sparse_matrix_outputs(isoform_cell_counts, paste0(output_prefix, "^isoform-sparseM") )


message("all done.")



