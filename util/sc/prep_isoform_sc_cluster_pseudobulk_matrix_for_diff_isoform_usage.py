#!//usr/bin/env Rscript


suppressPackageStartupMessages(library("argparse"))
suppressPackageStartupMessages(library("tidyverse"))

parser = ArgumentParser()
parser$add_argument("--sc_isoform_pseudobulk_matrix", help="sc isoform pseudobulk count matrix", required=TRUE, nargs=1)
parser$add_argument("--gene_trans_info", help="gene trans info file", required=TRUE, nargs=1)
parser$add_argument("--output_matrix", help="name for output matrix file", required=TRUE, nargs=1)

args = parser$parse_args()
input_isoform_cluster_count_matrix = args$sc_isoform_pseudobulk_matrix
gene_trans_info_filename = args$gene_trans_info

input_counts_matrix = read.csv(input_isoform_cluster_count_matrix, header=T, row.names=1, sep="\t")
transcript_ids = rownames(input_counts_matrix)
cluster_names = colnames(input_counts_matrix)

input_counts_matrix$transcript_id = transcript_ids


gene_trans_info = read.csv(gene_trans_info_filename, header=T, sep="\t")

counts_matrix = left_join(input_counts_matrix, gene_trans_info,
                          by='transcript_id')

counts_matrix = counts_matrix %>% select(gene_id, transcript_id, all_of(cluster_names))


write.table(counts_matrix, file=args$output_matrix, sep="\t", quote=F, row.names=FALSE)

quit(save = "no", status = 0, runLast = FALSE)

