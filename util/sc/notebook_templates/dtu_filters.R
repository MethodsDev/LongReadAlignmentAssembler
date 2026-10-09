library(tidyverse)

# LRAA can join neighboring same-strand genes into one gene component through
# read-through (cis-fusion) transcripts; e.g. PTPRCAP and CORO1B share
# PTPRCAP^g:chr11:-:comp-799. The DTU tests group by the component's gene symbol,
# so such a group compares isoforms of two different genes and reports a shift in
# relative gene expression as differential isoform usage.
#
# Each isoform keeps its own symbol as the text before '^' in its transcript id.
# Rows whose dominant and alternate isoforms carry more than one such symbol are
# dropped. Ids without '^' are unannotated and don't count as a second gene.

exclude_cross_gene_DTU_comparisons = function(dtu_df) {

    n_symbols = map2_int(dtu_df$dominant_transcript_ids, dtu_df$alternate_transcript_ids,
                         function(dominant_ids, alternate_ids) {
                             ids = unlist(str_split(c(dominant_ids, alternate_ids), ","))
                             symbols = str_match(ids, "^([^^]+)\\^")[, 2]
                             n_distinct(na.omit(symbols))
                         })

    cross_gene = n_symbols > 1

    message("excluding ", sum(cross_gene), " of ", nrow(dtu_df),
            " DTU rows that compare isoforms of different genes, across ",
            n_distinct(dtu_df$gene_symbol[cross_gene]), " gene groups")

    dtu_df[! cross_gene, ]
}
