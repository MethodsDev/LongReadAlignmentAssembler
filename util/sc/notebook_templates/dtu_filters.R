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


# The DTU tests also group by gene symbol, so a comparison can pair two models that share
# the symbol but not the locus: separate LRAA components, or one model sitting in the
# other's intron. GPM6A's top-ranked ALS comparison pairs comp-2379:iso-6 (the gene body)
# with comp-2423:iso-1, a two-exon model ~170 kb upstream; a usage shift between disjoint
# transcription units is not an isoform choice.
#
# Rows whose dominant and alternate models share less than min_exonic_overlap of the
# shorter model's exonic length are dropped (the same measure as transcript_pair_overlap /
# get_top_dtu_pair; a side carrying several comma-joined models counts its best-overlapping
# member). Models are looked up in gtf_parsed (exon rows: isoforms, or splice patterns in a
# splice-pattern gtf); rows with a model absent from it are kept and counted.

exclude_nonoverlapping_DTU_pairs = function(dtu_df, min_exonic_overlap = 0.05, gtf = gtf_parsed) {

    split_ids = function(x) unlist(str_split(x, ","))
    ids = unique(c(split_ids(dtu_df$dominant_transcript_ids), split_ids(dtu_df$alternate_transcript_ids)))
    ex = gtf %>% filter(feature == "exon", transcript_id %in% ids) %>%
        transmute(transcript_id, seqname, start = as.integer(start), end = as.integer(end))
    ex_by = split(ex, ex$transcript_id)

    exonic_frac = function(a, b) {
        A = ex_by[[a]]; B = ex_by[[b]]
        if (is.null(A) || is.null(B)) return(NA_real_)
        if (A$seqname[1] != B$seqname[1]) return(0)
        hit = 0
        for (i in seq_len(nrow(A))) hit = hit + sum(pmax(0, pmin(A$end[i], B$end) - pmax(A$start[i], B$start) + 1))
        hit / min(sum(A$end - A$start + 1), sum(B$end - B$start + 1))
    }

    pairs = dtu_df %>% distinct(dominant_transcript_ids, alternate_transcript_ids)
    pairs$.exonic_overlap = map2_dbl(pairs$dominant_transcript_ids, pairs$alternate_transcript_ids, function(d, a) {
        v = unlist(lapply(split_ids(d), function(x) sapply(split_ids(a), function(y) exonic_frac(x, y))))
        if (all(is.na(v))) NA_real_ else max(v, na.rm = TRUE)
    })

    d = dtu_df %>% left_join(pairs, by = c("dominant_transcript_ids", "alternate_transcript_ids"))
    drop = ! is.na(d$.exonic_overlap) & d$.exonic_overlap < min_exonic_overlap

    message("excluding ", sum(drop), " of ", nrow(d), " DTU rows whose two models share < ",
            min_exonic_overlap, " exonic overlap (", n_distinct(d$gene_symbol[drop]), " genes; ",
            sum(is.na(d$.exonic_overlap)), " rows with a model not in the gtf kept)")

    d[! drop, ] %>% select(-.exonic_overlap)
}
