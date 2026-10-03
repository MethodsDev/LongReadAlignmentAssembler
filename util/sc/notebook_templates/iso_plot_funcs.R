library(tidyverse)
library(Seurat)
library(pheatmap)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# define inputs

# sample_name = "PBMCs_LRAA-isoforms"

# sparse_matrix_data_dir = "../../LRAA_sc_PBMCs^isoform-sparseM/"

# umap_cluster_file = "../../LRAA_sc_PBMCs.genes-cell_cluster_assignments.wUMAP.wCAS.tsv"

# cluster_pseudobulk_matrix_filename = "LRAA_sc_PBMCs^isoform-sparseM.clusters_pseudobulk.matrix"

# gtf_filename = "../../PBMCs_pbio.LRAA.sc_merged.gtf.updated.gtf.segmented.gtf"

# diff_iso_usage_stats_filename = "LRAA_sc_PBMCs^isoform-sparseM.top_only_w_recip_delta_pi.diff_iso.tsv.signif_only"

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


parse_inputs = function(sample_name,
                        sparse_matrix_data_dir,
                        umap_cluster_file,
                        cluster_pseudobulk_matrix_filename,
                        gtf_filename,
                        diff_iso_usage_stats_filename,
                        cell_fractions_matrix_filename = NULL) {
    
    
    message("-loading sparse matrix data: ", sparse_matrix_data_dir)
    isoform_expr_data <<- Read10X(data.dir=sparse_matrix_data_dir,
                                gene.column = 1,
                                cell.column = 2,
                                unique.features = TRUE,
                                strip.suffix = FALSE)
    

    message("-reading umap: ", umap_cluster_file)
    umap_df <<- read.csv(umap_cluster_file, header=T, sep="\t")
    
    base_umap <<- umap_df %>% ggplot(aes(x=umap_1, y=umap_2)) + geom_point(color='gray', alpha=0.1) + theme_void()
    
    message("-parsing cluster pseudobulk matrix: ", cluster_pseudobulk_matrix_filename)
    cluster_counts_matrix <<- read.csv(cluster_pseudobulk_matrix_filename, header=T, row.names=1, sep="\t")
    
    message("-making CPM matrix")
    cluster_CPM_matrix <<- sweep(cluster_counts_matrix, 2, colSums(cluster_counts_matrix), "/") * 1e6

    # Cells per cluster, used to turn cell FRACTIONS into cell COUNTS. Needed because a
    # fraction-of-cells threshold is biased toward small clusters in exactly the same way
    # the isoform fraction is: 2 cells of 92 (2.2%) outranks 7 cells of 1248 (0.56%).
    cluster_cell_counts <<- umap_df %>% count(seurat_clusters) %>%
        { setNames(.$n, paste0("Cluster_", .$seurat_clusters)) }

    # Fraction of cells per cluster in which each feature is detected. Optional: when
    # absent, the cell-support filter in
    # get_expression_ggplot2_heatmap_w_exon_structures is unavailable and is skipped.
    cluster_cell_fraction_matrix <<- NULL
    if (! is.null(cell_fractions_matrix_filename)) {
        message("-parsing cell fractions expressed matrix: ", cell_fractions_matrix_filename)
        cluster_cell_fraction_matrix <<- read.csv(cell_fractions_matrix_filename,
                                                  header=T, row.names=1, sep="\t")
    }
    
    
    message("-parsing diff iso usage stats: ", diff_iso_usage_stats_filename)
    diff_iso_usage_stats <<- read.csv(diff_iso_usage_stats_filename, sep="\t", header=T)
    
    
    message("-parsing gtf file: ", gtf_filename)
    gtf_parsed <<- parse_gtf_file(gtf_filename)
    
}

parse_gtf_file = function(gtf_filename) {

    gtf_parsed_file = paste0(gtf_filename, ".rds")
    
    # Check if cached file exists and is newer than source GTF
    if (file.exists(gtf_parsed_file) && 
        file.mtime(gtf_parsed_file) > file.mtime(gtf_filename)) {
        
        gtf_parsed = readRDS(gtf_parsed_file)
        
    } else {
        
        # Read the GTF file
        gtf_parsed <- read_tsv( gtf_filename,
                                comment = "#", 
                                col_names = FALSE, 
                                col_types = cols(.default = "c"))
        
        # Assign standard GTF column names
        colnames(gtf_parsed)[1:9] <- c("seqname", "source", "feature", "start", "end", 
                                       "score", "strand", "frame", "attribute")
        
        # Function to parse attributes into a named list
        parse_attributes <- function(attr_str) {
            attrs <- str_split(attr_str, ";\\s*")[[1]]
            kv_pairs <- str_match(attrs, '^(\\S+)\\s+"([^"]+)"')
            kv_pairs <- kv_pairs[!is.na(kv_pairs[, 1]), , drop = FALSE]
            if (nrow(kv_pairs) == 0) return(named(list()))
            set_names(kv_pairs[, 3], kv_pairs[, 2])
        }
        
        # Safely parse attributes
        gtf_parsed$attr_list <- map(gtf_parsed$attribute, parse_attributes)
        
        # Turn the named list column into separate columns
        gtf_parsed = gtf_parsed %>% unnest_wider(attr_list)
        saveRDS(object = gtf_parsed, file=gtf_parsed_file)
        
    }
    
    return(gtf_parsed)
}


################
## UMAP display
################

get_isoform_umap = function(gene_of_interest, restrict_to_transcript_ids = NULL) {
    
    
    if (! is.null(restrict_to_transcript_ids)) {
        transcript_expr_data = data.frame(isoform_expr_data[
            rownames(isoform_expr_data) %in% restrict_to_transcript_ids,])
    } else {
        
        transcript_expr_data = data.frame(isoform_expr_data[grepl(gene_of_interest, rownames(isoform_expr_data)),])
    }
    
    transcript_expr_data$transcript_id = rownames(transcript_expr_data)
    transcript_expr_data = transcript_expr_data %>% gather(key=cell_barcode, value=read_count, -transcript_id)
    
    umap_df_w_expr_data = right_join(umap_df, transcript_expr_data,
                                     by='cell_barcode')
    
    return(umap_df_w_expr_data)
}


plot_isoform_umap = function(gene_of_interest, restrict_to_transcript_ids = NULL) {
    
    isoform_umap = get_isoform_umap(gene_of_interest, restrict_to_transcript_ids)
    
    #isoform_umap =  isoform_umap %>% filter(read_count > 0) %>% mutate(log_read_count = log1p(read_count)) 
    
    
    isoform_umap =  isoform_umap %>% filter(read_count > 0) 
    
    # Calculate 95th percentile for color scale
    max_color_value = quantile(isoform_umap$read_count, 0.95, na.rm = TRUE)
    
    base_umap + geom_point(data=isoform_umap, aes(color=read_count)) +
        facet_wrap(~transcript_id) +
        #theme_bw() +
        ggtitle(gene_of_interest)  +
        theme(legend.position="none") +
        scale_color_viridis_c(limits = c(0, max_color_value), oob = scales::squish) +
        geom_text(data = umap_df %>% 
                      group_by(seurat_clusters) %>% 
                      summarise(umap_1 = mean(umap_1), 
                                umap_2 = mean(umap_2)), 
                  aes(label = seurat_clusters), 
                  size = 5,
                  color = 'purple',
                  fontface = "bold") +
    theme_void()
    
}


#####################################
# Per-cell isoform USAGE FRACTION on the UMAP
#####################################

# Transcript ids belonging to a gene. Exact gene_id match first, then an ANCHORED
# pattern so that asking for "EIF1" cannot pull in "EIF1AX" / "EIF1B".
get_gene_transcript_ids = function(gene_of_interest) {
    
    ids = gtf_parsed %>%
        filter(gene_id == gene_of_interest, feature == "transcript") %>%
        pull(transcript_id) %>% unique()
    
    if (length(ids) == 0) {
        ids = gtf_parsed %>%
            filter(grepl(paste0("^", gene_of_interest, "($|[^A-Za-z0-9_])"), gene_id),
                   feature == "transcript") %>%
            pull(transcript_id) %>% unique()
    }
    
    return(ids)
}


# Load the cell-cell neighbor graph that Seurat built alongside the clustering and UMAP,
# so smoothing uses the SAME neighborhood structure the embedding came from rather than
# a fresh kNN computed on 2-D UMAP coordinates (which are a lossy projection).
#   "RNA_snn"  shared-nearest-neighbor, edge weight = neighbor overlap, ~68 neighbors/cell
#   "RNA_nn"   plain kNN, binary, k = 20
load_cell_neighbor_graph = function(seurat_obj_rds, graph_name = "RNA_snn") {
    
    message("-loading neighbor graph ", graph_name, " from ", seurat_obj_rds)
    obj = readRDS(seurat_obj_rds)
    
    if (! graph_name %in% Graphs(obj)) {
        stop("graph '", graph_name, "' not in object; available: ",
             paste(Graphs(obj), collapse=", "))
    }
    
    G = obj[[graph_name]]
    cell_neighbor_graph <<- G
    
    message("  ", nrow(G), " cells, median ", median(diff(G@p)), " neighbors/cell")
    
    invisible(G)
}


# Per-CELL isoform usage fraction: reads on an isoform divided by reads on the
# denominator set in that same cell.
#
#   denominator = "gene"     all isoforms of the gene. Matches how pi (and therefore
#                            delta_pi in the DTU results) is defined when the DTU test
#                            is run with --group_by_feature gene_symbol. Default.
#   denominator = "selected" only the transcripts being plotted, so a pair's fractions
#                            are complementary and sum to 1. Use when the question is
#                            strictly "which of these two", e.g. alt-termini pairs.
#
# min_reads_per_cell is a floor on the cell's OWN denominator reads for the unsmoothed
# view. It defaults to 0, meaning any cell with nonzero reads: cells with no reads carry
# no fraction and are always excluded (a 0/0 has no value), but nothing above that is
# required.
#
# A floor of 2 reads was the earlier default and it cost far too much. It discarded 86%
# of TARDBP's expressing cells (227 of 1,676), 66% of STMN2's and 27% of GPM6A's, and it
# silently dropped every cell whose EM assignment is a FRACTION below the floor -- 460
# cells for TARDBP alone sit in (0,1). Those cells hold real evidence.
#
# The cost of 0 is granularity, not correctness: a cell holding one read can only report
# a fraction of exactly 0 or 1, so the unsmoothed panel is coarse where coverage is thin
# (TARDBP: 48% at exactly 0, 20% at exactly 1). That is what the per-cell data actually
# says, and the smoothed panel beneath it is what resolves it. Raise the floor if you
# want the raw panel restricted to better-covered cells, at the price of hiding most of
# the expressing ones.
#
# 0 also puts the unsmoothed panel on the SAME cell set as the smoothed panel under
# cell_evidence = "gene", so the two differ in their values rather than in which cells
# they draw.
# smooth_graph: NULL for the raw per-cell view, or a cell-cell graph (see
# load_cell_neighbor_graph) to pool each cell with its neighbors.
#
# Smoothing POOLS COUNTS over the neighborhood and divides once:
#
#     sum_n w_cn * reads_isoform(n)  /  sum_n w_cn * reads_denominator(n)
#
# rather than averaging the neighbors' individual usage fractions. The distinction
# matters: averaging fractions gives a neighbor holding one read the same say as one
# holding fifty, which reintroduces exactly the 0/1 noise min_reads_per_cell exists to
# suppress. Pooling weights every neighbor by its depth automatically, is defined for a
# cell with no reads of its own, and is the same quantity as the cluster-level pi -- just
# measured over a neighborhood instead of a cluster. It is also invariant to row-scaling
# of the graph, so normalized and unnormalized weights give identical answers.
# cell_evidence controls WHICH cells may receive a smoothed value. Smoothing can either
# adjust cells that carry their own evidence, or paint cells that carry none at all --
# these are different claims and the choice should be explicit.
#
#   "gene"    (default) only cells with > 0 reads of the denominator in that cell. The
#             neighborhood refines a measurement the cell actually contributed to.
#   "any"     every cell whose NEIGHBORHOOD clears min_neighborhood_reads, including
#             cells with no reads of the gene at all. Smoothest field, but much of it is
#             pure imputation: for EIF1 that is 3,454 of 6,524 painted cells (53%) with
#             zero reads of the gene. Legitimate for showing where usage would be
#             expected, not for claiming the gene was observed there.
#   "isoform" only cells with > 0 reads of THAT isoform. Looks stricter and is actually
#             BIASED: it drops precisely the cells whose observed usage of the isoform
#             is zero, so the surviving field is shifted upward. Provided for
#             completeness; prefer "gene".
#
# min_neighborhood_reads interacts with cell_evidence, and the two must be set together:
#
#   With cell_evidence = "any" it is the ONLY guard against a value computed from almost
#   nothing, so it should be positive.
#
#   With cell_evidence = "gene" it mostly deletes real cells. Every displayed cell
#   already holds its own reads, and a cell in a thinly-covered region is sparse biology,
#   not noise. Hence the default of 0. A floor of 10 discarded 433 of STMN2's 713
#   expressing cells (median 1 own read, median neighborhood total 5) and 1,375 of EIF1's
#   3,070, while changing the spread of the retained values hardly at all -- EIF1 sd
#   0.134 at a floor of 10 against 0.145 at 0. It was gating display, not improving
#   estimates.
get_isoform_usage_fraction_umap = function(gene_of_interest,
                                           transcript_ids,
                                           denominator = c("gene", "selected"),
                                           min_reads_per_cell = 0,
                                           smooth_graph = NULL,
                                           min_neighborhood_reads = 0,
                                           weighted = TRUE,
                                           cell_evidence = c("gene", "any", "isoform"),
                                           ignore_unspliced = FALSE) {
    
    denominator   = match.arg(denominator)
    cell_evidence = match.arg(cell_evidence)
    
    denom_ids = if (denominator == "gene") get_gene_transcript_ids(gene_of_interest) else transcript_ids
    
    # Same convention as get_expression_ggplot2_heatmap_w_exon_structures: unspliced
    # models are identified by the ":iso-" naming. This matters on a splice-pattern
    # matrix, where the unspliced entries are a separate slice of the gene's signal
    # (~15% for STMN2) and the DTU test that produced delta_pi ran with
    # --ignore_unspliced. Leaving them in the denominator makes these fractions
    # disagree with the delta_pi values shown alongside them.
    if (ignore_unspliced && denominator == "gene") {
        denom_ids = denom_ids[! grepl(":iso-", denom_ids)]
    }
    
    denom_ids = union(denom_ids, transcript_ids)   # numerator must be inside the denominator
    
    present = rownames(isoform_expr_data) %in% denom_ids
    if (! any(present)) {
        stop("No transcripts of ", gene_of_interest, " found in the sparse matrix.")
    }
    
    mat = isoform_expr_data[present, , drop = FALSE]     # drop=FALSE: a single isoform must stay a matrix
    
    if (! is.null(smooth_graph)) {
        
        # Align cells: the graph and the isoform matrix must be on the same barcodes,
        # in the same order, before any matrix product.
        common = intersect(colnames(mat), rownames(smooth_graph))
        if (length(common) == 0) {
            stop("No cell barcodes shared between the isoform matrix and the graph.")
        }
        if (length(common) < ncol(mat)) {
            message("  smoothing over ", length(common), " of ", ncol(mat),
                    " cells present in both the matrix and the graph")
        }
        
        G = smooth_graph[common, common, drop = FALSE]
        if (! weighted) {
            G@x = rep(1, length(G@x))                   # binary neighborhood
        }
        
        mat_c = mat[, common, drop = FALSE]
        den_vec = Matrix::colSums(mat_c)                # denominator reads per cell
        
        sm_den = as.vector(G %*% den_vec)
        
        sm = lapply(transcript_ids, function(tid) {
            if (! tid %in% rownames(mat_c)) return(NULL)
            own_iso = as.vector(mat_c[tid, ])
            sm_num  = as.vector(G %*% own_iso)
            data.frame(transcript_id = tid,
                       cell_barcode  = common,
                       read_count    = own_iso,
                       own_denominator_reads = den_vec,
                       neighborhood_reads = sm_den,
                       usage_fraction = ifelse(sm_den > 0, sm_num / sm_den, NA_real_),
                       stringsAsFactors = FALSE)
        })
        
        usage_df = bind_rows(sm) %>%
            filter(!is.na(usage_fraction), neighborhood_reads >= min_neighborhood_reads)
        
        n_before = n_distinct(usage_df$cell_barcode)
        
        usage_df = switch(cell_evidence,
            any     = usage_df,
            gene    = usage_df %>% filter(own_denominator_reads > 0),
            isoform = usage_df %>% filter(read_count > 0))
        
        n_after = n_distinct(usage_df$cell_barcode)
        if (n_after < n_before) {
            message("  cell_evidence='", cell_evidence, "': ", n_before - n_after, " of ",
                    n_before, " cells dropped for carrying no qualifying reads of their own")
        }
        
        usage_df = usage_df %>% inner_join(umap_df, by = "cell_barcode")
        
        return(usage_df)
    }
    
    expr = as.data.frame(as.matrix(mat), check.names = FALSE)
    expr$transcript_id = rownames(expr)
    
    expr_long = expr %>%
        gather(key = "cell_barcode", value = "read_count", -transcript_id)
    
    per_cell = expr_long %>%
        group_by(cell_barcode) %>%
        summarise(denominator_reads = sum(read_count), .groups = "drop")
    
    usage_df = expr_long %>%
        filter(transcript_id %in% transcript_ids) %>%
        left_join(per_cell, by = "cell_barcode") %>%
        filter(denominator_reads > 0, denominator_reads >= min_reads_per_cell) %>%
        mutate(usage_fraction = read_count / denominator_reads) %>%
        inner_join(umap_df, by = "cell_barcode")
    
    return(usage_df)
}


# One UMAP panel per transcript, coloured by per-cell usage fraction on a shared 0..1
# scale so the panels are directly comparable.
plot_isoform_usage_fraction_umap = function(gene_of_interest,
                                            transcript_ids,
                                            denominator = c("gene", "selected"),
                                            min_reads_per_cell = 0,
                                            point_size = 0.8,
                                            label_clusters = TRUE,
                                            ncol = NULL,
                                            smooth_graph = NULL,
                                            min_neighborhood_reads = 0,
                                            weighted = TRUE,
                                            cell_evidence = c("gene", "any", "isoform"),
                                            ignore_unspliced = FALSE,
                                            transcript_labels = NULL) {
    
    denominator   = match.arg(denominator)
    cell_evidence = match.arg(cell_evidence)
    
    usage_df = get_isoform_usage_fraction_umap(gene_of_interest, transcript_ids,
                                               denominator, min_reads_per_cell,
                                               smooth_graph, min_neighborhood_reads,
                                               weighted, cell_evidence, ignore_unspliced)
    
    n_cells = usage_df %>% distinct(cell_barcode) %>% nrow()
    if (n_cells == 0) {
        stop("No cells pass the read floor for ", gene_of_interest, "; nothing to plot.")
    }
    
    # Draw low usage first so high-usage cells are not hidden under overplotting.
    usage_df = usage_df %>% arrange(usage_fraction)

    # transcript_labels: named by transcript id; relabels the panels and orders them as given
    if (! is.null(transcript_labels)) {
        usage_df = usage_df %>%
            mutate(transcript_id = factor(transcript_labels[transcript_id], levels = unname(transcript_labels)))
    }
    
    denom_label = if (denominator == "gene") "all isoforms of the gene" else "the plotted isoforms"
    
    evidence_label = switch(cell_evidence,
        gene    = "cells with their own reads of the gene",
        any     = "all covered cells, INCLUDING cells with no reads of the gene",
        isoform = "cells with their own reads of that isoform (upward biased)")
    
    subtitle = if (is.null(smooth_graph)) {
        if (min_reads_per_cell > 0) {
            paste0(n_cells, " cells with >= ", min_reads_per_cell, " reads over ", denom_label)
        } else {
            paste0(n_cells, " cells with any reads over ", denom_label)
        }
    } else {
        paste0("SNN-smoothed: ", n_cells, " ", evidence_label,
               "; neighborhood >= ", min_neighborhood_reads, " reads over ", denom_label,
               if (weighted) " (edge-weighted)" else " (unweighted)")
    }
    
    p = base_umap +
        geom_point(data = usage_df, aes(color = usage_fraction), size = point_size) +
        facet_wrap(~ transcript_id, ncol = ncol) +
        scale_color_viridis_c(limits = c(0, 1), name = "isoform\nusage fraction") +
        ggtitle(gene_of_interest, subtitle = subtitle)
    
    if (label_clusters) {
        p = p + geom_text(data = umap_df %>%
                              group_by(seurat_clusters) %>%
                              summarise(umap_1 = mean(umap_1),
                                        umap_2 = mean(umap_2), .groups = "drop"),
                          aes(label = seurat_clusters),
                          size = 4, color = 'purple', fontface = "bold")
    }
    
    p = p + theme_void()
    
    return(p)
}


# How much two transcripts occupy the same locus, as
#   span   : overlap of their outermost coordinates
#   exonic : shared EXONIC bases
# both expressed as a fraction of the SHORTER transcript, so a short isoform nested
# inside a long one scores 1 rather than being penalised for its length.
#
# Both numbers are needed. Span overlap alone accepts a transcript sitting inside
# another's INTRON -- real case in this data: GPM6A's c9c63b4f nests entirely within
# 92af509f's span (span 1.00) while sharing no exonic base with it (exonic 0.00).
# Returns NA when either transcript is absent from gtf_parsed.
transcript_pair_overlap = function(tx_a, tx_b, exons = NULL) {
    
    if (is.null(exons)) {
        exons = gtf_parsed %>% filter(feature == "exon") %>%
            transmute(transcript_id, seqname,
                      start = as.integer(start), end = as.integer(end))
    }
    
    A = exons %>% filter(transcript_id == tx_a)
    B = exons %>% filter(transcript_id == tx_b)
    if (nrow(A) == 0 || nrow(B) == 0) return(c(span = NA_real_, exonic = NA_real_))
    if (A$seqname[1] != B$seqname[1])  return(c(span = 0, exonic = 0))
    
    span_a = c(min(A$start), max(A$end)); span_b = c(min(B$start), max(B$end))
    span_hit = max(0, min(span_a[2], span_b[2]) - max(span_a[1], span_b[1]) + 1)
    span_frac = span_hit / min(diff(span_a) + 1, diff(span_b) + 1)
    
    # interval intersection, not coordinate expansion: exon counts are small and the
    # spans here reach 370 kb.
    exonic_hit = 0
    for (i in seq_len(nrow(A))) {
        exonic_hit = exonic_hit +
            sum(pmax(0, pmin(A$end[i], B$end) - pmax(A$start[i], B$start) + 1))
    }
    len_a = sum(A$end - A$start + 1); len_b = sum(B$end - B$start + 1)
    exonic_frac = exonic_hit / min(len_a, len_b)
    
    return(c(span = span_frac, exonic = exonic_frac))
}


# The two transcripts showing the strongest differential usage for a gene, taken from a
# SINGLE DTU row so the pair is one measured comparison rather than two ids collected
# across different cluster pairs.
#
# by = "reciprocal" (default) ranks on min(|delta_pi|, |alternate_delta_pi|). The DTU
#      test here runs with --reciprocal_delta_pi, so a pair is only convincingly trading
#      share when BOTH transcripts move; the weaker of the two shifts is the honest
#      measure of that. A large delta_pi paired with a negligible reciprocal is usually
#      one transcript moving against the rest of the locus, not a swap.
# by = "delta_pi" ranks on |delta_pi| alone, i.e. the dominant transcript's shift.
# by = "pvalue"   ranks on significance instead of effect size.
#
# Rows failing `significant` are excluded when that column exists; if a gene has none,
# all its rows are ranked and a message says so, because an insignificant top pair is a
# legitimate thing to look at and a silent empty plot is not.
#
# The pair must also share EXONIC sequence, so that the two transcripts are structural
# alternatives at one locus rather than separate transcription units that a shared gene
# symbol happens to group. Because the DTU test groups by gene_symbol, nothing upstream
# enforces this: GPM6A's strongest-scoring pair is 167 kb apart with no shared base, and
# a usage fraction between two disjoint units is not an isoform choice at all.
#
# min_exonic_overlap is the gate that matters and is the only one on by default. Nonzero
# exonic overlap implies the spans overlap, so it subsumes a span test while also
# rejecting a transcript nested inside another's INTRON (GPM6A c9c63b4f inside 92af509f:
# span 1.00, exonic 0.00). Keep it low. It is a "same locus?" test, not a similarity
# test, and staggered alternatives legitimately share little: TARDBP's pair shares 508
# exonic bases for a fraction of 0.38, and an earlier span-based bar of 0.5 wrongly
# rejected it at 0.31.
#
# min_span_overlap defaults to 0, i.e. off. Raise it to demand roughly co-extensive
# transcripts, which is a narrower question than sharing a locus.
#
# If no pair clears the bars, the best-scoring pair is returned with a WARNING rather
# than an error: a fragmented locus is worth seeing, but it must not pass silently.
get_top_dtu_pair = function(gene_of_interest,
                            stats = diff_iso_usage_stats,
                            by = c("reciprocal", "delta_pi", "pvalue"),
                            require_significant = TRUE,
                            min_span_overlap = 0,
                            min_exonic_overlap = 0.05) {
    
    by = match.arg(by)
    
    rows = stats %>% filter(gene_symbol == gene_of_interest)
    if (nrow(rows) == 0) {
        stop("No DTU rows for ", gene_of_interest)
    }
    
    if (require_significant && "significant" %in% colnames(rows)) {
        sig = rows %>% filter(significant == "True")
        if (nrow(sig) > 0) {
            rows = sig
        } else {
            message("  ", gene_of_interest, ": no significant DTU rows; ranking all ",
                    nrow(rows), " comparisons")
        }
    }
    
    # A side can carry TWO transcripts, comma-joined. The DTU code sums the top movers
    # per direction over a hard-coded [:2] slice, which --top_isoforms_each does not
    # constrain (it limits the candidate pool, not the slice). 168 of the 21,197 rows in
    # this dataset are like that, across 68 genes with significant rows, GPM6A included.
    # Such a row has no single pair to plot, and passing the joined string on would look
    # up a transcript id that cannot exist, so drop those rows rather than emit a panel
    # for a feature absent from the matrix.
    n_multi = sum(grepl(",", rows$dominant_transcript_ids) |
                  grepl(",", rows$alternate_transcript_ids))
    if (n_multi > 0) {
        message("  ", gene_of_interest, ": ", n_multi, " of ", nrow(rows),
                " comparisons dropped for carrying multiple transcripts on one side")
        rows = rows %>% filter(! grepl(",", dominant_transcript_ids),
                               ! grepl(",", alternate_transcript_ids))
        if (nrow(rows) == 0) {
            stop(gene_of_interest, ": every comparison pairs transcript SETS rather than ",
                 "single transcripts; no pair to plot")
        }
    }
    
    rows = rows %>% mutate(.dtu_score = switch(by,
        # equals abs(alternate_delta_pi): dominant is by definition the larger-magnitude
        # side, verified on all 21,197 rows. Written as pmin to stay correct if that
        # convention ever changes upstream.
        reciprocal = pmin(abs(delta_pi), abs(alternate_delta_pi)),
        delta_pi   = abs(delta_pi),
        pvalue     = -pvalue)) %>%
        arrange(desc(.dtu_score), pvalue)
    
    keep = rows
    if (min_span_overlap > 0 || min_exonic_overlap > 0) {
        
        exons = gtf_parsed %>% filter(feature == "exon") %>%
            transmute(transcript_id, seqname,
                      start = as.integer(start), end = as.integer(end))
        
        # one overlap computation per distinct pair, not per row
        pairs = rows %>% distinct(dominant_transcript_ids, alternate_transcript_ids)
        ov = t(mapply(transcript_pair_overlap,
                      pairs$dominant_transcript_ids, pairs$alternate_transcript_ids,
                      MoreArgs = list(exons = exons)))
        pairs$.span = ov[, "span"]; pairs$.exonic = ov[, "exonic"]
        
        keep = rows %>%
            left_join(pairs, by = c("dominant_transcript_ids", "alternate_transcript_ids")) %>%
            filter(!is.na(.span), !is.na(.exonic),
                   .span >= min_span_overlap, .exonic >= min_exonic_overlap)
        
        n_drop = nrow(rows) - nrow(keep)
        if (n_drop > 0) {
            message("  ", gene_of_interest, ": ", n_drop, " of ", nrow(rows),
                    " comparisons dropped for insufficient genomic overlap")
        }
        
        if (nrow(keep) == 0) {
            warning(gene_of_interest, ": no DTU pair clears the overlap bars (span >= ",
                    min_span_overlap, ", exonic >= ", min_exonic_overlap,
                    "); returning the top-scoring pair, which does NOT share the locus.",
                    call. = FALSE)
            keep = rows
        }
    }
    
    top = keep %>% slice(1)
    
    message("  ", gene_of_interest, " top pair: ", top$cluster_A, " vs ", top$cluster_B,
            "  delta_pi=", round(top$delta_pi, 3),
            "  alternate_delta_pi=", round(top$alternate_delta_pi, 3),
            "  p=", signif(top$pvalue, 3),
            if (!is.null(top$.span)) paste0("  span_overlap=", round(top$.span, 2),
                                            "  exonic_overlap=", round(top$.exonic, 2)) else "")
    
    return(c(top$dominant_transcript_ids, top$alternate_transcript_ids))
}


# Per-cell usage fraction for a gene's strongest DTU pair, as the raw view and -- when a
# graph is supplied -- the SNN-smoothed view stacked beneath it, so the two are read
# against each other on the same colour scale.
plot_top_dtu_pair_usage_umaps = function(gene_of_interest,
                                         stats = diff_iso_usage_stats,
                                         smooth_graph = NULL,
                                         denominator = c("gene", "selected"),
                                         ignore_unspliced = FALSE,
                                         min_reads_per_cell = 0,
                                         cell_evidence = c("gene", "any", "isoform"),
                                         by = c("reciprocal", "delta_pi", "pvalue"),
                                         require_significant = TRUE,
                                         min_span_overlap = 0,
                                         min_exonic_overlap = 0.05,
                                         point_size = 0.8) {
    
    denominator   = match.arg(denominator)
    cell_evidence = match.arg(cell_evidence)
    by            = match.arg(by)
    
    pair = get_top_dtu_pair(gene_of_interest, stats, by, require_significant,
                            min_span_overlap, min_exonic_overlap)
    
    p_raw = plot_isoform_usage_fraction_umap(gene_of_interest, pair,
                                             denominator = denominator,
                                             ignore_unspliced = ignore_unspliced,
                                             min_reads_per_cell = min_reads_per_cell,
                                             point_size = point_size)
    
    if (is.null(smooth_graph)) {
        return(p_raw)
    }
    
    p_smooth = plot_isoform_usage_fraction_umap(gene_of_interest, pair,
                                                denominator = denominator,
                                                ignore_unspliced = ignore_unspliced,
                                                smooth_graph = smooth_graph,
                                                cell_evidence = cell_evidence,
                                                point_size = point_size)
    
    return(plot_grid(p_raw, p_smooth, ncol = 1))
}

#####################################
# Gene structure and heatmap display
#####################################

get_gene_structure_matrix = function(gene_of_interest) {
    # Filter to exons only - prioritize exact match, fall back to pattern match
    exon_df <- gtf_parsed %>%
        filter(gene_id == gene_of_interest) %>%
        filter(feature == "exon")
    
    # If no exact match found, try pattern matching
    if (nrow(exon_df) == 0) {
        exon_df <- gtf_parsed %>%
            filter(grepl(paste0("^", gene_of_interest, "($|[^A-Za-z0-9_])"), gene_id)) %>%
            filter(feature == "exon")
    }
    
    exon_df <- exon_df %>%
        mutate(
            start = as.integer(start),
            end = as.integer(end),
            exon_coords = paste0(start, "-", end)
        )
    
    # Sort exon coordinates by start position
    exon_levels <- exon_df %>%
        distinct(exon_coords, start) %>%
        arrange(start) %>%
        pull(exon_coords)
    
    # Build binary presence matrix (long format)
    exon_binary_df <- exon_df %>%
        select(transcript_id, exon_coords) %>%
        distinct() %>%
        mutate(present = 1)
    
    # Set factor levels for exon coordinates (x-axis)
    exon_binary_df$exon_coords <- factor(exon_binary_df$exon_coords, levels = exon_levels)
    
    
    # Ensure binary_df is a regular data frame (not grouped or list-columned)
    exon_binary_df <- exon_binary_df %>% ungroup()
    
    # Count number of exons per transcript using tally() and sort by num exons
    transcript_levels <- exon_binary_df %>%
        filter(present == 1) %>%
        group_by(transcript_id) %>%
        tally(name = "n_exons") %>%
        arrange(desc(n_exons)) %>%
        pull(transcript_id)
    
    # Set transcript_id factor levels accordingly
    exon_binary_df$transcript_id <- factor(exon_binary_df$transcript_id, levels = transcript_levels)
    
    return(exon_binary_df)
    
}

gene_structure_heatmap = function(gene_of_interest) {
    
    exon_binary_df = get_gene_structure_matrix(gene_of_interest)
    
    # Plot heatmap
    p = ggplot(exon_binary_df, aes(x = exon_coords, y = transcript_id, fill = factor(present))) +
        geom_tile(color = "grey80") +
        scale_fill_manual(values = c("1" = "black"), guide = "none") +
        labs(
            title = "Transcript Exon Usage Heatmap",
            x = "Exon Coordinates (start-end)",
            y = "Transcript ID"
        ) +
        theme_minimal(base_size = 10) +
        theme(
            axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
            panel.grid = element_blank()
        )
    
    return(p)
    
}


get_expression_pheatmap_w_exon_structures = function(gene_of_interest) {
    
    # Prioritize exact match, fall back to pattern match
    transcript_ids = gtf_parsed %>% filter(gene_id == gene_of_interest) %>% 
        filter(feature == "transcript") %>% select(transcript_id) %>% unique() %>% pull(transcript_id)
    
    # If no exact match found, try pattern matching
    if (length(transcript_ids) == 0) {
        transcript_ids = gtf_parsed %>% filter(grepl(paste0("^", gene_of_interest, "($|[^A-Za-z0-9_])"), gene_id)) %>% 
            filter(feature == "transcript") %>% select(transcript_id) %>% unique() %>% pull(transcript_id)
    }
    
    isoform_expr = cluster_CPM_matrix[rownames(cluster_CPM_matrix) %in% transcript_ids,]
    
    # Load expression matrix
    expr_mat <- as.matrix(isoform_expr)
    
    # Log transform if desired
    log_expr <- log2(expr_mat + 1)
    
    # Load exon structure info (from earlier script)
    # binary_df has transcript_id, exon_coords, and present (1/0)
    # Ensure transcript IDs match (e.g., same format)
    
    exon_binary_df = get_gene_structure_matrix(gene_of_interest)
    
    # Create exon structure annotation matrix (transcripts x exon_coords)
    exon_mat <- exon_binary_df %>%
        filter(present == 1) %>%
        select(transcript_id, exon_coords) %>%
        mutate(value = "●") %>%
        pivot_wider(names_from = exon_coords, values_from = value, values_fill = "") %>%
        column_to_rownames("transcript_id") %>%
        as.data.frame()
    
    
    
    # Get exon column names and sort by numeric start coordinate
    sorted_exon_cols <- exon_mat %>%
        colnames() %>%
        as_tibble() %>%
        dplyr::rename(coord = value) %>%
        mutate(start = as.integer(str_extract(coord, "^\\d+"))) %>%
        arrange(start) %>%
        pull(coord)
    
    # Reorder the exon matrix columns
    exon_mat_sorted <- exon_mat[, rev(sorted_exon_cols)]
    
    
    # If exon presence is encoded as "●" and "" or 1/0
    # Convert to "yes"/"no" to be explicit
    exon_mat_color <- exon_mat_sorted %>%
        mutate(across(everything(), ~ ifelse(. != "", "yes", "no")))  # or . == 1
    
    # Make all columns factors
    exon_mat_color <- exon_mat_color %>%
        mutate(across(everything(), as.factor))
    
    # Build color mapping for all exon annotation columns
    annotation_colors <- list()
    
    for (col in colnames(exon_mat_color)) {
        annotation_colors[[col]] <- c("yes" = "black", "no" = "white")
    }
    
    # Filter so exon_mat rows match log_expr
    common_ids <- intersect(rownames(log_expr), rownames(exon_mat))
    log_expr <- log_expr[common_ids, ]
    exon_mat_sorted <- exon_mat_sorted[common_ids, ]
    
    
    # Generate heatmap with exon structure annotations as row labels
    pheatmap(log_expr,
             show_rownames = TRUE,
             show_colnames = TRUE,
             cluster_rows = TRUE,
             cluster_cols = TRUE,
             annotation_row = exon_mat_color,
             annotation_colors = annotation_colors,
             fontsize_row = 6,
             main = "Transcript Expression with Exon Structure",
             annotation_legend = FALSE,
    )
    
}


get_expression_ggplot2_heatmap_w_exon_structures = function(
        gene_of_interest, 
        min_isoform_frac_expr_any_cluster = 0,
        ignore_unspliced = FALSE,
        transcript_ids = NULL,
        min_cells_expressed = 0,
        min_cell_frac_expressed = 0) {
    
    library(patchwork)
    library(ggdendro)
    
    message("Transcript_ids: ", transcript_ids)

    # Prioritize exact match, fall back to pattern match
    all_transcript_ids_for_gene = gtf_parsed %>% filter(gene_id == gene_of_interest) %>% 
            filter(feature == "transcript") %>% select(transcript_id) %>% unique() %>% pull(transcript_id)
    
    # If no exact match found, try pattern matching
    if (length(all_transcript_ids_for_gene) == 0) {
        all_transcript_ids_for_gene = gtf_parsed %>% filter(grepl(paste0("^", gene_of_interest, "($|[^A-Za-z0-9_])"), gene_id)) %>% 
                filter(feature == "transcript") %>% select(transcript_id) %>% unique() %>% pull(transcript_id)
    }

    all_expr_for_gene = cluster_CPM_matrix[rownames(cluster_CPM_matrix) %in% all_transcript_ids_for_gene, ]

    # Determine the transcript set to use for computing isoform fractions (denominator)
    fraction_basis_transcript_ids = all_transcript_ids_for_gene
    if (ignore_unspliced) {
        fraction_basis_transcript_ids = fraction_basis_transcript_ids[ ! grepl(":iso-", fraction_basis_transcript_ids)]
    }
    
    # Calculate isoform fractions based on all (or all spliced) transcripts
    fraction_basis_expr = all_expr_for_gene[rownames(all_expr_for_gene) %in% fraction_basis_transcript_ids, ]
    all_isoform_frac_expr <- sweep(fraction_basis_expr, 2, colSums(fraction_basis_expr), "/")
    all_isoform_frac_expr[is.na(all_isoform_frac_expr)] = 0
    
    
    if (is.null(transcript_ids)) {
        transcript_ids = fraction_basis_transcript_ids
    }
        
    message("Transcript_ids: ", transcript_ids)
    
    isoform_expr = all_expr_for_gene[rownames(all_expr_for_gene) %in% transcript_ids,]
    
    # Extract isoform fractions for the transcripts we're displaying
    isoform_frac_expr <- all_isoform_frac_expr[rownames(all_isoform_frac_expr) %in% transcript_ids, ]
    isoform_frac_expr[is.na(isoform_frac_expr)] = 0
    
    # Two INDEPENDENT criteria, because they ask different questions.
    #
    # min_cells_expressed / min_cell_frac_expressed is an ABUNDANCE test: is this isoform
    # actually detected in a decent number of cells somewhere? That is the right basis for
    # a DTU view, where the isoforms of interest are the ones TRADING share between
    # clusters. They need not ever dominate: STMN2's cryptic isoform sits in 53 of 299
    # cells in Cluster_8 and peaks at an isoform fraction of only 0.39.
    #
    # min_isoform_frac_expr_any_cluster is a DOMINANCE test: does this isoform win its
    # locus in some cluster? Useful for trimming a long tail, but it is depth-blind on its
    # own -- a cluster holding 2 assigned read-equivalents lets 1.0/2.0 = 0.50 outrank an
    # isoform carrying 62 in a deep cluster -- so when both are set, a cluster must clear
    # the cell bar before its isoform fraction is allowed to count.

    cell_ok = NULL
    if (min_cells_expressed > 0 || min_cell_frac_expressed > 0) {
        if (is.null(cluster_cell_fraction_matrix)) {
            stop("min_cells_expressed / min_cell_frac_expressed require parse_inputs() ",
                 "to have been given cell_fractions_matrix_filename")
        }
        cell_frac = cluster_cell_fraction_matrix[rownames(isoform_frac_expr),
                                                 colnames(isoform_frac_expr), drop=FALSE]
        cell_frac[is.na(cell_frac)] = 0
        n_cells = cluster_cell_counts[colnames(isoform_frac_expr)]
        cells_expressed = sweep(as.matrix(cell_frac), 2, n_cells, "*")
        # stricter of the absolute and the proportional bar, per cluster
        bar = pmax(min_cells_expressed,
                   matrix(rep(min_cell_frac_expressed * n_cells, each=nrow(cell_frac)),
                          nrow=nrow(cell_frac)))
        cell_ok = cells_expressed >= bar
    }

    if (min_isoform_frac_expr_any_cluster > 0 || ! is.null(cell_ok)) {

        qualifies = if (is.null(cell_ok)) {
            isoform_frac_expr >= min_isoform_frac_expr_any_cluster
        } else if (min_isoform_frac_expr_any_cluster > 0) {
            (isoform_frac_expr >= min_isoform_frac_expr_any_cluster) & cell_ok
        } else {
            cell_ok
        }

        transcript_ids = rownames(isoform_frac_expr)[rowSums(qualifies) > 0]
        
        if (length(transcript_ids) < 2) {
            stop("Too few isoforms left after filtering on isoform fraction / cell support") 
        }
        
        isoform_expr = isoform_expr[rownames(isoform_expr) %in% transcript_ids,]
        
        # Filter isoform_frac_expr to match (but don't recalculate - keep original denominators)
        isoform_frac_expr = isoform_frac_expr[rownames(isoform_frac_expr) %in% transcript_ids,]
        
    }
    
    
    # Load expression matrix
    expr_mat <- as.matrix(isoform_expr)
    
    # Log transform if desired
    log_expr <- log2(expr_mat + 1)
    
    
    col_order <- hclust(dist(t(log_expr)))$order
    ordered_clusters <- colnames(log_expr)[col_order]
    
    # Convert to long format
    expr_long <- as.data.frame(log_expr) %>%
        rownames_to_column("transcript_id") %>%
        pivot_longer(-transcript_id, names_to = "cluster", values_to = "expression")
    
    frac_expr_long = isoform_frac_expr %>%
        rownames_to_column("transcript_id") %>%
        pivot_longer(-transcript_id, names_to = "cluster", values_to = "iso_frac")
    
    
    expr_long$cluster <- factor(expr_long$cluster, levels = ordered_clusters)
    frac_expr_long$cluster <- factor(frac_expr_long$cluster, levels = ordered_clusters)
    
    
    # Cluster rows using base hclust on Euclidean distances
    row_order <- hclust(dist(log_expr))$order
    expr_long$transcript_id <- factor(expr_long$transcript_id, 
                                      levels = rownames(log_expr)[row_order])
    
    frac_expr_long$transcript_id <- factor(frac_expr_long$transcript_id, 
                                           levels = rownames(log_expr)[row_order])
    
    
    exon_binary_df = get_gene_structure_matrix(gene_of_interest)
    
    exon_binary_df = exon_binary_df %>% filter(transcript_id %in% transcript_ids)
    
    # Create exon structure annotation matrix (transcripts x exon_coords)
    exon_mat <- exon_binary_df %>%
        filter(present == 1) %>%
        select(transcript_id, exon_coords) %>%
        mutate(value = "●") %>%
        pivot_wider(names_from = exon_coords, values_from = value, values_fill = "") %>%
        column_to_rownames("transcript_id") %>%
        as.data.frame()
    
    
    
    # Long format for exon annotations
    exon_anno_long <- exon_mat %>%
        rownames_to_column("transcript_id") %>%
        pivot_longer(-transcript_id, names_to = "exon_coords", values_to = "present") %>%
        filter(present != "")  # keep only present exons
    
    # Ensure same transcript order as clustering
    exon_anno_long$transcript_id <- factor(exon_anno_long$transcript_id, 
                                           levels = levels(expr_long$transcript_id))
    
    
    
    
    # Get dendrogram data
    row_dend <- hclust(dist(log_expr))
    dend <- as.dendrogram(row_dend)
    ddata <- dendro_data(dend, type = "rectangle")
    
    # Get row order for consistency
    row_order <- row_dend$order
    transcript_order <- rownames(log_expr)[row_order]
    
    # Create mapping from x (numeric position) to label
    label_map <- ddata$labels %>%
        mutate(leaf_position = row_number()) %>%
        select(leaf_position, label)
    
    # Create a factor for labels ordered to match heatmap transcript_id
    label_map <- label_map %>%
        mutate(transcript_id = factor(label, levels = transcript_order),
               new_pos = as.numeric(transcript_id))
    
    # Replace x and xend in segments using leaf positions
    segments_fixed <- ddata$segments %>%
        left_join(label_map %>% select(leaf_position, new_pos), by = c("x" = "leaf_position")) %>%
        mutate(x = ifelse(!is.na(new_pos), new_pos, x)) %>%
        select(-new_pos) %>%
        left_join(label_map %>% select(leaf_position, new_pos), by = c("xend" = "leaf_position")) %>%
        mutate(xend = ifelse(!is.na(new_pos), new_pos, xend)) %>%
        select(-new_pos)
    
    p_dend <- ggplot(segments_fixed) +
        geom_segment(aes(x = y, xend = yend, y = x, yend = xend)) +
        scale_x_reverse() +  # ← This flips the tree so root is on the left
        theme_void() +
        theme(plot.margin = margin(t = 10, r = 5, b = 10, l = 10))
    
    
    # Expression heatmap
    p_expr <- ggplot(expr_long, aes(x = cluster, y = transcript_id, fill = expression)) +
        geom_tile() +
        scale_fill_viridis_c() +
        theme_minimal() +
        labs(title = "Iso Expression", x = "Cluster", y = NULL) +
        theme(axis.text.x = element_text(angle = 90, hjust = 1),
              axis.text.y = element_blank(),
              axis.ticks.y = element_blank())

    # Isoform fraction heatmap
    p_frac_expr = ggplot(frac_expr_long, aes(x = cluster, y = transcript_id, fill = iso_frac)) +
        geom_tile() +
        scale_fill_viridis_c() +
        #scale_fill_viridis_c(limits = c(0, 1), oob = scales::squish) +
        #scale_fill_gradient(
        #    low = "black", 
        #    high = "red", 
        #    limits = c(0, NA), 
        #    oob = scales::squish
        #) +
        #theme_minimal() +
        labs(title = "Iso Fraction", x = "Cluster", y = NULL) +
        theme(axis.text.x = element_text(angle = 90, hjust = 1),
              axis.text.y = element_blank(),
              axis.ticks.y = element_blank())

    
    # restore legends and use relative units for scaling
    p_expr <- p_expr + theme(
        legend.position = "bottom",
        legend.key.size = unit(1, "line"),
        legend.text = element_text(size = rel(0.8), angle = 45),
        legend.title = element_text(angle = 0)
    )
    p_frac_expr <- p_frac_expr + theme(
        legend.position = "bottom",
        legend.key.size = unit(1, "line"),
        legend.text = element_text(size = rel(0.8), angle = 45),
        legend.title = element_text(angle = 0)
    )
    
    
    # Exon structure tile map
    p_exon <- ggplot(exon_anno_long, aes(x = exon_coords, y = transcript_id)) +
        geom_tile(fill = "black") +
        theme_minimal() +
        labs(title = "Exon Structure", x = "Exon (start-end)", y = NULL) +
        theme(axis.text.x = element_text(angle = 90, hjust = 1),
              axis.text.y = element_text(size = 6),
              plot.margin = margin(t = 10, r = 10, b = 10, l = 0))
    
    
    
    g_dend <- ggplotGrob(p_dend)
    g_exon <- ggplotGrob(p_exon)
    g_expr <- ggplotGrob(p_expr)
    g_frac_expr <- ggplotGrob(p_frac_expr)
    
    
    common_heights <- grid::unit.pmax(g_dend$heights,
                                      g_exon$heights,
                                      g_expr$heights,
                                      g_frac_expr$heights)
    
    g_dend$heights <- common_heights
    g_exon$heights <- common_heights
    g_expr$heights <- common_heights
    g_frac_expr$heights <- common_heights
    
    
    
    combined_plot <- 
        wrap_elements(g_dend) +
        wrap_elements(g_exon) +
        wrap_elements(g_expr) +
        wrap_elements(g_frac_expr) +
        plot_layout(ncol = 4, widths = c(0.5, 1.2, 0.5, 0.5))
    
    return(list(
        plot = combined_plot, 
        transcript_ids = transcript_ids)
    )
    
}



######
# combine gene structure and iso expr / freac heatmaps together with umaps

library(cowplot)

# usage_fraction_umap = TRUE swaps the read-count umap for per-cell usage fractions of the
# same transcripts (see plot_isoform_usage_fraction_umap), smoothed when smooth_graph is given.
make_diff_iso_usage_compound_plot = function(gene_of_interest, min_iso_fraction = 0, ignore_unspliced=FALSE,
                                             transcript_ids=NULL,
                                             min_cells_expressed = 0,
                                             min_cell_frac_expressed = 0,
                                             usage_fraction_umap = FALSE,
                                             denominator = c("selected", "gene"),
                                             smooth_graph = NULL) {
    
    p_exon_expr_info = get_expression_ggplot2_heatmap_w_exon_structures(gene_of_interest, min_iso_fraction, 
                                                                        ignore_unspliced, transcript_ids,
                                                                        min_cells_expressed,
                                                                        min_cell_frac_expressed)
    
    p_exon_expr = p_exon_expr_info$plot
    
    if (usage_fraction_umap) {
        p_umap = plot_isoform_usage_fraction_umap(gene_of_interest, p_exon_expr_info$transcript_ids,
                                                  denominator = match.arg(denominator),
                                                  smooth_graph = smooth_graph,
                                                  cell_evidence = "gene",
                                                  ignore_unspliced = ignore_unspliced)
    } else {
        p_umap = plot_isoform_umap(gene_of_interest, p_exon_expr_info$transcript_ids)
    }
    
    p_both = plot_grid(p_exon_expr, p_umap, ncol=1)
    
    return(p_both)
    
}


################
# get gene-level expr view of umap

get_gene_umap = function(gene_symbol_val, gene_component) {
  
  gene_isoform_umap =  get_isoform_umap(gene_symbol_val)
  
  gene_isoform_umap = gene_isoform_umap %>% filter(grepl(gene_component, transcript_id)) %>% filter(read_count > 0) 
  
  gene_isoform_umap = gene_isoform_umap %>% group_by(cell_barcode) %>% mutate(read_count = sum(read_count))
  
  # Calculate 95th percentile for color scale
  max_color_value = quantile(gene_isoform_umap$read_count, 0.95, na.rm = TRUE)
  
  p_gene = base_umap + geom_point(data=gene_isoform_umap, aes(color=read_count)) +
    scale_color_viridis_c(limits = c(0, max_color_value), oob = scales::squish) +
    theme_void()

  return(p_gene)
}


###################3
# plot cluster pairs where transcripts show significant DTU, plotting delta_pi according to cell clusters

plot_dtu_pair_heatmap <- function(DTU_results, tx_dom, tx_alt) {
  
  # 1) Subset to the transcript pair in either orientation
  pair_df <- DTU_results %>%
    filter(
      (dominant_transcript_ids == tx_dom & alternate_transcript_ids == tx_alt) |
      (dominant_transcript_ids == tx_alt & alternate_transcript_ids == tx_dom)
    ) %>%
    select(
      dominant_transcript_ids,
      alternate_transcript_ids,
      cluster_A, cluster_B,
      delta_pi, alternate_delta_pi
    )
  
  if (nrow(pair_df) == 0) {
    stop("No rows found in DTU_results for the supplied transcript pair.")
  }
  
  # 2) Rows where dominant == tx_dom & alternate == tx_alt
  forward <- pair_df %>%
    filter(dominant_transcript_ids == tx_dom,
           alternate_transcript_ids == tx_alt)
  
  # A->B tiles (tx_dom -> tx_alt): delta_pi at (cluster_A, cluster_B)
  ab_forward <- forward %>%
    transmute(
      cluster_x = cluster_A,
      cluster_y = cluster_B,
      value     = delta_pi
    )
  
  # B->A tiles (tx_alt -> tx_dom): alternate_delta_pi at (cluster_B, cluster_A)
  ba_forward <- forward %>%
    transmute(
      cluster_x = cluster_B,
      cluster_y = cluster_A,
      value     = alternate_delta_pi
    )
  
  # 3) Rows where dominant == tx_alt & alternate == tx_dom
  reverse <- pair_df %>%
    filter(dominant_transcript_ids == tx_alt,
           alternate_transcript_ids == tx_dom)
  
  # In these rows, cluster_A is for tx_alt, cluster_B is for tx_dom.
  # We still want x = clusters for tx_dom, y = clusters for tx_alt.
  #
  # Each row contributes TWO mirror tiles, exactly as the forward branch does: the
  # tx_dom-perspective value at (cluster_A, cluster_B) and the tx_alt-perspective value
  # at (cluster_B, cluster_A). In a reverse row the roles of the two delta columns are
  # swapped, because `dominant_*` there refers to tx_alt.
  # Previously BOTH reverse tiles were emitted at (cluster_B, cluster_A), so one silently
  # overplotted the other and its mirror cell was left empty.
  
  # A->B (tx_dom -> tx_alt) here is alternate_delta_pi
  ab_reverse <- reverse %>%
    transmute(
      cluster_x = cluster_A,
      cluster_y = cluster_B,
      value     = alternate_delta_pi
    )
  
  # B->A (tx_alt -> tx_dom) here is delta_pi
  ba_reverse <- reverse %>%
    transmute(
      cluster_x = cluster_B,
      cluster_y = cluster_A,
      value     = delta_pi
    )
  
  # 4) Combine everything
  heat_df <- bind_rows(ab_forward, ba_forward, ab_reverse, ba_reverse)
  
  # 5) Build one unified cluster ordering, used for BOTH axes.
  # Numeric where the names allow it ("Cluster_10" after "Cluster_9", not after
  # "Cluster_1"), falling back to plain sorting otherwise. The fallback matters: with
  # as.integer() alone, any non-numeric cluster name becomes NA and sort() silently
  # DROPS it, which would delete that cluster's tiles from the plot without warning.
  all_clusters <- unique(c(as.character(heat_df$cluster_x),
                           as.character(heat_df$cluster_y)))
  
  cluster_nums <- suppressWarnings(as.integer(gsub("^Cluster_", "", all_clusters)))
  cluster_levels <- if (any(is.na(cluster_nums))) {
    sort(all_clusters)
  } else {
    all_clusters[order(cluster_nums)]
  }
  
  heat_df <- heat_df %>%
    mutate(
      cluster_x = factor(cluster_x, levels = cluster_levels),
      cluster_y = factor(cluster_y, levels = cluster_levels)
    )
  
  # 6) Diagonal coordinates
  diag_df <- data.frame(
    cluster_x = factor(cluster_levels, levels = cluster_levels),
    cluster_y = factor(cluster_levels, levels = cluster_levels)
  )
  
  # 7) Plot
  p <- ggplot(heat_df, aes(x = cluster_x, y = cluster_y, fill = value)) +
    geom_tile(color = "grey80") +
    geom_path(
      data = diag_df,
      aes(x = cluster_x, y = cluster_y, group = 1),
      inherit.aes = FALSE,
      color = "black",
      linewidth = 1.2,
      lineend = "round"
    ) +
    # Both axes are pinned to the SAME complete level set. Without drop = FALSE ggplot
    # drops levels unused on a given axis, so a cluster that only ever appears as a y
    # value vanishes from x and then gets re-appended at the end by the diagonal layer --
    # which is how the rows and columns ended up in different orders.
    scale_x_discrete(limits = cluster_levels, drop = FALSE) +
    scale_y_discrete(limits = cluster_levels, drop = FALSE) +
    coord_fixed() +
    scale_fill_gradient2(
      low = "purple",
      mid = "black",
      high = "yellow",
      midpoint = 0,
      name = expression(Delta*pi)
    ) +
    labs(
      x = paste0("Clusters for ", tx_dom),
      y = paste0("Clusters for ", tx_alt)
    ) +
    theme_bw() +
    theme(
      panel.grid = element_blank(),
      axis.text.x = element_text(angle = 45, hjust = 1)
    )
  
  return(p)
}


#####################################
# Well-supported alt-termini DTU examples
#####################################

# Two isoforms that share a splice pattern and differ only at one terminus are told apart
# by the EM mostly through reads compatible with both, so a significant DTU call can rest
# on how those reads were apportioned. These helpers keep the alt-termini calls whose
# isoforms are each measured rather than apportioned, then check the survivors against
# where the reads actually start or end.
#
# Model-level support, required of BOTH isoforms of a pair:
#   min_uniq_FSM_reads      uniquely assigned FSM reads, summed over the per-cluster quant runs
#   min_uniq_read_frac      share of the isoform's assigned reads that were assigned uniquely
#   min_terminus_reads      reads supporting the terminus that differs (TSS_read_count / PolyA_read_count)
#   PolyA sites             a PAS hexamer and no internal-priming flag
# and the differing termini at least min_termini_separation bp apart, so that one broad
# site called as two does not count as a switch.
#
# Expects gtf_parsed from parse_inputs, and dtu_results as read by parse_inputs (any
# gene_id / gene_symbol grouping; cross-gene rows should already be excluded).
get_alt_termini_support = function(dtu_results,
                                   cluster_quant_expr_tarball,
                                   transcript_id_mapping_tsv,
                                   min_uniq_FSM_reads = 5,
                                   min_uniq_read_frac = 0.1,
                                   min_terminus_reads = 20,
                                   min_termini_separation = 100) {

    quant_expr_dir = file.path(tempdir(), basename(cluster_quant_expr_tarball))
    untar(cluster_quant_expr_tarball, exdir = quant_expr_dir)

    # quant.expr ids lack the gene-symbol prefix; map them onto the ids used everywhere else
    transcript_id_to_symbol_id = read_tsv(transcript_id_mapping_tsv, col_types = cols(.default = "c")) %>%
        select(transcript_id, new_transcript_id) %>% distinct()

    isoform_read_support = list.files(quant_expr_dir, pattern = "quant.expr$", recursive = TRUE, full.names = TRUE) %>%
        map_dfr(~ read_tsv(.x, comment = "#", col_types = cols(.default = "c"))) %>%
        mutate(across(c(uniq_reads, all_reads, uniq_FSM_reads), as.numeric)) %>%
        group_by(transcript_id) %>%
        summarize(uniq_reads = sum(uniq_reads), all_reads = sum(all_reads),
                  uniq_FSM_reads = sum(uniq_FSM_reads), .groups = "drop") %>%
        inner_join(transcript_id_to_symbol_id, by = "transcript_id") %>%
        select(transcript_id = new_transcript_id, uniq_reads, all_reads, uniq_FSM_reads)

    isoform_termini = gtf_parsed %>% filter(feature == "transcript") %>%
        transmute(transcript_id,
                  TSS_pos = if_else(strand == "+", as.integer(start), as.integer(end)),
                  PolyA_pos = if_else(strand == "+", as.integer(end), as.integer(start)),
                  TSS_read_count = coalesce(as.numeric(TSS_read_count), 0),
                  PolyA_read_count = coalesce(as.numeric(PolyA_read_count), 0),
                  PolyA_called = coalesce(PolyA == "True", FALSE),
                  PAS = coalesce(PAS, "none"),
                  internal_priming = coalesce(InternalPriming == "True", FALSE)) %>%
        left_join(isoform_read_support, by = "transcript_id") %>%
        mutate(across(c(uniq_reads, all_reads, uniq_FSM_reads), ~ coalesce(.x, 0)),
               uniq_read_frac = uniq_reads / pmax(all_reads, 1),
               PolyA_ok = ! PolyA_called | (PAS != "none" & ! internal_priming))

    dtu_results %>%
        filter(as.character(significant) %in% c("True", "TRUE"),
               dominant_splice_hashcodes == alternate_splice_hashcodes,
               ! grepl(",", dominant_transcript_ids), ! grepl(",", alternate_transcript_ids)) %>%
        mutate(across(c(dominant_pi_A, dominant_pi_B, alternate_pi_A, alternate_pi_B),
                      ~ suppressWarnings(as.numeric(.x)))) %>%
        select(gene_symbol, cluster_A, cluster_B, pvalue, delta_pi, alternate_delta_pi,
               dominant_transcript_ids, alternate_transcript_ids,
               dominant_pi_A, dominant_pi_B, alternate_pi_A, alternate_pi_B) %>%
        inner_join(isoform_termini %>% rename_with(~ paste0("dom_", .x)), by = c("dominant_transcript_ids" = "dom_transcript_id")) %>%
        inner_join(isoform_termini %>% rename_with(~ paste0("alt_", .x)), by = c("alternate_transcript_ids" = "alt_transcript_id")) %>%
        mutate(TSS_separation = abs(dom_TSS_pos - alt_TSS_pos),
               PolyA_separation = abs(dom_PolyA_pos - alt_PolyA_pos),
               alt_terminus = case_when(TSS_separation >= min_termini_separation & PolyA_separation < min_termini_separation ~ "TSS",
                                        PolyA_separation >= min_termini_separation & TSS_separation < min_termini_separation ~ "PolyA",
                                        TSS_separation >= min_termini_separation ~ "both",
                                        TRUE ~ "neither"),
               dom_terminus_reads = if_else(alt_terminus == "TSS", dom_TSS_read_count, dom_PolyA_read_count),
               alt_terminus_reads = if_else(alt_terminus == "TSS", alt_TSS_read_count, alt_PolyA_read_count),
               # the model's within-pair share of the dominant isoform, comparable to the read-level share
               model_dom_share_A = dominant_pi_A / (dominant_pi_A + alternate_pi_A),
               model_dom_share_B = dominant_pi_B / (dominant_pi_B + alternate_pi_B),
               well_supported = alt_terminus %in% c("TSS", "PolyA") &
                   dom_uniq_FSM_reads >= min_uniq_FSM_reads & alt_uniq_FSM_reads >= min_uniq_FSM_reads &
                   dom_uniq_read_frac >= min_uniq_read_frac & alt_uniq_read_frac >= min_uniq_read_frac &
                   dom_terminus_reads >= min_terminus_reads & alt_terminus_reads >= min_terminus_reads &
                   dom_PolyA_ok & alt_PolyA_ok)
}


# Runs util/sc/diff_iso_usage/alt_termini_read_check.py over the candidate pairs and returns
# its table. The result is cached in output_tsv and recomputed only when the candidate set
# changes, since it reads the full BAM.
run_alt_termini_read_check = function(candidates, gtf, bam, cell_clusters, genome_fa, output_tsv,
                                      site_window = 50,
                                      lraa_root = Sys.getenv("LRAA_ROOT", "~/GITHUB/MDL/LongReadAlignmentAssembler")) {

    candidate_cols = c("gene_symbol", "alt_terminus", "dominant_transcript_ids", "alternate_transcript_ids",
                       "cluster_A", "cluster_B")

    candidates_tsv = paste0(output_tsv, ".candidates.tsv")
    new_candidates_tsv = tempfile(fileext = ".tsv")
    write_tsv(candidates %>% select(all_of(candidate_cols)) %>% distinct(), new_candidates_tsv)

    if (file.exists(output_tsv) && file.exists(candidates_tsv) &&
        unname(tools::md5sum(candidates_tsv)) == unname(tools::md5sum(new_candidates_tsv))) {
        message("-reusing read check in ", output_tsv)
    } else {
        file.copy(new_candidates_tsv, candidates_tsv, overwrite = TRUE)
        script = file.path(path.expand(lraa_root), "util/sc/diff_iso_usage/alt_termini_read_check.py")
        status = system2(script, c("--candidates", candidates_tsv, "--gtf", gtf, "--bam", bam,
                                   "--cell_clusters", cell_clusters, "--genome_fa", genome_fa,
                                   "--site_window", site_window, "--output", output_tsv))
        if (status != 0) {
            stop("alt_termini_read_check.py failed with status ", status)
        }
    }

    read_tsv(output_tsv, show_col_types = FALSE) %>%
        select(all_of(candidate_cols), dom_terminus_pos, alt_terminus_pos, n_reads,
               read_frac_at_dom, read_frac_at_alt, read_frac_elsewhere,
               reads_dom_A, reads_alt_A, read_dom_share_A, reads_dom_B, reads_alt_B, read_dom_share_B,
               top_read_end_peaks, dom_downstream_A_of_20, alt_downstream_A_of_20)
}


# A call is borne out by the reads when the share of reads ending at the dominant terminus
# moves between the pair's clusters in the model's direction, by at least
# min_read_share_delta and by at least min_read_vs_model_delta_ratio of the model's own
# shift (the EM must not be inflating it), with min_reads_per_cluster reads at either
# terminus in each cluster and no more than max_read_frac_elsewhere of the reads ending
# away from both termini (a smear there means the termini are not where the reads are).
confirm_alt_termini_by_reads = function(support_w_reads,
                                        min_read_share_delta = 0.15,
                                        min_read_vs_model_delta_ratio = 0.5,
                                        min_reads_per_cluster = 20,
                                        max_read_frac_elsewhere = 0.3) {
    support_w_reads %>%
        mutate(read_share_delta = read_dom_share_A - read_dom_share_B,
               model_share_delta = model_dom_share_A - model_dom_share_B,
               read_confirmed = coalesce(
                   sign(read_share_delta) == sign(model_share_delta) &
                   abs(read_share_delta) >= min_read_share_delta &
                   abs(read_share_delta) >= min_read_vs_model_delta_ratio * abs(model_share_delta) &
                   reads_dom_A + reads_alt_A >= min_reads_per_cluster &
                   reads_dom_B + reads_alt_B >= min_reads_per_cluster &
                   read_frac_elsewhere <= max_read_frac_elsewhere,
                   FALSE))
}


# Showcase examples among the read-confirmed comparisons: those whose reads shift by at
# least min_showcase_read_share_delta between the two clusters, one per gene (its largest
# shift), genes ordered by that shift. Ranked by the read shift rather than the p-value,
# since the p-value grows with depth while the shift is what a usage-fraction umap shows.
select_best_alt_termini_examples = function(confirmed, min_showcase_read_share_delta = 0.25) {
    confirmed %>%
        filter(read_confirmed) %>%
        group_by(gene_symbol) %>%
        mutate(n_confirmed_comparisons = n()) %>%
        filter(abs(read_share_delta) >= min_showcase_read_share_delta) %>%
        arrange(desc(abs(read_share_delta)), pvalue) %>% slice(1) %>% ungroup() %>%
        arrange(desc(abs(read_share_delta))) %>%
        transmute(alt_terminus, gene_symbol, dominant_transcript_ids, alternate_transcript_ids,
                  cluster_A, cluster_B,
                  read_share_delta = round(read_share_delta, 2),
                  read_dom_share_A, read_dom_share_B,
                  model_dom_share_A = round(model_dom_share_A, 2), model_dom_share_B = round(model_dom_share_B, 2),
                  delta_pi, alternate_delta_pi, pvalue, n_confirmed_comparisons,
                  TSS_separation, PolyA_separation,
                  dom_uniq_FSM_reads, alt_uniq_FSM_reads,
                  dom_uniq_read_frac = round(dom_uniq_read_frac, 2),
                  alt_uniq_read_frac = round(alt_uniq_read_frac, 2),
                  dom_terminus_reads, alt_terminus_reads, dom_PAS, alt_PAS, read_frac_elsewhere)
}


# Cluster colors for the comparison barplots: the first four categorical slots of the dataviz
# reference palette (colorblind-separable; aqua and yellow fall under 3:1 contrast on white,
# so every bar also carries its value as a label).
COMPARISON_CLUSTER_COLORS = c("#2a78d6", "#eb6834", "#1baf7a", "#eda100")


# Cell-type label per cluster, from whichever annotation column umap_df carries.
get_cluster_labels = function(max_chars = 24) {
    label_col = intersect(c("cell_type_simplified", "top_cell_type_annot", "cell_type_annot", "cas_cell_type_label_1"),
                          colnames(umap_df))[1]
    umap_df %>%
        count(seurat_clusters, cell_type = .data[[label_col]]) %>%
        group_by(seurat_clusters) %>%
        summarize(cell_type = cell_type[which.max(n)], n_cells = sum(n), .groups = "drop") %>%
        mutate(cluster = paste0("Cluster_", seurat_clusters),
               cluster_label = paste0(str_trunc(cell_type, max_chars), " (", seurat_clusters, "; ", n_cells, " cells)"))
}


# The pair's top significant DTU comparisons (largest |delta_pi|), one panel each: the two
# isoforms on the x axis, dodged bars for the two clusters, height = the isoform's fraction
# of the gene's reads in that cluster (pi, as tested).
plot_alt_termini_cluster_shifts = function(example, dtu_results, n_comparisons = 2, isoform_names = NULL) {

    dom_id = example$dominant_transcript_ids
    alt_id = example$alternate_transcript_ids
    iso_label = function(id) if (! is.null(isoform_names) && id %in% names(isoform_names)) isoform_names[[id]] else sub("^.*:", "", id)

    cluster_labels = get_cluster_labels()

    comparisons = dtu_results %>%
        filter(as.character(significant) %in% c("True", "TRUE"),
               (dominant_transcript_ids == dom_id & alternate_transcript_ids == alt_id) |
               (dominant_transcript_ids == alt_id & alternate_transcript_ids == dom_id)) %>%
        mutate(across(c(dominant_pi_A, dominant_pi_B, alternate_pi_A, alternate_pi_B), ~ suppressWarnings(as.numeric(.x)))) %>%
        arrange(desc(abs(delta_pi)), pvalue) %>%
        head(n_comparisons) %>%
        left_join(cluster_labels %>% transmute(cluster, short_A = paste0(str_trunc(cell_type, 22), " (", seurat_clusters, ")")),
                  by = c("cluster_A" = "cluster")) %>%
        left_join(cluster_labels %>% transmute(cluster, short_B = paste0(str_trunc(cell_type, 22), " (", seurat_clusters, ")")),
                  by = c("cluster_B" = "cluster")) %>%
        mutate(comparison = paste0(short_A, " vs ", short_B, "\nadj. p = ", signif(adj_pvalue, 2)),
               comparison = factor(comparison, levels = comparison))

    # one row per comparison x cluster x isoform, oriented to the showcased pair's isoforms
    bars = comparisons %>%
        mutate(same_orientation = dominant_transcript_ids == dom_id) %>%
        transmute(comparison, cluster_A, cluster_B,
                  dom_pi_A = if_else(same_orientation, dominant_pi_A, alternate_pi_A),
                  dom_pi_B = if_else(same_orientation, dominant_pi_B, alternate_pi_B),
                  alt_pi_A = if_else(same_orientation, alternate_pi_A, dominant_pi_A),
                  alt_pi_B = if_else(same_orientation, alternate_pi_B, dominant_pi_B)) %>%
        pivot_longer(c(dom_pi_A, dom_pi_B, alt_pi_A, alt_pi_B), names_to = "which", values_to = "pi") %>%
        mutate(cluster = if_else(str_ends(which, "_A"), cluster_A, cluster_B),
               isoform = factor(if_else(str_starts(which, "dom"), iso_label(dom_id), iso_label(alt_id)),
                                levels = c(iso_label(dom_id), iso_label(alt_id)))) %>%
        inner_join(cluster_labels %>% select(cluster, cluster_label), by = "cluster")

    # color follows the cluster, in order of first appearance, across both panels
    cluster_levels = unique(bars$cluster_label)
    bars = bars %>% mutate(cluster_label = factor(cluster_label, levels = cluster_levels))

    dodge = position_dodge(width = 0.8)

    ggplot(bars, aes(x = isoform, y = pi, fill = cluster_label)) +
        geom_col(position = dodge, width = 0.75, color = "white", linewidth = 0.5) +
        geom_text(aes(label = sprintf("%.2f", pi)), position = dodge, vjust = -0.4, size = 3, color = "#0b0b0b") +
        facet_wrap(~ comparison, nrow = 1) +
        scale_fill_manual(values = setNames(COMPARISON_CLUSTER_COLORS[seq_along(cluster_levels)], cluster_levels),
                          name = NULL) +
        scale_y_continuous(limits = c(0, 1.08), breaks = seq(0, 1, 0.25), expand = expansion(mult = c(0, 0))) +
        labs(x = NULL, y = "isoform fraction of gene reads (pi)",
             title = paste0(example$gene_symbol, ": top significant cluster comparisons")) +
        guides(fill = guide_legend(ncol = 2)) +
        theme_minimal(base_size = 10) +
        theme(legend.position = "bottom", panel.grid.major.x = element_blank(), panel.grid.minor = element_blank(),
              strip.text = element_text(face = "bold"))
}


# Names an alt-termini pair's isoforms by the terminus that tells them apart, e.g.
# "iso-8: distal PolyA" / "iso-4: proximal PolyA", or "upstream TSS" / "downstream TSS".
get_alt_termini_isoform_names = function(example) {
    ids = c(example$dominant_transcript_ids, example$alternate_transcript_ids)
    tx = gtf_parsed %>% filter(feature == "transcript", transcript_id %in% ids) %>%
        distinct(transcript_id, .keep_all = TRUE) %>%
        mutate(start = as.integer(start), end = as.integer(end),
               TSS_pos = if_else(strand == "+", start, end), PolyA_pos = if_else(strand == "+", end, start))
    tx = tx[match(ids, tx$transcript_id), ]
    plus = tx$strand[1] == "+"
    # position along the transcript's direction: larger = further downstream
    downstream_rank = function(pos) rank(if (plus) pos else -pos)
    if (example$alt_terminus == "TSS") {
        where = if_else(downstream_rank(tx$TSS_pos) == 2, "downstream TSS", "upstream TSS")
    } else {
        where = if_else(downstream_rank(tx$PolyA_pos) == 2, "distal PolyA", "proximal PolyA")
    }
    setNames(paste0(sub("^.*:", "", ids), ": ", where), ids)
}


# A compact figure for one alt-termini example: the SNN-smoothed usage-fraction umaps of
# the pair (each isoform's share of the pair's reads per cell) over the barplots of its
# top significant cluster comparisons.
plot_alt_termini_umaps_and_shifts = function(example, dtu_results, smooth_graph, file = NULL,
                                             width = 10, height = 9) {

    isoform_names = get_alt_termini_isoform_names(example)

    p_umap = plot_isoform_usage_fraction_umap(example$gene_symbol, names(isoform_names),
                                              denominator = "selected",
                                              smooth_graph = smooth_graph,
                                              cell_evidence = "gene",
                                              transcript_labels = isoform_names,
                                              ncol = 2)

    n_cells = n_distinct(p_umap$layers[[2]]$data$cell_barcode)
    terminus = if (example$alt_terminus == "TSS") "TSS" else "PolyA"
    panel_margin = margin(t = 18, r = 6, b = 6, l = 18)   # room for the A / B panel letters

    p_umap = p_umap +
        labs(title = paste0(example$gene_symbol, ": alternative ", terminus, " usage per cell"),
             subtitle = paste0("share of the pair's reads per cell, SNN-smoothed (", format(n_cells, big.mark = ","), " cells)")) +
        theme(plot.margin = panel_margin)

    p_bars = plot_alt_termini_cluster_shifts(example, dtu_results, isoform_names = isoform_names) +
        labs(title = "Top significant cluster comparisons") +
        theme(plot.margin = panel_margin)

    p = plot_grid(p_umap, p_bars, ncol = 1, rel_heights = c(1.15, 1), labels = c("A", "B"))

    if (! is.null(file)) {
        ggsave(p, file = file, width = width, height = height)
    }

    p
}


plot_alt_termini_example = function(example, file_prefix, smooth_graph = NULL, dtu_results = NULL) {

    p = make_diff_iso_usage_compound_plot(example$gene_symbol, 0,
                                          ignore_unspliced = FALSE,
                                          transcript_ids = c(example$dominant_transcript_ids, example$alternate_transcript_ids),
                                          usage_fraction_umap = TRUE, denominator = "selected",
                                          smooth_graph = smooth_graph)
    height = 8

    if (! is.null(dtu_results)) {
        p = plot_grid(p, plot_alt_termini_cluster_shifts(example, dtu_results), ncol = 1, rel_heights = c(2, 1.1))
        height = 12
    }

    ggsave(p, file = paste0(example$gene_symbol, ".", file_prefix, ".pdf"), width = 11, height = height)

    p
}
