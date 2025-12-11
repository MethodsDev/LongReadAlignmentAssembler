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
                        diff_iso_usage_stats_filename) {
    
    
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
    
    
    message("-parsing diff iso usage stats: ", diff_iso_usage_stats_filename)
    diff_iso_usage_stats <<- read.csv(diff_iso_usage_stats_filename, sep="\t", header=T)
    
    
    message("-parsing gtf file: ", gtf_filename)
    gtf_parsed <<- parse_gtf_file(gtf_filename)
    
}

parse_gtf_file = function(gtf_filename) {

    gtf_parsed_file = paste0(gtf_filename, ".rds")
    
    if (file.exists(gtf_parsed_file)) {
        
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
        transcript_ids = NULL) {
    
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
    
    if (min_isoform_frac_expr_any_cluster > 0) {
        
        transcript_ids = names(which(rowSums(isoform_frac_expr >= min_isoform_frac_expr_any_cluster) > 0))
        
        if (length(transcript_ids) < 2) {
            stop("Too few isoforms left after filtering based on min isoform fraction") 
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

make_diff_iso_usage_compound_plot = function(gene_of_interest, min_iso_fraction = 0, ignore_unspliced=FALSE,
                                             transcript_ids=NULL) {
    
    p_exon_expr_info = get_expression_ggplot2_heatmap_w_exon_structures(gene_of_interest, min_iso_fraction, 
                                                                        ignore_unspliced, transcript_ids)
    
    p_exon_expr = p_exon_expr_info$plot
    
    p_umap = plot_isoform_umap(gene_of_interest, p_exon_expr_info$transcript_ids)
    
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
  
  # A->B (tx_dom -> tx_alt) here is alternate_delta_pi at (cluster_B, cluster_A)
  ab_reverse <- reverse %>%
    transmute(
      cluster_x = cluster_B,
      cluster_y = cluster_A,
      value     = alternate_delta_pi
    )
  
  # B->A (tx_alt -> tx_dom) here is delta_pi at (cluster_A, cluster_B)
  ba_reverse <- reverse %>%
    transmute(
      cluster_x = cluster_A,
      cluster_y = cluster_B,
      value     = delta_pi
    )
  
  # 4) Combine everything
  heat_df <- bind_rows(ab_forward, ba_forward, ab_reverse, ba_reverse)
  
  # 5) Build unified cluster ordering (numeric order)
  all_clusters <- unique(c(as.character(heat_df$cluster_x),
                           as.character(heat_df$cluster_y)))
  
  cluster_levels <- all_clusters %>%
    gsub("^Cluster_", "", .) %>%   # strip prefix
    as.integer() %>%
    sort() %>%
    paste0("Cluster_", .)          # rebuild ordered names
  
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
    coord_fixed() +
    scale_fill_viridis_c(name = expression(Delta*pi)) +
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
