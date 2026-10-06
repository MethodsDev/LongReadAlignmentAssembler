#!/usr/bin/env Rscript

# Differential TSS or PolyA site usage across cell clusters with DEXSeq on
# pseudo-replicates.
#
# The statistical approach is not new here: testing differential usage of
# polyA sites (peaks) between single-cell populations with DEXSeq on pseudo-bulk
# replicates formed by aggregating each population's cells is the approach of
# Sierra (Patrick et al. 2020, Genome Biology 21:167,
# doi:10.1186/s13059-020-02071-7; its DUTest) and SCAPE (Zhou et al. 2022,
# Nucleic Acids Research 50:e66, doi:10.1093/nar/gkac167; cells shuffled into
# six pseudo-replicates for DEXSeq). This script applies it to LRAA's TSS and
# PolyA sites, with per-cell counts from long-read read ends
# (util/sc/site_read_support_to_sparse_matrix.py), and adds stageR site
# confirmation and seed-stability. DEXSeq: Anders, Reyes & Huber 2012, Genome
# Research 22:2008, doi:10.1101/gr.133744.111. stageR: Van den Berge et al.
# 2017, Genome Biology 18:151, doi:10.1186/s13059-017-1277-0.
#
# Each cluster's cells are dealt at random into K pseudo-replicates and each
# replicate's read ends are summed per site (sites play the part of DEXSeq's
# exon bins, genes its groups). Per gene, DEXSeq compares
#     ~ sample + exon + cluster:exon   against   ~ sample + exon
# so a gene is called when the share of its read ends at some site depends on
# the cluster. perGeneQValue gives the gene-level FDR; stageR (method "dtu")
# then confirms which sites in the significant genes carry the change, at the
# same overall FDR.
#
# Pseudo-replicates from one sample measure cell-to-cell sampling noise only,
# not sample-to-sample variation, so p-values are optimistic in absolute terms
# (Squair et al. 2021). Dealing cells is random, so the whole test is repeated
# under --n_seeds seeds and a gene is reported as stable when it is significant
# under at least --min_stable_seeds of them.
#
# Pairwise cluster contrasts are then run, under the first seed, on the stable
# genes only: for each pair of clusters a DEXSeq LRT on that pair's
# pseudo-replicates, reusing the dispersions fitted on all clusters (a pair has
# too few replicates to fit its own).
#
# Outputs (--output_prefix.<kind>.*):
#   seed<N>.sites.tsv.gz  per site: DEXSeq p, padj; gene q; stageR site padj
#   genes.tsv             per gene: q under each seed, seeds significant, stable
#   sites.tsv             per site: stageR padj under each seed, seeds confirmed
#   cluster_usage.tsv.gz  per site x cluster: read ends, gene read ends, usage
#   pairwise.tsv.gz       per stable gene x cluster pair x site: DEXSeq p, padj
#                         (BH over all pairwise tests), usage in each cluster

suppressPackageStartupMessages({
    library(argparse)
    library(Matrix)
    library(data.table)
    library(DEXSeq)
    library(stageR)
    library(BiocParallel)
})

parser <- ArgumentParser()
parser$add_argument("--counts_prefix", required=TRUE, help="output prefix given to count_site_read_ends.py")
parser$add_argument("--sites", required=TRUE, help="site table from prep_site_table.py")
parser$add_argument("--kind", required=TRUE, choices=c("TSS", "PolyA"))
parser$add_argument("--pseudoreps", type="integer", default=3, help="pseudo-replicates per cluster")
parser$add_argument("--n_seeds", type="integer", default=5)
parser$add_argument("--min_stable_seeds", type="integer", default=4)
parser$add_argument("--min_cluster_cells", type="integer", default=30,
                    help="clusters with fewer cells are left out")
parser$add_argument("--min_site_reads", type="integer", default=10, help="site read ends over all cells")
parser$add_argument("--min_site_usage", type="double", default=0.1,
                    help="a site is tested only if it takes at least this share of its gene's read ends in at least one cluster ...")
parser$add_argument("--min_gene_cluster_reads", type="integer", default=20,
                    help="... with at least this many gene read ends in that cluster")
parser$add_argument("--fdr", type="double", default=0.05)
parser$add_argument("--cores", type="integer", default=8)
parser$add_argument("--no_pairwise", action="store_true")
parser$add_argument("--max_genes", type="integer", default=0, help="test only this many genes, chosen at random (for trial runs)")
parser$add_argument("--output_prefix", required=TRUE)
args <- parser$parse_args()

msg <- function(...) message(format(Sys.time(), "%H:%M:%S "), sprintf(...))
bp <- if (args$cores > 1) MulticoreParam(args$cores) else SerialParam()
out <- function(suffix) paste0(args$output_prefix, ".", args$kind, ".", suffix)

## ---- counts: sites x cells

counts <- readMM(paste0(args$counts_prefix, ".site_counts.mtx.gz"))
counts <- as(counts, "CsparseMatrix")
site_ids <- fread(paste0(args$counts_prefix, ".sites.tsv.gz"), header=FALSE)$V1
cells <- fread(paste0(args$counts_prefix, ".barcodes.tsv.gz"), header=FALSE, col.names=c("barcode", "cluster"))
rownames(counts) <- site_ids
colnames(counts) <- cells$barcode

sites <- fread(args$sites)
stopifnot(identical(sites$site_id, site_ids))
sites <- sites[kind == args$kind & competing == TRUE & gene_key != ""]

cluster_sizes <- cells[, .N, by=cluster]
keep_clusters <- cluster_sizes[N >= args$min_cluster_cells, cluster]
cells <- cells[cluster %in% keep_clusters]
clusters <- cells[, unique(cluster)]
clusters <- clusters[order(as.integer(sub("Cluster_", "", clusters)))]
msg("%s: %d competing sites, %d cells in %d clusters", args$kind, nrow(sites), nrow(cells), length(clusters))

counts <- counts[sites$site_id, cells$barcode]

## ---- per-cluster usage, and the site filter (blind to the pseudo-replicates, so the same under every seed)

pseudobulk <- function(m, groups) {
    g <- factor(groups, levels=unique(groups))
    ind <- sparseMatrix(i=seq_along(g), j=as.integer(g), x=1, dims=c(length(g), nlevels(g)),
                        dimnames=list(NULL, levels(g)))
    as.matrix(m %*% ind)
}

cluster_counts <- pseudobulk(counts, cells$cluster)[, clusters, drop=FALSE]

usage_table <- function(cm, site_tab) {
    gene_tot <- rowsum(cm, site_tab$gene_key)[site_tab$gene_key, , drop=FALSE]
    dt <- data.table(site_id=rep(rownames(cm), ncol(cm)),
                     gene_key=rep(site_tab$gene_key, ncol(cm)),
                     cluster=rep(colnames(cm), each=nrow(cm)),
                     reads=as.vector(cm), gene_reads=as.vector(gene_tot))
    dt[, usage := ifelse(gene_reads > 0, reads / gene_reads, NA_real_)]
    dt
}

usage_all <- usage_table(cluster_counts, sites)
site_ok <- usage_all[, .(passes=any(gene_reads >= args$min_gene_cluster_reads & usage >= args$min_site_usage, na.rm=TRUE)),
                     by=site_id]
site_tot <- rowSums(cluster_counts)
keep <- sites$site_id %in% site_ok[passes == TRUE, site_id] & site_tot[sites$site_id] >= args$min_site_reads
sites <- sites[keep]
sites <- sites[, if (.N >= 2) .SD, by=gene_key]
if (args$max_genes > 0) {
    set.seed(1)
    trial_genes <- sample(unique(sites$gene_key), min(args$max_genes, uniqueN(sites$gene_key)))
    sites <- sites[gene_key %in% trial_genes]
}
# DEXSeq strips ':' from its ids, so it is given plain tags and the results are mapped back
sites[, feature_tag := paste0("s", .I)]
sites[, group_tag := paste0("g", match(gene_key, unique(gene_key)))]
site_of_tag <- setNames(sites$site_id, sites$feature_tag)
gene_of_tag <- setNames(sites$gene_key, sites$group_tag)
counts <- counts[sites$site_id, ]
cluster_counts <- cluster_counts[sites$site_id, , drop=FALSE]
msg("after filtering: %d sites in %d genes with >= 2 sites", nrow(sites), uniqueN(sites$gene_key))

# usage over the tested sites only: the denominator is the gene's read ends at tested sites
usage <- usage_table(cluster_counts, sites)
fwrite(usage, out("cluster_usage.tsv.gz"), sep="\t")

## ---- DEXSeq per seed

sample_table_for_seed <- function(seed) {
    set.seed(seed)
    rep_of <- integer(nrow(cells))
    for (cl in clusters) {
        idx <- which(cells$cluster == cl)
        rep_of[idx] <- sample(rep_len(seq_len(args$pseudoreps), length(idx)))
    }
    paste0(cells$cluster, ".r", rep_of)
}

run_dexseq <- function(cm, sample_clusters) {
    sample_data <- data.frame(row.names=colnames(cm), cluster=factor(sample_clusters, levels=unique(sample_clusters)))
    dxd <- DEXSeqDataSet(cm, sample_data, design=~sample + exon + cluster:exon,
                         featureID=sites$feature_tag[match(rownames(cm), sites$site_id)],
                         groupID=sites$group_tag[match(rownames(cm), sites$site_id)])
    # The default median-of-ratios size factors need at least one site with no zero
    # count in any pseudo-replicate. Genome-wide there always is; with only a few genes
    # (a small test set, a --max_genes trial) there may be none, so fall back to DESeq2's
    # "poscounts" geometric means (zeros left out), which DESeq2 recommends for sparse
    # counts. Untouched whenever the default works.
    tryCatch(estimateSizeFactors(dxd), error = function(e) {
        cts <- featureCounts(dxd)
        geo_means <- apply(cts, 1, function(x) if (all(x == 0)) 0 else exp(sum(log(x[x > 0])) / length(x)))
        msg("size factors: every site has a zero in some pseudo-replicate; using poscounts geometric means")
        sizeFactors(dxd) <- DESeq2::estimateSizeFactorsForMatrix(cts, geoMeans = geo_means)
        dxd
    })
}

seed_results <- list()
dxd_first <- NULL
pb_first <- NULL
for (seed in seq_len(args$n_seeds)) {
    msg("seed %d: DEXSeq", seed)
    samples <- sample_table_for_seed(seed)
    pb <- pseudobulk(counts, samples)
    storage.mode(pb) <- "integer"
    pb <- pb[, order(match(sub("\\.r[0-9]+$", "", colnames(pb)), clusters), colnames(pb))]
    dxd <- run_dexseq(pb, sub("\\.r[0-9]+$", "", colnames(pb)))
    dxd <- estimateDispersions(dxd, BPPARAM=bp)
    dxd <- testForDEU(dxd, BPPARAM=bp)
    res <- DEXSeqResults(dxd, independentFiltering=FALSE)
    gene_q <- perGeneQValue(res)

    site_res <- data.table(gene_key=gene_of_tag[res$groupID], site_id=site_of_tag[res$featureID],
                           dispersion=res$dispersion, pvalue=res$pvalue, padj=res$padj)
    site_res[, gene_q := gene_q[res$groupID]]
    names(gene_q) <- gene_of_tag[names(gene_q)]

    # stageR: screen on the gene q, confirm sites within the significant genes
    # (stageR splits ids at ':', so it too gets the plain tags)
    p_conf <- matrix(site_res$pvalue, ncol=1, dimnames=list(res$featureID, "site"))
    p_conf[is.na(p_conf)] <- 1
    p_screen <- gene_q
    stage <- stageRTx(pScreen=p_screen, pConfirmation=p_conf, pScreenAdjusted=TRUE,
                      tx2gene=data.frame(res$featureID, site_res$gene_key))
    stage <- stageWiseAdjustment(stage, method="dtu", alpha=args$fdr, allowNA=TRUE)
    adj <- as.data.table(getAdjustedPValues(stage, order=FALSE, onlySignificantGenes=FALSE))
    site_res[, stageR_site_padj := adj$transcript[match(res$featureID, adj$txID)]]

    fwrite(site_res, out(sprintf("seed%d.sites.tsv.gz", seed)), sep="\t")
    seed_results[[seed]] <- site_res
    msg("seed %d: %d of %d genes at q < %g", seed, sum(gene_q < args$fdr, na.rm=TRUE), length(gene_q), args$fdr)

    if (seed == 1) {
        dxd_first <- dxd
        pb_first <- pb
    } else {
        rm(dxd)
    }
    gc()
}

## ---- stability over seeds

genes <- seed_results[[1]][, .(n_sites=.N), by=gene_key]
for (seed in seq_along(seed_results)) {
    g <- unique(seed_results[[seed]][, .(gene_key, gene_q)])
    genes[, (paste0("q_seed", seed)) := g$gene_q[match(gene_key, g$gene_key)]]
}
qcols <- grep("^q_seed", names(genes), value=TRUE)
genes[, n_seeds_significant := rowSums(as.matrix(.SD) < args$fdr, na.rm=TRUE), .SDcols=qcols]
genes[, median_q := apply(as.matrix(.SD), 1, median, na.rm=TRUE), .SDcols=qcols]
genes[, stable := n_seeds_significant >= args$min_stable_seeds]
genes[, c("gene_symbol", "chrom", "strand") := tstrsplit(gene_key, "|", fixed=TRUE)]
fwrite(genes, out("genes.tsv"), sep="\t")

site_tab <- seed_results[[1]][, .(gene_key, site_id)]
for (seed in seq_along(seed_results)) {
    s <- seed_results[[seed]]
    site_tab[, (paste0("stageR_padj_seed", seed)) := s$stageR_site_padj[match(site_id, s$site_id)]]
}
pcols <- grep("^stageR_padj_seed", names(site_tab), value=TRUE)
site_tab[, n_seeds_confirmed := rowSums(as.matrix(.SD) < args$fdr, na.rm=TRUE), .SDcols=pcols]
fwrite(site_tab, out("sites.tsv"), sep="\t")
msg("stable genes (q < %g in >= %d of %d seeds): %d of %d", args$fdr, args$min_stable_seeds, args$n_seeds,
    sum(genes$stable), nrow(genes))

## ---- pairwise contrasts on the stable genes, first seed

if (!args$no_pairwise && any(genes$stable)) {
    stable_sites <- sites[gene_key %in% genes[stable == TRUE, gene_key], site_id]
    disp <- setNames(mcols(dxd_first)$dispersion, as.character(mcols(dxd_first)$featureID))
    sample_cluster <- sub("\\.r[0-9]+$", "", colnames(pb_first))
    pairs <- combn(clusters, 2, simplify=FALSE)
    msg("pairwise: %d stable genes (%d sites) x %d cluster pairs", genes[stable == TRUE, .N], length(stable_sites), length(pairs))

    pair_res <- bplapply(pairs, function(pr) {
        cols <- sample_cluster %in% pr
        cm <- pb_first[stable_sites, cols, drop=FALSE]
        cm <- cm[, order(match(sample_cluster[cols], pr))]
        pd <- run_dexseq(cm, sub("\\.r[0-9]+$", "", colnames(cm)))
        # a pair's 2 K replicates are too few to fit dispersions; take the all-cluster fit
        dispersions(pd) <- disp[as.character(mcols(pd)$featureID)]
        pd <- testForDEU(pd)
        r <- DEXSeqResults(pd, independentFiltering=FALSE)
        data.table(gene_key=gene_of_tag[r$groupID], site_id=site_of_tag[r$featureID], cluster_A=pr[1], cluster_B=pr[2],
                   pvalue=r$pvalue)
    }, BPPARAM=bp)
    pair_res <- rbindlist(pair_res)
    pair_res[, padj := p.adjust(pvalue, method="BH")]

    ua <- usage[, .(site_id, cluster_A=cluster, reads_A=reads, gene_reads_A=gene_reads, usage_A=usage)]
    ub <- usage[, .(site_id, cluster_B=cluster, reads_B=reads, gene_reads_B=gene_reads, usage_B=usage)]
    pair_res <- merge(pair_res, ua, by=c("site_id", "cluster_A"))
    pair_res <- merge(pair_res, ub, by=c("site_id", "cluster_B"))
    pair_res[, delta_usage := usage_B - usage_A]
    setcolorder(pair_res, c("gene_key", "site_id", "cluster_A", "cluster_B"))
    fwrite(pair_res, out("pairwise.tsv.gz"), sep="\t")
    msg("pairwise: %d site tests, %d at padj < %g", nrow(pair_res), sum(pair_res$padj < args$fdr, na.rm=TRUE), args$fdr)
}

msg("done")
