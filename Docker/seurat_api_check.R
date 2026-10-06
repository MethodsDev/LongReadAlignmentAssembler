# Build-time check that the Seurat stack in lraa-sc-base works as one set.
#
# Seurat, SeuratObject and ggplot2 are installed separately, and a mismatch does not
# fail at library() -- it fails inside the calls below, at analysis time. Seurat 5.0.0
# (from the GitHub seurat5 branch) against current CRAN SeuratObject 5.4 and ggplot2
# 4.0 broke FindMarkers, AddModuleScore, FeaturePlot and VlnPlot (slot= is defunct)
# and DotPlot with a named feature list (facet_grid(facets=) is defunct). Run by
# Dockerfile.sc-base; a non-zero exit fails the image build.
#
# Adapted from a check written against lraa-sc:0.44.2. GetAssayData(slot=) is tested
# as layer= here: slot= is defunct in SeuratObject 5.x and cannot pass. AddModuleScore
# uses ctrl = 5 because its default draws 100 control genes per expression bin, more
# than this 200-gene matrix holds.

suppressPackageStartupMessages({library(Seurat); library(ggplot2)})
for (p in c("Seurat", "SeuratObject", "ggplot2", "Matrix")) cat(sprintf("%-13s %s\n", p, packageVersion(p)))
set.seed(1)
m <- matrix(rpois(200 * 300, 1), 200, 300, dimnames = list(paste0("G", 1:200), paste0("C", 1:300)))
so <- CreateSeuratObject(as(m, "dgCMatrix"))
so <- NormalizeData(so, verbose = FALSE); so <- FindVariableFeatures(so, verbose = FALSE); so <- ScaleData(so, verbose = FALSE)
so <- RunPCA(so, npcs = 10, verbose = FALSE); so <- RunUMAP(so, dims = 1:10, verbose = FALSE)
so$grp <- sample(c("a", "b"), 300, TRUE); pdf(NULL)
check <- function(name, expr) {
  r <- tryCatch({force(expr); "OK"}, error = function(e) paste("FAIL:", conditionMessage(e)))
  cat(sprintf("%-24s %s\n", name, substr(gsub("\n", " ", r), 1, 160))); startsWith(r, "OK")
}
ok <- c(
  check("presto installed", stopifnot(requireNamespace("presto", quietly = TRUE))),
  check("AddModuleScore", AddModuleScore(so, features = list(paste0("G", 1:10)), name = "s", ctrl = 5)),
  check("FindMarkers", FindMarkers(so, ident.1 = "a", group.by = "grp", logfc.threshold = 0, min.pct = 0)),
  check("GetAssayData(layer=)", GetAssayData(so, layer = "data")),
  check("LayerData", LayerData(so, layer = "counts")),
  check("DotPlot(named list)", print(DotPlot(so, features = list(A = paste0("G", 1:3), B = paste0("G", 4:6)), group.by = "grp"))),
  check("DimPlot", print(DimPlot(so, group.by = "grp"))),
  check("FeaturePlot", print(FeaturePlot(so, features = "G1"))),
  check("VlnPlot", print(VlnPlot(so, features = "G1", group.by = "grp"))))
quit(status = if (all(ok)) 0 else 1)
