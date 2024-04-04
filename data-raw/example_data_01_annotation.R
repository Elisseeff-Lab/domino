# example_data_01_annotation

library(DropletUtils)
library(EnsDb.Hsapiens.v86)
library(scater)
library(scran)

temp_dir <- tempdir()

# Preparing SCE object for use in vignettes
# Data
pbmc_url <- "https://cf.10xgenomics.com/samples/cell/pbmc3k/pbmc3k_filtered_gene_bc_matrices.tar.gz"
pbmc_tar <- paste0(temp_dir, "/pbmc3k_filtered_gene_bc_matrices.tar.gz")
download.file(url = pbmc_url, destfile = pbmc_tar)
pbmc_dir <- paste0(temp_dir, "/pbmc3k_filtered_gene_bc_matrices")
untar(tarfile = pbmc_tar, exdir = pbmc_dir)
fname <- file.path(pbmc_dir, "/filtered_gene_bc_matrices/hg19")

sce <- read10xCounts(fname, col.names = TRUE)
rownames(sce) <- uniquifyFeatureNames(rowData(sce)$ID, rowData(sce)$Symbol)
location <- mapIds(EnsDb.Hsapiens.v86, keys=rowData(sce)$ID, column="SEQNAME", keytype="GENEID")

# QC
stats <- perCellQCMetrics(sce, subsets=list(Mito=which(location=="MT")))
high.mito <- isOutlier(stats$subsets_Mito_percent, type="higher")

colData(sce) <- cbind(colData(sce), stats)
sce$discard <- high.mito
sce <- sce[,!high.mito]

# Remove features with no expression
feat_keep <- apply(
  assay(sce, "counts"), 
  MARGIN = 1,
  FUN = function(x) {sum(x) > 0}
)
sce <- sce[feat_keep,]

# Normalization
set.seed(1000)
clusters <- quickCluster(sce)
sce <- computeSumFactors(sce, cluster=clusters)
sce <- logNormCounts(sce)

# Variance modeling
set.seed(1001)
dec <- modelGeneVarByPoisson(sce)
top <- getTopHVGs(dec, prop=0.1)

# Dimensionality reduction
set.seed(10000)
sce <- denoisePCA(sce, subset.row=top, technical=dec)

set.seed(100000)
sce <- runTSNE(sce, dimred="PCA")

set.seed(1000000)
sce <- runUMAP(sce, dimred="PCA")

# Clustering
g <- buildSNNGraph(sce, k=10, use.dimred = 'PCA')
clust <- igraph::cluster_walktrap(g)$membership
colLabels(sce) <- factor(clust)

#Cell type annotation
annot <- c(
  "3" = "naive_CD4_T_cell",
  "2" = "memory_CD4_T_cell",
  "1" = "CD8_T_cell",
  "7" = "CD8_T_cell",
  "5" = "CD14_monocyte",
  "6" = "CD14_monocyte",
  "8" = "B_cell",
  "4" = "CD16_monocyte",
  "10" = "NK_cell",
  "9" = "dendritic_cell",
  "12" = "Platelet",
  "11" = "CD8_T_cell"
)

sce$cell_type <- annot[as.character(sce$label)]

saveRDS(sce, paste0(temp_dir, "/pbmc3k_sce.rds"))
