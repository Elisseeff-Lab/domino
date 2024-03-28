# Preparing example data from 10x PBMC3k dataset for use in vignettes

library(DropletUtils)
library(EnsDb.Hsapiens.v86)
library(scater)
library(scran)
library(dominoSignal)

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

saveRDS(sce, paste0(temp_dir, "pbmc_sce.rds"))

# Load SCENIC outputs
data_url <- "https://zenodo.org/records/10884027/files"

scenic_dir <- paste0(temp_dir, "/scenic")
if (!dir.exists(scenic_dir)) {
  dir.create(scenic_dir)
}
download.file(url = paste0(data_url, "/scenic_auc_pbmc_3k.csv"),
              destfile = paste0(scenic_dir, "/auc_pbmc_3k.csv"))
download.file(url = paste0(data_url, "/scenic_regulons_pbmc_3k.csv"),
              destfile = paste0(scenic_dir, "/regulons_pbmc_3k.csv"))
auc <- read.table(paste0(scenic_dir, "/auc_pbmc_3k.csv"),
                header = TRUE, row.names = 1,
                stringsAsFactors = FALSE, sep = ",")
regulons <- read.csv(paste0(scenic_dir, "/regulons_pbmc_3k.csv"))

# Create list of SCENIC regulons
regulons <- regulons[-1:-2,]
colnames(regulons) <- c("TF", "MotifID", "AUC", "NES", "MotifSimilarityQvalue", "OrthologousIdentity", "Annotation", "Context", "TargetGenes", "RankAtMax")
regulon_list <- create_regulon_list_scenic(regulons = regulons)

# Load AUC and adjust names
auc_in <- as.data.frame(t(auc))
# Remove pattern "..." from the end of all rownames:
rownames(auc_in) <- gsub("\\.\\.\\.$", "", rownames(auc_in))

# Load CellPhoneDB Database
cellphone_url <- "https://github.com/ventolab/cellphonedb-data/archive/refs/tags/v4.0.0.tar.gz"
cellphone_tar <- paste0(temp_dir, "/cellphoneDB_v4.tar.gz")
download.file(url = cellphone_url, destfile = cellphone_tar)
cellphone_dir <- paste0(temp_dir, "/cellphoneDB_v4")
untar(tarfile = cellphone_tar, exdir = cellphone_dir)
cellphone_data <- paste0(cellphone_dir, "/cellphonedb-data-4.0.0/data")

interactions <- read.csv(paste0(cellphone_data, "/interaction_input.csv"), stringsAsFactors = FALSE)
complexes <- read.csv(paste0(cellphone_data, "/complex_input.csv"), stringsAsFactors = FALSE)
genes <- read.csv(paste0(cellphone_data, "/gene_input.csv"), stringsAsFactors = FALSE)
proteins <- read.csv(paste0(cellphone_data, "/protein_input.csv"), stringsAsFactors = FALSE)

# Make RL Map
rl_map <- create_rl_map_cellphonedb(
    genes = genes, 
    proteins = proteins, 
    interactions = interactions,
    complexes = complexes,
    database_name = "CellPhoneDB_v4.0" # database version used
)

# Prepare inputs for pbmc object:
counts = assay(sce, "counts")
z_scores = t(scale(t(assay(sce, "logcounts"))))
clusters = factor(sce$cell_type)

# Create domino object
pbmc_dom <- create_domino(
    rl_map = rl_map,             # receptor-ligand map data frame
    features = auc_in,           # TF scores (AUC matrix)
    counts = counts,             # counts matrix
    z_scores = z_scores,         # scaled expression data
    clusters = clusters,         # vector of cell cluster assignments
    tf_targets = regulon_list,   # list of TFs and their regulons
    use_clusters = TRUE,         # assess receptor activation and ligand expression on a per-cluster basis
    use_complexes = TRUE,        # include receptors and genes that function as a complex in results
    remove_rec_dropout = FALSE   # whether to remove zeroes from correlation calculations
)

# Build domino object with parameters
pbmc_dom <- build_domino(
    dom = pbmc_dom,
    min_tf_pval = .001, # Threshold for p-value of DE for TFs
    max_tf_per_clust = 25,
    max_rec_per_tf = 25,
    rec_tf_cor_threshold = .25, # Minimum correlation between receptor and TF
    min_rec_percentage = 0.1 # Minimum percent of cells that must express receptor
)

# Save domino object for generating test data
saveRDS(pbmc_dom, file="pbmc_domino_built.rds")
