# generate objects for comparison in testing scripts

library(Seurat)
library(domino2, lib = "../domino_libraries/libraries/v0_2_1/")

# load data for generation of test results from zenodo repository
# Zenodo host of outputs from SCENIC analysis
data_url <- "https://zenodo.org/records/10161143/files"
temp_dir <- tempdir()

pbmc_dir <- paste0(temp_dir, "/pbmc")
if (!dir.exists(pbmc_dir)) {
  dir.create(pbmc_dir)
}
# Seurat object of preprocessed PBMC3K data
download.file(url = paste0(data_url, "/pbmc_seurat.rds"),
              destfile = paste0(pbmc_dir, "/pbmc_seurat.rds"))
pbmc <- readRDS(paste0(pbmc_dir, "/pbmc_seurat.rds"))

# SCENIC input files
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

# CellPhoneDB Database
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

# subset the pbmc data to fewer cells to meet package requirements
RNA_features <- c("FAS", "FASLG", "ITGAM", "ITGB2", "FCER2", "LCK", "CD8A", "CD8B", "C5", "C5AR1")
TF_features <- c("CLOCK", "TAF1", "BCL6")
name_features <- c("FAS", "FASLG", "integrin_aMb2_complex", "FCER2", "LCK", "CD8_receptor", "C5", "C5AR1")
cell_types_dwn <- c("CD8_T_cell", "NK_cell", "B_cell")
n <- 120
cell_list <- list()
set.seed(123)
for(i in seq_along(cell_types_dwn)){
  cell <- cell_types_dwn[i]
  cell_barcodes <- colnames(pbmc)[pbmc$cell_type == cell]
  dwn_barcodes <- sample(cell_barcodes, n, replace = FALSE)
  cell_list[[cell]] <- dwn_barcodes
}
barcodes_dwn <- unlist(cell_list)
cluster_dwn <- factor(
    pbmc$cell_type[barcodes_dwn],
    levels = cell_types_dwn
)
clusters_test <- cluster_dwn
RNA_count_test <- pbmc@assays$RNA@counts[
  rownames(pbmc@assays$RNA@counts) %in% RNA_features, 
  colnames(pbmc@assays$RNA@counts) %in% barcodes_dwn]
RNA_zscore_test <- pbmc@assays$RNA@scale.data[
  rownames(pbmc@assays$RNA@scale.data) %in% RNA_features, 
  colnames(pbmc@assays$RNA@scale.data) %in% barcodes_dwn]

# subset CellPhoneDB inputs
complexes_test <- complexes[complexes$complex_name %in% name_features,]
genes_test <- genes[genes$gene_name %in% RNA_features,]
proteins_test <- proteins[proteins$uniprot %in% genes_test$uniprot,]
interactions_test <- interactions[
  (interactions$partner_a %in% proteins_test$uniprot | interactions$partner_a %in% complexes_test$complex_name) & 
    (interactions$partner_b%in% proteins_test$uniprot | interactions$partner_b %in% complexes_test$complex_name),]


# subset SCENIC inputs
auc <- t(auc)
rownames(auc) <- gsub("\\.\\.\\.$", "", rownames(auc))
auc_test <- auc[TF_features, barcodes_dwn]

regulons <- regulons[-1:-2, ]
colnames(regulons) <- c("TF", "MotifID", "AUC", "NES", "MotifSimilarityQvalue", "OrthologousIdentity", "Annotation", "Context", "TargetGenes", "RankAtMax")
regulons_test <- regulons[regulons$TF %in% TF_features,]

# Make rl_map
rl_map_test <- domino2::create_rl_map_cellphonedb(
  genes = genes_test,
  proteins = proteins_test,
  interactions = interactions_test,
  complexes = complexes_test
)

# Get regulon list
regulon_list_test <- domino2::create_regulon_list_scenic(
  regulons = regulons_test
)

# Create test domino object
pbmc_dom_test <- create_domino(
  rl_map = rl_map_test,
  features = auc_test,
  counts = RNA_count_test,
  z_scores = RNA_zscore_test,
  clusters = clusters_test,
  tf_targets = regulon_list_test,
  use_clusters = TRUE,
  use_complexes = TRUE,
  remove_rec_dropout = FALSE
)

# Create built domino object
pbmc_dom_built_test <- build_domino(
  dom = pbmc_dom_test,
  min_tf_pval = .05,
  max_tf_per_clust = Inf,
  max_rec_per_tf = Inf,
  rec_tf_cor_threshold = .15,
  min_rec_percentage = 0.03 
)

# Save all test files to internal sysdata object
usethis::use_data(pbmc_dom_built_test, complexes_test, genes_test, proteins_test, interactions_test,
    pbmc_dom_test, regulon_list_test, rl_map_test, regulons_test, clusters_test,
    RNA_count_test, RNA_zscore_test, auc_test, internal = TRUE)