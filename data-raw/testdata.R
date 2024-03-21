# generate objects for comparison in testing scripts

library(SingleCellExperiment)
library(dominoSignal)

# load data for generation of test results from zenodo repository
# Zenodo host of outputs from SCENIC analysis
data_url <- "https://zenodo.org/records/10850479/files"
temp_dir <- tempdir()

pbmc_dir <- paste0(temp_dir, "/pbmc")
if (!dir.exists(pbmc_dir)) {
  dir.create(pbmc_dir)
}

# SingleCellExperiment object of preprocessed PBMC3K data
download.file(url = paste0(data_url, "/pbmc_3k_sce.rds"),
              destfile = paste0(pbmc_dir, "/pbmc_3k_sce.rds"))
pbmc <- readRDS(paste0(pbmc_dir, "/pbmc_3k_sce.rds"))

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
RNA_features <- c(
  "TNF", "TNFSF13", "FASLG", "FAS", 
  "ITGAM", "ITGB2", "ITGAV", "ITGB3", "FCER2", 
  "IL7", "IL7R", "IL2RG",
  "CD22", "PTPRC",
  "IGF1", "IGF1R")
TF_features <- c("ZNF257", "ATF4", "RUNX1")
name_features <- c(
  "TNF", "TNFSF13", "FASLG", "FAS", 
  "integrin_aMb2_complex", "integrin_aVb3_complex", "FCER2", 
  "IL7", "IL7_receptor", 
  "CD22", "PTPRC",
  "IGF1", "IGF1R")
cell_types_dwn <- c("CD8_T_cell", "CD14_monocyte", "B_cell")
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
clusters_tiny <- factor(rep(names(cell_list), lengths(cell_list)))
names(clusters_tiny) <- clusters_tiny

counts <- assay(pbmc, "counts")
z_scores <-  t(scale(t(assay(pbmc, "logcounts"))))

RNA_count_tiny <- counts[rownames(assay(pbmc, "counts")) %in% RNA_features,
                         colnames(assay(pbmc, "counts")) %in% barcodes_dwn]
RNA_zscore_tiny <- z_scores[rownames(assay(pbmc, "logcounts")) %in% RNA_features,
                            colnames(assay(pbmc, "logcounts")) %in% barcodes_dwn]

# subset CellPhoneDB inputs
complexes_tiny <- complexes[complexes$complex_name %in% name_features,]
genes_tiny <- genes[genes$gene_name %in% RNA_features,]
proteins_tiny <- proteins[proteins$uniprot %in% genes_tiny$uniprot,]
interactions_tiny <- interactions[
  (interactions$partner_a %in% proteins_tiny$uniprot | interactions$partner_a %in% complexes_tiny$complex_name) & 
    (interactions$partner_b%in% proteins_tiny$uniprot | interactions$partner_b %in% complexes_tiny$complex_name),]


# subset SCENIC inputs
auc <- t(auc)
rownames(auc) <- gsub("\\.\\.\\.$", "", rownames(auc))
auc_tiny <- auc[TF_features, barcodes_dwn]

regulons <- regulons[-1:-2, ]
colnames(regulons) <- c("TF", "MotifID", "AUC", "NES", "MotifSimilarityQvalue", "OrthologousIdentity", "Annotation", "Context", "TargetGenes", "RankAtMax")
regulons_tiny <- regulons[regulons$TF %in% TF_features,]

# Make rl_map
rl_map_tiny <- dominoSignal::create_rl_map_cellphonedb(
  genes = genes_tiny,
  proteins = proteins_tiny,
  interactions = interactions_tiny,
  complexes = complexes_tiny
)

# Get regulon list
regulon_list_tiny <- dominoSignal::create_regulon_list_scenic(
  regulons = regulons_tiny
)

# Create test domino object
pbmc_dom_tiny <- create_domino(
  rl_map = rl_map_tiny,
  features = auc_tiny,
  counts = RNA_count_tiny,
  z_scores = RNA_zscore_tiny,
  clusters = clusters_tiny,
  tf_targets = regulon_list_tiny,
  use_clusters = TRUE,
  use_complexes = TRUE,
  remove_rec_dropout = FALSE
)

# Create built domino object
pbmc_dom_built_tiny <- build_domino(
  dom = pbmc_dom_tiny,
  min_tf_pval = .05,
  max_tf_per_clust = Inf,
  max_rec_per_tf = Inf,
  rec_tf_cor_threshold = .1,
  min_rec_percentage = 0.01
)

# Save all test files to internal sysdata object
usethis::use_data(pbmc_dom_built_tiny, complexes_tiny, genes_tiny, proteins_tiny, interactions_tiny,
    pbmc_dom_tiny, regulon_list_tiny, rl_map_tiny, regulons_tiny, clusters_tiny,
    RNA_count_tiny, RNA_zscore_tiny, auc_tiny, internal = TRUE)
