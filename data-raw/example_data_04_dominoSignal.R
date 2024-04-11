# example_data_04_dominoSignal

library(DropletUtils)
library(EnsDb.Hsapiens.v86)
library(scater)
library(scran)
library(dominoSignal)

temp_dir <- tempdir()

# obtain expression data
sce <- readRDS("pbmc3k_sce.rds")

# obtain SCENIC data

scenic_dir <- "scenic"
if (!dir.exists(scenic_dir)) {
  stop("there is no scenic dir, have you ran example_data_04 script?")
}
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
counts <- assay(sce, "counts")
z_scores <- t(scale(t(assay(sce, "logcounts"))))
clusters <- factor(sce$cell_type)
names(clusters) <- colnames(sce)

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
