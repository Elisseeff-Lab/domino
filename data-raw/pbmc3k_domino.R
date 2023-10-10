## code to prepare `pbmc_dom` dataset goes here

# Load pbmc data:
load("../data/pbmc.rda")

# Load SCENIC outputs
scenic_dir = "../scenic/"
regulons <- read.csv(paste0(scenic_dir, "/pbmc_regulons.csv"))
auc <- read.table(paste0(scenic_dir, "/pbmc_auc.csv"),
                  header = TRUE, row.names = 1,
                  stringsAsFactors = FALSE, sep = ",")

# Create list of SCENIC regulons
regulons <- regulons[-1:-2,]
colnames(regulons) <- c("TF", "MotifID", "AUC", "NES", "MotifSimilarityQvalue", "OrthologousIdentity", "Annotation", "Context", "TargetGenes", "RankAtMax")
regulon_list <- create_regulon_list_scenic(regulons = regulons)

# Load AUC and adjust names
auc_in <- as.data.frame(t(auc))
# Remove pattern "..." from the end of all rownames:
rownames(auc_in) <- gsub("\\.\\.\\.$", "", rownames(auc_in))

# Load CellPhoneDB Database
cellphonedb_2_path <- "../cellphonedb"
complexes <- read.csv(paste0(cellphonedb_2_path, "/complex_input.csv"), stringsAsFactors = FALSE)
genes <- read.csv(paste0(cellphonedb_2_path, "/gene_input.csv"), stringsAsFactors = FALSE)
interactions <- read.csv(paste0(cellphonedb_2_path, "/interaction_input.csv"), stringsAsFactors = FALSE)
proteins <- read.csv(paste0(cellphonedb_2_path, "/protein_input.csv"), stringsAsFactors = FALSE)

# Make RL Map
rl_map <- create_rl_map_cellphonedb(
  genes = genes, 
  proteins = proteins, 
  interactions = interactions,
  complexes = complexes,
  database_name = "CellPhoneDB_v4.0" # database version used
)

# Prepare pbmc data
counts <- pbmc@assays$RNA@counts
z_scores <- as.matrix(pbmc@assays$RNA@scale.data)
clusters <- as.factor(pbmc$cell_type)

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

# Save
usethis::use_data(pbmc_dom, overwrite = TRUE)
