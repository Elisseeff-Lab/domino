# Code for preparing example objects for use in vignettes
# Accessible to users for exploration as well
library(Seurat)
pbmc.data <- Read10X(data.dir = "../inst/extdata/pbmc3k_filtered_gene_bc_matrices")
pbmc <- CreateSeuratObject(counts = pbmc.data, project = "pbmc3k", min.cells = 3, min.features = 200)
pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)
pbmc <- ScaleData(pbmc, features = rownames(pbmc))
pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))
pbmc <- FindNeighbors(pbmc, dims = 1:10)
pbmc <- FindClusters(pbmc, resolution = 0.5)

# Annotate clusters with cell phenotypes
cell_dict <- data.frame(
    cluster = c(0:8),
    cell_type = c("naive_CD4_T_cell", "CD14_monocyte", "memory_CD4_T_cell", "B_cell", "CD8_T_cell", "CD16_monocyte", "NK_cell", "dendritic_cell", "platelet")
)

pbmc$cell_type <- 
    plyr::mapvalues(
    pbmc$seurat_clusters,
    from = cell_dict$cluster,
    to = cell_dict$cell_type
    )

# Load SCENIC outputs
scenic_dir = "../inst/extdata/scenic/"
regulons <- read.csv(paste0(scenic_dir, "/regulons_pbmc_3k.csv"))
auc <- read.table(paste0(scenic_dir, "/auc_pbmc_3k.csv"),
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
cellphonedb_2_path <- "../inst/extdata/cellphoneDB_v4"
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
usethis::use_data(pbmc_dom, compress = "xz")