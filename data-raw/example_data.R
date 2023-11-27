# Code for preparing example objects for use in vignettes
# Accessible to users for exploration as well
library(Seurat)
library(domino2)

# Zenodo host of outputs from SCENIC analysis
data_url <- "https://zenodo.org/records/10161143/files"
temp_dir <- tempdir()

# install 10X Genomics PBMC3K data
pbmc_url <- "https://cf.10xgenomics.com/samples/cell/pbmc3k/pbmc3k_filtered_gene_bc_matrices.tar.gz"
pbmc_tar <- paste0(temp_dir, "/pbmc3k_filtered_gene_bc_matrices.tar.gz")
download.file(url = pbmc_url, destfile = pbmc_tar)
pbmc_dir <- paste0(temp_dir, "/pbmc3k_filtered_gene_bc_matrices")
untar(tarfile = pbmc_tar, exdir = pbmc_dir)
pbmc.data <- Read10X(data.dir = paste0(pbmc_dir, "/filtered_gene_bc_matrices/hg19"))

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

# Save Seurat object for generating test data
saveRDS(pbmc, "inst/extdata/pbmc_seurat.rds")

# Load SCENIC outputs
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

# Save domino object for generating test data
saveRDS(pbmc_dom, "inst/extdata/pbmc_domino_built.rds")
