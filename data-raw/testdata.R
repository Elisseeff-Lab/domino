# generate objects for comparison in testing scripts

library(domino2, lib = "../domino_libraries/libraries/v0_2_1/")

# code to prepare `pbmc3k` dataset
pbmc.data <- Seurat::Read10X(data.dir = "../inst/extdata/pbmc3k_filtered_gene_bc_matrices/")
pbmc <- Seurat::CreateSeuratObject(counts = pbmc.data, project = "pbmc3k", min.cells = 3, min.features = 200)
pbmc[["percent.mt"]] <- Seurat::PercentageFeatureSet(pbmc, pattern = "^MT-")
pbmc <- Seurat::NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)
pbmc <- Seurat::FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)
pbmc <- Seurat::ScaleData(pbmc, features = rownames(pbmc))
pbmc <- Seurat::RunPCA(pbmc, features = Seurat::VariableFeatures(object = pbmc), verbose = FALSE)
pbmc <- Seurat::FindNeighbors(pbmc, dims = 1:10)
pbmc <- Seurat::FindClusters(pbmc, resolution = 0.5)
pbmc <- Seurat::RunUMAP(pbmc, dims = 1:10)
# Annotate clusters with cell phenotypes
cell_dict <- data.frame(cluster = c(0:8), cell_type = c("naive_CD4_T_cell", "CD14_monocyte", "memory_CD4_T_cell", "B_cell", "CD8_T_cell", "CD16_monocyte", "NK_cell", "dendritic_cell", "platelet"))
pbmc$cell_type <- plyr::mapvalues(pbmc$seurat_clusters, from = cell_dict$cluster, to = cell_dict$cell_type)

# subset the pbmc data to fewer cells to meet package requirements
RNA_features <- c("FAS", "FASLG", "ITGAM", "ITGB2", "FCER2", "LCK", "CD8A", "CD8B", "C5", "C5AR1")
TF_features <- c("CLOCK", "TAF1", "BCL6")
name_features <- c("FAS", "FASLG", "integrin_aMb2_complex", "FCER2", "LCK", "CD8_receptor", "C5", "C5AR1")
cell_types_dwn <- c("CD8_T_cell", "NK_cell", "B_cell")
n <- 100
cell_list <- list()
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

# cellphoneDB database
cellphonedb_4_path <- "../inst/extdata/cellphoneDB_v4/"
complexes <- read.csv(paste0(cellphonedb_4_path, "/complex_input.csv"), stringsAsFactors = FALSE)
complexes_test <- complexes[complexes$complex_name %in% name_features,]
genes <- read.csv(paste0(cellphonedb_4_path, "/gene_input.csv"), stringsAsFactors = FALSE)
genes_test <- genes[genes$gene_name %in% RNA_features,]
proteins <- read.csv(paste0(cellphonedb_4_path, "/protein_input.csv"), stringsAsFactors = FALSE)
proteins_test <- proteins[proteins$uniprot %in% genes_test$uniprot,]
interactions <- read.csv(paste0(cellphonedb_4_path, "/interaction_input.csv"), stringsAsFactors = FALSE)
interactions_test <- interactions[
  (interactions$partner_a %in% proteins_test$uniprot | interactions$partner_a %in% complexes_test$complex_name) & 
    (interactions$partner_b%in% proteins_test$uniprot | interactions$partner_b %in% complexes_test$complex_name),]

# read in and subset SCENIC inputs
auc <- read.table(
  "../inst/extdata/scenic/auc_pbmc_3k.csv",
  header = TRUE, row.names = 1, stringsAsFactors = FALSE, sep = ","
)
auc <- t(auc)
rownames(auc) <- gsub("\\.\\.\\.$", "", rownames(auc))
auc_test <- auc[TF_features, barcodes_dwn]

regulons <- read.csv("../inst/extdata/scenic/regulons_pbmc_3k.csv")
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