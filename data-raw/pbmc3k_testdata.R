set.seed(123)

# code to prepare `pbmc3k` dataset
pbmc.data <- Seurat::Read10X(data.dir = "inst/extdata/pbmc3k_filtered_gene_bc_matrices/")
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
# save expression features used by domino
save(clusters_test, file = "data/cluster_test.rda")
save(RNA_count_test, file = "data/RNA_count_test.rda")
save(RNA_zscore_test, file = "data/RNA_zscore_test.rda")

# cellphoneDB database
cellphonedb_4_path <- "inst/extdata/cellphoneDB_v4/"
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
# save CellPhoneDB inputs
save(complexes_test, file = "data/complexes_test.rda")
save(genes_test, file = "data/genes_test.rda")
save(proteins_test, file = "data/proteins_test.rda")
save(interactions_test, file = "data/interactions_test.rda")

# read in and subset SCENIC inputs
auc <- read.table(
  "inst/extdata/auc_pbmc_3k.csv",
  header = TRUE, row.names = 1, stringsAsFactors = FALSE, sep = ","
)
auc <- t(auc)
rownames(auc) <- gsub("\\.\\.\\.$", "", rownames(auc))
auc_test <- auc[TF_features, barcodes_dwn]
save(auc_test, file="data/auc_test.rda")

regulons <- read.csv("inst/extdata/regulons_pbmc_3k.csv")
regulons <- regulons[-1:-2, ]
colnames(regulons) <- c("TF", "MotifID", "AUC", "NES", "MotifSimilarityQvalue", "OrthologousIdentity", "Annotation", "Context", "TargetGenes", "RankAtMax")
regulons_test <- regulons[regulons$TF %in% TF_features,]
save(regulons_test, file="data/regulons_test.rda")
