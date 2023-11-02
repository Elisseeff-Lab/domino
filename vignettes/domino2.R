## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  echo = TRUE,
  message = FALSE,
  warning = FALSE,
  fig.cap="",
  tidy = 'styler'
)


## ----libraries, include=FALSE-------------------------------------------------
library(domino2)
library(Seurat)
library(loomR)
library(plyr)
library(circlize)
library(ComplexHeatmap)
library(knitr)

set.seed(123)

## ----eval=TRUE, message = FALSE-----------------------------------------------
# Preprocessing of PBMC 3K tutorial data set using Seurat functions
pbmc.data <- Read10X(data.dir = "../inst/extdata/pbmc3k_filtered_gene_bc_matrices")
pbmc <- CreateSeuratObject(counts = pbmc.data, project = "pbmc3k", min.cells = 3, min.features = 200)
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")
pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)
pbmc <- ScaleData(pbmc, features = rownames(pbmc))
pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))
pbmc <- FindNeighbors(pbmc, dims = 1:10)
pbmc <- FindClusters(pbmc, resolution = 0.5)
pbmc <- RunUMAP(pbmc, dims = 1:10)
# Annotate clusters with cell phenotypes
cell_dict <- data.frame(
  cluster = c(0:8),
  cell_type = c("naive_CD4_T_cell", "CD14_monocyte", "memory_CD4_T_cell", "B_cell", "CD8_T_cell", 
    "CD16_monocyte", "NK_cell", "dendritic_cell", "platelet")
)
pbmc$cell_type <- 
  plyr::mapvalues(
    pbmc$seurat_clusters,
    from = cell_dict$cluster,
    to = cell_dict$cell_type
  )

## ----eval = FALSE-------------------------------------------------------------
#  # save Seurat object as RDS
#  save(pbmc, file = "../data/pbmc.rda")
#  
#  # save loom counts matrix
#  pbmc_counts <- pbmc@assays$RNA@counts
#  pbmc_loom <- loomR::create(filename = "../inst/extdata/pbmc3k_counts.loom", data = pbmc_counts)
#  pbmc_loom$close_all() # Remember to manually close connection to loom files!

## ----Save Counts Matrix tsv File, eval=FALSE----------------------------------
#  pbmc_counts <- pbmc@assays$RNA@counts
#  write.table(t(as.matrix(pbmc_counts)), "../data/pbmc3k_counts.tsv",
#      sep = "\t", col.names = NA)

## ----eval = FALSE-------------------------------------------------------------
#  if(!require(remotes)){
#      install.packages('remotes')
#  }
#  remotes::install_github('FertigLab/domino2@v0.2.2')

## ----Load SCENIC Results------------------------------------------------------
scenic_dir = "../inst/extdata/scenic/"
regulons <- read.csv(paste0(scenic_dir, "/regulons_pbmc_3k.csv"))
auc <- read.table(paste0(scenic_dir, "/auc_pbmc_3k.csv"),
                  header = TRUE, row.names = 1,
                  stringsAsFactors = FALSE, sep = ",")

## ----Create Regulon List------------------------------------------------------
regulons <- regulons[-1:-2,]
colnames(regulons) <- c("TF", "MotifID", "AUC", "NES", "MotifSimilarityQvalue", "OrthologousIdentity", 
  "Annotation", "Context", "TargetGenes", "RankAtMax")
regulon_list <- create_regulon_list_scenic(regulons = regulons)

## ----Orient AUC Matrix--------------------------------------------------------
auc_in <- as.data.frame(t(auc))
# Remove pattern "..." from the end of all rownames:
rownames(auc_in) <- gsub("\\.\\.\\.$", "", rownames(auc_in))

## ----Load CellPhoneDB Database------------------------------------------------
cellphonedb_2_path <- "../inst/extdata/cellphoneDB_v4/"
complexes <- read.csv(paste0(cellphonedb_2_path, "/complex_input.csv"), stringsAsFactors = FALSE)
genes <- read.csv(paste0(cellphonedb_2_path, "/gene_input.csv"), stringsAsFactors = FALSE)
interactions <- read.csv(paste0(cellphonedb_2_path, "/interaction_input.csv"), stringsAsFactors = FALSE)
proteins <- read.csv(paste0(cellphonedb_2_path, "/protein_input.csv"), stringsAsFactors = FALSE)

rl_map <- create_rl_map_cellphonedb(
  genes = genes, 
  proteins = proteins, 
  interactions = interactions,
  complexes = complexes,
  database_name = "CellPhoneDB_v4.0" # database version used
)

knitr::kable(head(rl_map))

## ----Appending interactions---------------------------------------------------
# Integrin complexes are not annotated as receptors in CellPhoneDB_v4.0
# collagen-integrin interactions between cells may be missed unless tables from the CellPhoneDB reference are edited or the interactions are manually added

col_int_df <- data.frame(
  "int_pair" = "a11b1 complex & COLA1_HUMAN",
  "name_A" = "a11b1 complex", "uniprot_A" = "P05556,Q9UKX5", "gene_A" = "ITB1,ITA11", "type_A" = "R",
  "name_B" = "COLA1_HUMAN", "uniprot_B" = "P02452,P08123", "gene_B" = "COL1A1,COL1A2", "type_B" = "L",
  "annotation_strategy" = "manual", "source" = "manual", "database_name" = "manual"
)
rl_map_append <- rbind(col_int_df, rl_map)
knitr::kable(head(rl_map_append))

## ----load cell features-------------------------------------------------------
counts <- pbmc@assays$RNA@counts
z_scores <- as.matrix(pbmc@assays$RNA@scale.data)
clusters <- as.factor(pbmc$cell_type)

## ----Create Domino, results='hide', warning=FALSE-----------------------------
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

## ----Create Domino Seurat, results='hide', warning=FALSE, eval = FALSE--------
#  pbmc_dom <- create_domino(
#    rl_map = rl_map,             # receptor-ligand map data frame
#    features = auc_in,           # TF scores (AUC matrix)
#    ser = pbmc,                  # Seurat object containing counts, scaled counts, and cell cluster assignments
#    tf_targets = regulon_list,   # list of TFs and their regulons
#    use_clusters = TRUE,         # assess receptor activation and ligand expression on a per-cluster basis
#    use_complexes = TRUE,        # include receptors and genes that function as a complex in results
#    remove_rec_dropout = FALSE   # whether to remove zeroes from correlation calculations
#  )

## ----Build Domino-------------------------------------------------------------
pbmc_dom <- build_domino(
  dom = pbmc_dom,
  min_tf_pval = .001, # Threshold for p-value of DE for TFs
  max_tf_per_clust = 25,
  max_rec_per_tf = 25,
  rec_tf_cor_threshold = .25, # Minimum correlation between receptor and TF
  min_rec_percentage = 0.1 # Minimum percent of cells that must express receptor
)

## ----Build Domino All Significant Interactions, eval=FALSE--------------------
#  pbmc_dom_all <- build_domino(
#    dom = pbmc_dom,
#    min_tf_pval = .001,
#    max_tf_per_clust = Inf,
#    max_rec_per_tf = Inf,
#    rec_tf_cor_threshold = .25,
#    min_rec_percentage = 0.1
#  )

## -----------------------------------------------------------------------------
feat_heatmap(pbmc_dom, norm = TRUE, bool = FALSE)

## -----------------------------------------------------------------------------
signaling_network(pbmc_dom, edge_weight = 1, max_thresh = 5)

## -----------------------------------------------------------------------------
gene_network(pbmc_dom, clust = 'dendritic_cell', layout = 'grid')

## -----------------------------------------------------------------------------
gene_network(pbmc_dom, clust = 'dendritic_cell', OutgoingSignalingClust = "CD14_monocyte",
              layout = 'grid')

