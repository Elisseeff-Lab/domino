# generate objects for comparison in testing scripts

library(SingleCellExperiment)
library(dominoSignal)

# load data for generation of test results from zenodo repository
# Zenodo host of outputs from SCENIC analysis
data_url <- "https://zenodo.org/records/10951634/files"
temp_dir <- tempdir()

pbmc_dir <- paste0(temp_dir, "/pbmc")
if (!dir.exists(pbmc_dir)) {
  dir.create(pbmc_dir)
}

# SingleCellExperiment object of preprocessed PBMC3K data
download.file(url = paste0(data_url, "/pbmc3k_sce.rds"),
              destfile = paste0(pbmc_dir, "/pbmc3k_sce.rds"))
pbmc <- readRDS(paste0(pbmc_dir, "/pbmc3k_sce.rds"))

# SCENIC input files
scenic_dir <- paste0(temp_dir, "/scenic")
if (!dir.exists(scenic_dir)) {
  dir.create(scenic_dir)
}
download.file(url = paste0(data_url, "/auc_pbmc_3k.csv"),
              destfile = paste0(scenic_dir, "/auc_pbmc_3k.csv"))
download.file(url = paste0(data_url, "/regulons_pbmc_3k.csv"),
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
  "CCL20", "CXCR3", "CCR6",
  "IL7", "IL7R", "IL2RG",
  "TGFB3", "TGFBR3",
  "ITGA6", "ITGB4", "NRG1"
  )
TF_features <- c(
  "DBP",
  "FLI1", "ZNF431",
  "ZNF324", "CREM", "FOSL1"
)
name_features <- c(
  "CCL20", "CXCR3", "CCR6",
  "IL7", "IL7_receptor",
  "TGFB3", "TGFBR3",
  "integrin_a6b4_complex", "NRG1"
  )
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
names(clusters_tiny) <- barcodes_dwn

counts <- assay(pbmc, "counts")
z_scores <- t(scale(t(assay(pbmc, "logcounts"))))

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

# create a second tiny domino object for list comparison functions
clusters_tiny_alt <- clusters_tiny[c(121:240, 1:120, 241:360)]
names(clusters_tiny_alt) <- names(clusters_tiny)
pbmc_dom_tiny_alt <- create_domino(
  rl_map = rl_map_tiny,
  features = auc_tiny,
  counts = RNA_count_tiny,
  z_scores = RNA_zscore_tiny,
  # reassigned clusters
  clusters = clusters_tiny_alt,
  tf_targets = regulon_list_tiny,
  use_clusters = TRUE,
  use_complexes = TRUE,
  remove_rec_dropout = FALSE
)
pbmc_dom_built_tiny_alt <- build_domino(
  dom = pbmc_dom_tiny_alt,
  min_tf_pval = .05,
  max_tf_per_clust = Inf,
  max_rec_per_tf = Inf,
  rec_tf_cor_threshold = .1,
  min_rec_percentage = 0.01
)
dom_ls_tiny <- list(
  dom1 = pbmc_dom_built_tiny,
  dom2 = pbmc_dom_built_tiny_alt
)

# example linkage summary for comparitive functions
linkage_sum_tiny <- new("linkage_summary",
  subject_meta = data.frame(
    "subject_names" = paste0("P",1:6),
    "group" = c(rep("G1", 3), rep("G2", 3))
  ), 
  subject_names = factor(
    paste0("P",1:6), levels = paste0("P",1:6)
  ), 
  subject_linkages = list(
    "P1" = list(
      "C1" = list(
        "tfs" = c("TF1", "TF2", "TF3", "TF4"),
        "rec" = c("R1", "R2", "R3", "R4"),
        "incoming_lig" = c("L1", "L2", "L3", "L4"),
        "tfs_rec" = c("TF1 <- R1", "TF2 <- R2", "TF3 <- R3", "TF4 <- R4"),
        "rec_lig" = c("R1 <- L1", "R2 <- L2", "R3 <- L3", "R4 <- L4")
      ),
      "C2" = list(
        "tfs" = c("TF2", "TF3", "TF4"),
        "rec" = c("R2", "R3", "R4"),
        "incoming_lig" = c("L2", "L3", "L4"),
        "tfs_rec" = c("TF2 <- R2", "TF3 <- R3", "TF4 <- R4"),
        "rec_lig" = c("R2 <- L2", "R3 <- L3", "R4 <- L4")
      )
    ),
    "P2" = list(
      "C1" = list(
        "tfs" = c("TF1", "TF2"),
        "rec" = c("R1", "R2"),
        "incoming_lig" = c("L1", "L2"),
        "tfs_rec" = c("TF1 <- R1", "TF2 <- R2"),
        "rec_lig" = c("R1 <- L1", "R2 <- L2")
      ),
      "C2" = list(
        "tfs" = c("TF3", "TF4"),
        "rec" = c("R3", "R4"),
        "incoming_lig" = c("L3", "L4"),
        "tfs_rec" = c("TF3 <- R3", "TF4 <- R4"),
        "rec_lig" = c("R3 <- L3", "R4 <- L4")
      )
    ),
    "P3" = list(
      "C1" = list(
        "tfs" = c("TF1", "TF2"),
        "rec" = c("R1", "R2"),
        "incoming_lig" = c("L1", "L2"),
        "tfs_rec" = c("TF1 <- R1", "TF2 <- R2"),
        "rec_lig" = c("R1 <- L1", "R2 <- L2")
      ),
      "C2" = list(
        "tfs" = c("TF3"),
        "rec" = c("R3"),
        "incoming_lig" = c("L3"),
        "tfs_rec" = c("TF3 <- R3"),
        "rec_lig" = c("R3 <- L3")
      )
    ),
    "P4" = list(
      "C1" = list(
        "tfs" = c("TF2", "TF3", "TF4"),
        "rec" = c("R2", "R3", "R4"),
        "incoming_lig" = c("L2", "L3", "L4"),
        "tfs_rec" = c("TF2 <- R2", "TF3 <- R3", "TF4 <- R4"),
        "rec_lig" = c("R2 <- L2", "R3 <- L3", "R4 <- L4")
      ),
      "C2" = list(
        "tfs" = c("TF1", "TF2", "TF3", "TF4"),
        "rec" = c("R1", "R2", "R3", "R4"),
        "incoming_lig" = c("L1", "L2", "L3", "L4"),
        "tfs_rec" = c("TF1 <- R1", "TF2 <- R2", "TF3 <- R3", "TF4 <- R4"),
        "rec_lig" = c("R1 <- L1", "R2 <- L2", "R3 <- L3", "R4 <- L4")
      )
    ),
    "P5" = list(
      "C1" = list(
        "tfs" = c("TF3"),
        "rec" = c("R3"),
        "incoming_lig" = c("L3"),
        "tfs_rec" = c("TF3 <- R3"),
        "rec_lig" = c("R3 <- L3")
      ),
      "C2" = list(
        "tfs" = c("TF1", "TF2", "TF3", "TF4"),
        "rec" = c("R1", "R2", "R3", "R4"),
        "incoming_lig" = c("L1", "L2", "L3", "L4"),
        "tfs_rec" = c("TF1 <- R1", "TF2 <- R2", "TF3 <- R3", "TF4 <- R4"),
        "rec_lig" = c("R1 <- L1", "R2 <- L2", "R3 <- L3", "R4 <- L4")
      )
    ),
    "P6" = list(
      "C1" = list(
        "tfs" = c(),
        "rec" = c(),
        "incoming_lig" = c(),
        "tfs_rec" = c(),
        "rec_lig" = c()
      ),
      "C2" = list(
        "tfs" = c("TF1", "TF2", "TF3"),
        "rec" = c("R1", "R2", "R3"),
        "incoming_lig" = c("L1", "L2", "L3"),
        "tfs_rec" = c("TF1 <- R1", "TF2 <- R2", "TF3 <- R3"),
        "rec_lig" = c("R1 <- L1", "R2 <- L2", "R3 <- L3")
      )
    )
  )
)

# test result object from differential linkage tests
tiny_differential_linkage_c1 <- test_differential_linkages(
  linkage_summary = linkage_sum_tiny, cluster = "C1", group.by = "group", 
  linkage = "rec", subject_names = linkage_sum_tiny@subject_names, test_name = "fishers.exact"
)
tiny_differential_linkage_c2 <- test_differential_linkages(
  linkage_summary = linkage_sum_tiny, cluster = "C2", group.by = "group", 
  linkage = "rec", subject_names = linkage_sum_tiny@subject_names, test_name = "fishers.exact"
)

# Save all test files to internal sysdata object
usethis::use_data(
  pbmc_dom_built_tiny, complexes_tiny, genes_tiny, proteins_tiny, interactions_tiny,
  pbmc_dom_tiny, regulon_list_tiny, rl_map_tiny, regulons_tiny, clusters_tiny,
  RNA_count_tiny, RNA_zscore_tiny, auc_tiny,
  dom_ls_tiny, linkage_sum_tiny, tiny_differential_linkage_c1, tiny_differential_linkage_c2
)
