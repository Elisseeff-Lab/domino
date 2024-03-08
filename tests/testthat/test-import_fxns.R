## unit tests use internal data stored in R/sysdata.R

## large scale test that processing functions generate results matching the state of the code in v0.2.1
test_that("creation of rl_map from CellPhoneDB v4 input", {
  # use rl_map created in domino2 v0.2.1
  rl_map <- create_rl_map_cellphonedb(
    genes = genes_tiny,
    proteins = proteins_tiny,
    interactions = interactions_tiny,
    complexes = complexes_tiny
  )
  expect_equal(rl_map, rl_map_tiny)
})

test_that("formating of SCENIC regulons output as a list", {
  # use regulon list created in domino2 v0.2.1
  regulon_ls <- create_regulon_list_scenic(regulons = regulons_tiny)
  expect_equal(regulon_ls, regulon_list_tiny)
})

test_that("creation of a domino object from SCENIC and CellPhoneDB inputs", {
  # domino object created in domino2 v0.2.1
  # rl_map created in domino2 v0.2.1 from CellPhoneDB v4
  # scenic inputs from pySCENIC v0.11.0
  # minimal expression data from 360 cells in 3 cell types,
  # expression features for 16 genes
  # 360 cell barcodes have annotations of cell type assignment in clusters_tiny
  pbmc_dom <- create_domino(
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
  expect_equal(pbmc_dom@linkages$complexes, pbmc_dom_tiny@linkages$complexes)
  expect_equal(pbmc_dom@linkages$rec_lig, pbmc_dom_tiny@linkages$rec_lig)
  expect_equal(pbmc_dom@linkages$tf_targets, pbmc_dom_tiny@linkages$tf_targets)
  expect_equal(pbmc_dom@clust_de, pbmc_dom_tiny@clust_de)
  expect_equal(pbmc_dom@cor, pbmc_dom_tiny@cor)
})

test_that("building a domino object under set parameters", {
  # built domino object created in domino2 v0.2.1
  # domino object created in domino2 v0.2.1
  pbmc_dom_built <- build_domino(
    dom = pbmc_dom_tiny,
    min_tf_pval = .05,
    max_tf_per_clust = Inf,
    max_rec_per_tf = Inf,
    rec_tf_cor_threshold = .1,
    min_rec_percentage = 0.01
  )
  expect_equal(pbmc_dom_built@linkages$clust_tf,
               pbmc_dom_built_tiny@linkages$clust_tf)

  expect_equal(pbmc_dom_built@linkages$clust_rec,
               pbmc_dom_built_tiny@linkages$clust_rec)

  expect_equal(pbmc_dom_built@linkages$clust_tf_rec,
               pbmc_dom_built_tiny@linkages$clust_tf_rec)

  expect_equal(pbmc_dom_built@linkages$clust_incoming_lig,
               pbmc_dom_built_tiny@linkages$clust_incoming_lig)

  expect_equal(pbmc_dom_built@linkages$tf_rec,
               pbmc_dom_built_tiny@linkages$tf_rec)

  expect_equal(pbmc_dom_built@cl_signaling_matrices,
               pbmc_dom_built_tiny@cl_signaling_matrices)

  expect_equal(pbmc_dom_built@signaling,
               pbmc_dom_built_tiny@signaling)
})

test_that("create_rl_map_cellphonedb fails on wrong input arg type.", {

  expect_error(create_rl_map_cellphonedb(
    genes = list(), proteins = proteins_tiny,
    interactions = interactions_tiny, complexes = complexes_tiny
  ))

  expect_error(create_rl_map_cellphonedb(
    genes = genes_tiny, proteins = list(),
    interactions = interactions_tiny, complexes = complexes_tiny
  ))

  expect_error(create_rl_map_cellphonedb(
    genes = genes_tiny, proteins = proteins_tiny,
    interactions = list(), complexes = complexes_tiny
  ))

  expect_error(create_rl_map_cellphonedb(
    genes = genes_tiny, proteins = proteins_tiny,
    interactions = interactions_tiny, complexes = list()
  ))

  expect_error(create_rl_map_cellphonedb(
    genes = genes_tiny, proteins = proteins_tiny,
    interactions = interactions_tiny, complexes = complexes_tiny,
    database_name = list()
  ))

  expect_error(create_rl_map_cellphonedb(
    genes = genes_tiny, proteins = proteins_tiny,
    interactions = interactions_tiny, complexes = complexes_tiny,
    database_name = c("length", ">1")
  ))
})

test_that("create_rl_map_cellphonedb fails on wrong input arg type.", {
  #bad rl map
  bad_rl_map <- "rl_map"
  expect_error(create_domino(bad_rl_map,
                             features_tiny,
                             counts = RNA_count_tiny,
                             z_scores = RNA_zscore_tiny,
                             clusters = clusters_tiny),
  "rl_map must be a data.frame with column names gene_A, gene_B, type_A, and type_B")

  bad_rl_map <- rl_map_tiny
  colnames(bad_rl_map) <- paste(colnames(bad_rl_map), "qq")
  expect_error(create_domino(bad_rl_map,
                             features_tiny,
                             counts = RNA_count_tiny,
                             z_scores = RNA_zscore_tiny,
                             clusters = clusters_tiny),
  "rl_map must be a data.frame with column names gene_A, gene_B, type_A, and type_B")

  #bad features
  bad_features <- matrix()
  expect_error(create_domino(rl_map_tiny,
                             bad_features,
                             counts = RNA_count_tiny,
                             z_scores = RNA_zscore_tiny,
                             clusters = clusters_tiny),
  "features must be either a file path or a named matrix with cells as columns and features as rows")

  #seurat or counts, zscores and clusters
  expect_error(create_domino(rl_map_tiny,
                             auc_tiny),
  "Either a Seurat object OR counts, z scores, and clusters must be provided")

  expect_error(create_domino(rl_map_tiny,
                             auc_tiny,
                             counts = RNA_count_tiny,
                             z_scores = RNA_zscore_tiny),
  "Either a Seurat object OR counts, z scores, and clusters must be provided")

  #bad rec_min threshold
  expect_error(create_domino(rl_map_tiny,
                             auc_tiny,
                             counts = RNA_count_tiny,
                             z_scores = RNA_zscore_tiny,
                             clusters = clusters_tiny,
                             rec_min_thresh = 20
                             ),
  "rec_min_thresh must be a number between 0 and 1")

  expect_error(create_domino(rl_map_tiny,
                             auc_tiny,
                             counts = RNA_count_tiny,
                             z_scores = RNA_zscore_tiny,
                             clusters = clusters_tiny,
                             rec_min_thresh = -20
                             ),
  "rec_min_thresh must be a number between 0 and 1")

  expect_error(create_domino(rl_map_tiny,
                             auc_tiny,
                             counts = RNA_count_tiny,
                             z_scores = RNA_zscore_tiny,
                             clusters = clusters_tiny,
                             tf_selection_method = "non-existent"
                             ),
  "tf_selection_method must be one of all, clusters, or variable")
})