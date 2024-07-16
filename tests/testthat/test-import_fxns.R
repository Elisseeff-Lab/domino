test_that("create_rl_map_cellphonedb fails on wrong input arg type.", {
  data(CellPhoneDB)
  genes_tiny <- CellPhoneDB$genes_tiny
  proteins_tiny <- CellPhoneDB$proteins_tiny
  interactions_tiny <- CellPhoneDB$interactions_tiny
  complexes_tiny <- CellPhoneDB$complexes_tiny

  expect_error(create_rl_map_cellphonedb(
    genes = list(), proteins = proteins_tiny,
    interactions = interactions_tiny, complexes = complexes_tiny
  ), "Class of genes must be one of: character,data.frame")

  expect_error(create_rl_map_cellphonedb(
    genes = genes_tiny, proteins = list(),
    interactions = interactions_tiny, complexes = complexes_tiny
  ), "Class of proteins must be one of: character,data.frame")

  expect_error(create_rl_map_cellphonedb(
    genes = genes_tiny, proteins = proteins_tiny,
    interactions = list(), complexes = complexes_tiny
  ), "Class of interactions must be one of: character,data.frame")

  expect_error(create_rl_map_cellphonedb(
    genes = genes_tiny, proteins = proteins_tiny,
    interactions = interactions_tiny, complexes = list()
  ), "Class of complexes must be one of: character,data.frame")

  expect_error(create_rl_map_cellphonedb(
    genes = genes_tiny, proteins = proteins_tiny,
    interactions = interactions_tiny, complexes = complexes_tiny,
    database_name = list()
  ), "Class of database_name must be one of: character")

  expect_error(create_rl_map_cellphonedb(
    genes = genes_tiny, proteins = proteins_tiny,
    interactions = interactions_tiny, complexes = complexes_tiny,
    database_name = c("length", ">1")
  ), "Length of database_name must be one of: 1")
})


test_that("create_domino fails on wrong input arg type.", {
  data(PBMC)
  data(SCENIC)
  data(CellPhoneDB)

  rl_map_tiny <- create_rl_map_cellphonedb(genes = CellPhoneDB$genes_tiny,
    proteins = CellPhoneDB$proteins_tiny,
    interactions = CellPhoneDB$interactions_tiny,
    complexes = CellPhoneDB$complexes_tiny)

  auc_tiny <- SCENIC$auc_tiny

  RNA_count_tiny <- PBMC$RNA_count_tiny
  RNA_zscore_tiny <- PBMC$RNA_zscore_tiny
  clusters_tiny <- PBMC$clusters_tiny

  #bad rl map
  bad_rl_map <- "rl_map"
  expect_error(create_domino(bad_rl_map,
                             counts = RNA_count_tiny,
                             z_scores = RNA_zscore_tiny,
                             clusters = clusters_tiny),
  "Class of rl_map must be one of: data.frame")

  bad_rl_map <- rl_map_tiny
  colnames(bad_rl_map) <- paste(colnames(bad_rl_map), "qq")
  expect_error(create_domino(bad_rl_map,
                             counts = RNA_count_tiny,
                             z_scores = RNA_zscore_tiny,
                             clusters = clusters_tiny),
  "Required variables gene_A, gene_B, type_A, type_B not found in rl_map")

  #bad features
  bad_features <- matrix()
  expect_error(create_domino(rl_map_tiny,
                             bad_features,
                             counts = RNA_count_tiny,
                             z_scores = RNA_zscore_tiny,
                             clusters = clusters_tiny),
  "No rownames found in features")

  #seurat or counts, zscores and clusters
  expect_error(create_domino(rl_map_tiny,
                             auc_tiny),
  "Class of counts must be one of: matrix,data.frame")

  expect_error(create_domino(rl_map_tiny,
                             auc_tiny,
                             counts = RNA_count_tiny,
                             z_scores = RNA_zscore_tiny),
  "Class of clusters must be one of: factor")

  #bad rec_min threshold
  expect_error(create_domino(rl_map_tiny,
                             auc_tiny,
                             counts = RNA_count_tiny,
                             z_scores = RNA_zscore_tiny,
                             clusters = clusters_tiny,
                             rec_min_thresh = 20
                             ),
  "All values in rec_min_thresh must be between 0 and 1")

  expect_error(create_domino(rl_map_tiny,
                             auc_tiny,
                             counts = RNA_count_tiny,
                             z_scores = RNA_zscore_tiny,
                             clusters = clusters_tiny,
                             rec_min_thresh = -20
                             ),
  "All values in rec_min_thresh must be between 0 and 1")

  expect_error(create_domino(rl_map_tiny,
                             auc_tiny,
                             counts = RNA_count_tiny,
                             z_scores = RNA_zscore_tiny,
                             clusters = clusters_tiny,
                             tf_selection_method = "non-existent"
                             ),
  "All values in tf_selection_method must be one of: clusters, variable, all")
})