test_that("rename_clusters function works correctly", {
  # Test domino object:
  data(CellPhoneDB)
  rl_map_tiny <- create_rl_map_cellphonedb(
    genes = CellPhoneDB$genes_tiny,
    proteins = CellPhoneDB$proteins_tiny,
    interactions = CellPhoneDB$interactions_tiny,
    complexes = CellPhoneDB$complexes_tiny)

  data(SCENIC)
  regulon_list_tiny <- create_regulon_list_scenic(
    regulons = SCENIC$regulons_tiny)
  
  data(PBMC)
  pbmc_dom_tiny <- create_domino(
    rl_map = rl_map_tiny, features = SCENIC$auc_tiny,
    counts = PBMC$RNA_count_tiny, z_scores = PBMC$RNA_zscore_tiny,
    clusters = PBMC$clusters_tiny, tf_targets = regulon_list_tiny,
    use_clusters = TRUE, use_complexes = TRUE, remove_rec_dropout = FALSE)

  dom <- build_domino(
    dom = pbmc_dom_tiny, min_tf_pval = .05, max_tf_per_clust = Inf,
    max_rec_per_tf = Inf, rec_tf_cor_threshold = .1, min_rec_percentage = 0.01
  )

  # Define the cluster conversion, Z/new_clust does not match data intentionally
  clust_conv <- c("W", "X", "Y", "Z")
  names(clust_conv) <- c(levels(dom@clusters[1:3]), "new_clust")

  # Run the function
  dom_renamed <- rename_clusters(dom, clust_conv)

  # Check that the clusters were renamed correctly
  expect_equal(levels(dom_renamed@clusters), c("W", "X", "Y"))
  expect_equal(colnames(dom_renamed@clust_de), c("W", "X", "Y"))
  expect_equal(names(dom_renamed@linkages$clust_tf), c("W", "X", "Y"))
  expect_equal(names(dom_renamed@linkages$clust_rec), c("X", "Y")) # only CD8 T Cell and CD14 monocyte present
  expect_equal(names(dom_renamed@linkages$clust_incoming_lig), c("X", "Y")) # only CD14 monocyte and CD8 T Cell present
  expect_equal(names(dom_renamed@linkages$clust_tf_rec), c("W", "X", "Y"))
  expect_equal(colnames(dom_renamed@signaling), paste0("L_", c("W", "X", "Y")))
  expect_equal(rownames(dom_renamed@signaling), paste0("R_", c("W", "X", "Y")))
  expect_equal(names(dom_renamed@cl_signaling_matrices), c("W", "X", "Y"))
  expect_equal(colnames(dom_renamed@cl_signaling_matrices[["W"]]),
               paste0("L_", c("W", "X", "Y")))
  expect_equal(colnames(dom_renamed@cl_signaling_matrices[["X"]]),
               paste0("L_", c("W", "X", "Y")))
  expect_equal(colnames(dom_renamed@cl_signaling_matrices[["Y"]]),
               paste0("L_", c("W", "X", "Y")))
})
