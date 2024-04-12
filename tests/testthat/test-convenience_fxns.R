test_that("rename_clusters function works correctly", {
  # Test domino object:
  dom <- dominoSignal:::pbmc_dom_built_tiny
  
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
