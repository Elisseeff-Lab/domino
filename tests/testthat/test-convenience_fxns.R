test_that("rename_clusters function works correctly", {
  # Create a mock domino object
  dom <- domino()
  dom@clusters <- factor(c("A", "B"))
  mock_cl <- matrix(1:4, nrow = 2,
                    dimnames = list(c("A", "B"), c("gene1", "gene2")))
  dom@misc <- list(create = TRUE, build = TRUE, cl_rec_percent = mock_cl)
  dom@clust_de <- mock_cl
  dom@linkages <- list(
    clust_tf = c("A", "B"),
    clust_rec = c("A", "B"),
    clust_incoming_lig = c("A", "B"),
    clust_tf_rec = c("A", "B")
  )

  mock_sm <- matrix(1:4, nrow = 2,
                    dimnames = list(paste0("R_", c("A", "B")),
                                    paste0("L_", c("A", "B"))))
  dom@signaling <- mock_sm
  dom@cl_signaling_matrices <- list(A = mock_sm, B = mock_sm)

  # Define the cluster conversion, C=Z does not match data intentionally
  clust_conv <- c("A" = "X", "B" = "Y", "C" = "Z")

  # Run the function
  dom_renamed <- rename_clusters(dom, clust_conv)

  # Check that the clusters were renamed correctly
  expect_equal(dom_renamed@clusters, factor(c("X", "Y")))
  expect_equal(colnames(dom_renamed@clust_de), c("X", "Y"))
  expect_equal(names(dom_renamed@linkages$clust_tf), c("X", "Y"))
  expect_equal(names(dom_renamed@linkages$clust_rec), c("X", "Y"))
  expect_equal(names(dom_renamed@linkages$clust_incoming_lig), c("X", "Y"))
  expect_equal(names(dom_renamed@linkages$clust_tf_rec), c("X", "Y"))
  expect_equal(colnames(dom_renamed@signaling), paste0("L_", c("X", "Y")))
  expect_equal(rownames(dom_renamed@signaling), paste0("R_", c("X", "Y")))
  expect_equal(names(dom_renamed@cl_signaling_matrices), c("X", "Y"))
  expect_equal(colnames(dom_renamed@cl_signaling_matrices[["X"]]),
               paste0("L_", c("X", "Y")))
  expect_equal(colnames(dom_renamed@cl_signaling_matrices[["Y"]]),
               paste0("L_", c("X", "Y")))
})