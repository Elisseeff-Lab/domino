context("domino2")

## large scale test that processing functions generate results matching the state of the code in v0.2.1
test_that("creation of rl_map from CellPhoneDB v4 input", {
  # # load the rl_map created in domino2 v0.2.1
  load(test_path("testdata", "sysdata.rda"))
  rl_map <- create_rl_map_cellphonedb(
    genes = genes_test,
    proteins = proteins_test,
    interactions = interactions_test,
    complexes = complexes_test
  )
  expect_equal(rl_map, rl_map_test)
})

test_that("formating of SCENIC regulons output as a list", {
  # # load the regulon list created in domino2 v0.2.1
  load(test_path("testdata", "sysdata.rda"))
  regulon_ls <- create_regulon_list_scenic(regulons = regulons_test)
  expect_equal(regulon_ls, regulon_list_test)
})

test_that("creation of a domino object from SCENIC and CellPhoneDB inputs", {
  # # load the domino object created in domino2 v0.2.1
  # # load the rl_map created in domino2 v0.2.1 from CellPhoneDB v4
  # # load Scenic inputs to create a domino result
  # # load scRNAseq expression data
  # # named vector of cell type cluster annotations for cell barcodes
  load(test_path("testdata", "sysdata.rda"))
  pbmc_dom <- create_domino(
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
  expect_equal(pbmc_dom@linkages$complexes, pbmc_dom_test@linkages$complexes)
  expect_equal(pbmc_dom@linkages$rec_lig, pbmc_dom_test@linkages$rec_lig)
  expect_equal(pbmc_dom@linkages$tf_targets, pbmc_dom_test@linkages$tf_targets)
  expect_equal(pbmc_dom@clust_de, pbmc_dom_test@clust_de)
  expect_equal(pbmc_dom@cor, pbmc_dom_test@cor)
})

test_that("building a domino object under set parameters", {
  # # load build domino object created in domino2 v0.2.1
  # # load domino object created in domino2 v0.2.1
  load(test_path("testdata", "sysdata.rda"))
  pbmc_dom_built <- build_domino(
    dom = pbmc_dom_test,
    min_tf_pval = .05,
    max_tf_per_clust = Inf,
    max_rec_per_tf = Inf,
    rec_tf_cor_threshold = .15,
    min_rec_percentage = 0.03
  )
  expect_equal(pbmc_dom_built@linkages$clust_tf,
               pbmc_dom_built_test@linkages$clust_tf)

  expect_equal(pbmc_dom_built@linkages$clust_rec,
               pbmc_dom_built_test@linkages$clust_rec)

  expect_equal(pbmc_dom_built@linkages$clust_tf_rec,
               pbmc_dom_built_test@linkages$clust_tf_rec)

  expect_equal(pbmc_dom_built@linkages$clust_incoming_lig,
               pbmc_dom_built_test@linkages$clust_incoming_lig)

  expect_equal(pbmc_dom_built@linkages$tf_rec,
               pbmc_dom_built_test@linkages$tf_rec)

  expect_equal(pbmc_dom_built@cl_signaling_matrices,
               pbmc_dom_built_test@cl_signaling_matrices)

  expect_equal(pbmc_dom_built@signaling,
               pbmc_dom_built_test@signaling)
})
