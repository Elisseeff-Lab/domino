## regression tests use internal data stored in R/sysdata.R
## see example of creation in d8b0e11

## large scale test that processing functions generate results matching the state of the code in v0.2.1
test_that("creation of rl_map from CellPhoneDB v4 input", {
  # use rl_map created in domino2 v0.2.1
  rl_map <- create_rl_map_cellphonedb(
    genes = v0.2.1$genes_tiny,
    proteins = v0.2.1$proteins_tiny,
    interactions = v0.2.1$interactions_tiny,
    complexes = v0.2.1$complexes_tiny
  )
  expect_equal(rl_map, v0.2.1$rl_map_tiny)
})

test_that("formating of SCENIC regulons output as a list", {
  # use regulon list created in domino2 v0.2.1
  regulon_ls <- create_regulon_list_scenic(regulons = v0.2.1$regulons_tiny)
  expect_equal(regulon_ls, v0.2.1$regulon_list_tiny)
})

test_that("creation of a domino object from SCENIC and CellPhoneDB inputs", {
  # domino object created in domino2 v0.2.1
  # rl_map created in domino2 v0.2.1 from CellPhoneDB v4
  # scenic inputs from pySCENIC v0.11.0
  # minimal expression data from 360 cells in 3 cell types,
  # expression features for 16 genes
  # 360 cell barcodes have annotations of cell type assignment in clusters_tiny
  pbmc_dom <- create_domino(
    rl_map = v0.2.1$rl_map_tiny,
    features = v0.2.1$auc_tiny,
    counts = v0.2.1$RNA_count_tiny,
    z_scores = v0.2.1$RNA_zscore_tiny,
    clusters = v0.2.1$clusters_tiny,
    tf_targets = v0.2.1$regulon_list_tiny,
    use_clusters = TRUE,
    use_complexes = TRUE,
    remove_rec_dropout = FALSE
  )
  expect_equal(pbmc_dom@linkages$complexes, v0.2.1$pbmc_dom_tiny@linkages$complexes)
  expect_equal(pbmc_dom@linkages$rec_lig, v0.2.1$pbmc_dom_tiny@linkages$rec_lig)
  expect_equal(pbmc_dom@linkages$tf_targets, v0.2.1$pbmc_dom_tiny@linkages$tf_targets)
  expect_equal(pbmc_dom@clust_de, v0.2.1$pbmc_dom_tiny@clust_de)
  expect_equal(pbmc_dom@cor, v0.2.1$pbmc_dom_tiny@cor)
})

test_that("building a domino object under set parameters", {
  # built domino object created in domino2 v0.2.1
  # domino object created in domino2 v0.2.1
  pbmc_dom_built <- build_domino(
    dom = v0.2.1$pbmc_dom_tiny,
    min_tf_pval = .05,
    max_tf_per_clust = Inf,
    max_rec_per_tf = Inf,
    rec_tf_cor_threshold = .1,
    min_rec_percentage = 0.01
  )
  expect_equal(pbmc_dom_built@linkages$clust_tf,
               v0.2.1$pbmc_dom_built_tiny@linkages$clust_tf)

  expect_equal(pbmc_dom_built@linkages$clust_rec,
               v0.2.1$pbmc_dom_built_tiny@linkages$clust_rec)

  expect_equal(pbmc_dom_built@linkages$clust_tf_rec,
               v0.2.1$pbmc_dom_built_tiny@linkages$clust_tf_rec)

  expect_equal(pbmc_dom_built@linkages$clust_incoming_lig,
               v0.2.1$pbmc_dom_built_tiny@linkages$clust_incoming_lig)

  expect_equal(pbmc_dom_built@linkages$tf_rec,
               v0.2.1$pbmc_dom_built_tiny@linkages$tf_rec)

  expect_equal(pbmc_dom_built@cl_signaling_matrices,
               v0.2.1$pbmc_dom_built_tiny@cl_signaling_matrices)

  expect_equal(pbmc_dom_built@signaling,
               v0.2.1$pbmc_dom_built_tiny@signaling)
})