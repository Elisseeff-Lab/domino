context("domino2")

test_that("standard dataset on v0.2.1",{
  expect_equal(1, 1)
})

## large scale test that processing functions generate results matching the state of the code in v0.2.1
test_that("creation of rl_map from CellPhoneDB v4 input", {
  # # load the rl_map created in domino2 v0.2.1
  # load(test_path("../../data/rl_map_test.rda"))
  # # construct rl_map under current package
  # load(test_path("../../data/complexes_test.rda"))
  # load(test_path("../../data/genes_test.rda"))
  # load(test_path("../../data/interactions_test.rda"))
  # load(test_path("../../data/proteins_test.rda"))
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
  # load(test_path("../../data/regulon_list_test.rda"))
  # # load SCENIC regulon output
  # load(test_path("../../data/regulons_test.rda"))
  load(test_path("testdata", "sysdata.rda"))
  regulon_ls <- create_regulon_list_scenic(regulons = regulons_test)
  expect_equal(regulon_ls, regulon_list_test)
})
test_that("creation of a domino object from SCENIC and CellPhoneDB inputs", {
  # # load the domino object created in domino2 v0.2.1
  # load(test_path("../../data/pbmc_dom_test.rda"))
  # # load the rl_map created in domino2 v0.2.1 from CellPhoneDB v4
  # load(test_path("../../data/rl_map_test.rda"))
  # # load Scenic inputs to create a domino result
  # load(test_path("../../data/auc_test.rda"))
  # load(test_path("../../data/regulon_list_test.rda"))
  # # load scRNAseq expression data
  # load(test_path("../../data/RNA_count_test.rda"))
  # load(test_path("../../data/RNA_zscore_test.rda"))
  # # named vector of cell type cluster annotations for cell barcodes
  # load(test_path("../../data/cluster_test.rda"))
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
  expect_equal(pbmc_dom, pbmc_dom_test)
})
test_that("building a domino object under set parameters", {
  # # load build domino object created in domino2 v0.2.1
  # load(test_path("../../data/pbmc_dom_built_test.rda"))
  # # load domino object created in domino2 v0.2.1
  # load(test_path("../../data/pbmc_dom_test.rda"))
  load(test_path("testdata", "sysdata.rda"))
  pbmc_dom_built <- build_domino(
    dom = pbmc_dom_test,
    min_tf_pval = .05,
    max_tf_per_clust = Inf,
    max_rec_per_tf = Inf,
    rec_tf_cor_threshold = .15,
    min_rec_percentage = 0.03 
  )
  expect_equal(pbmc_dom_built, pbmc_dom_built_test)
})

