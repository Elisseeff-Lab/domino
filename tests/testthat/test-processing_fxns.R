## unit tests use internal data stored in R/sysdata.R

test_that("build_domino runs correctly", {
  # use rl_map created in domino2 v0.2.1
  data(pbmc_dom_built_tiny, pbmc_dom_tiny)
  dom_build <- build_domino(
    pbmc_dom_tiny,
    min_tf_pval = .05,
    max_tf_per_clust = Inf,
    max_rec_per_tf = Inf,
    rec_tf_cor_threshold = .1,
    min_rec_percentage = 0.01
  )
  expect_equal(dom_build, pbmc_dom_built_tiny)
})
