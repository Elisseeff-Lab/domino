# generate domino objects for comparison in testing scripts

library(domino2, lib = "../domino_libraries/libraries/v0_2_1/")

test_files <- c(
  "data/RNA_count_test.rda",
  "data/RNA_zscore_test.rda",
  "data/cluster_test.rda",
  "data/auc_test.rda",
  "data/regulons_test.rda",
  "data/complexes_test.rda",
  "data/genes_test.rda",
  "data/interactions_test.rda",
  "data/proteins_test.rda"
)

for(file in test_files) load(file)

rl_map_test <- domino2::create_rl_map_cellphonedb(
  genes = genes_test,
  proteins = proteins_test,
  interactions = interactions_test,
  complexes = complexes_test
)
save(rl_map_test, file = "data/rl_map_test.rda")

regulon_list_test <- create_regulon_list_scenic(
  regulons = regulons_test
)
save(regulon_list_test, file = "data/regulon_list_test.rda")

pbmc_dom_test <- create_domino(
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
save(pbmc_dom_test, file = "data/pbmc_dom_test.rda")

pbmc_dom_built_test <- build_domino(
  dom = pbmc_dom_test,
  min_tf_pval = .05,
  max_tf_per_clust = Inf,
  max_rec_per_tf = Inf,
  rec_tf_cor_threshold = .15,
  min_rec_percentage = 0.03 
)
save(pbmc_dom_built_test, file = "data/pbmc_dom_built_test.rda")


# usethis::use_data(pbmc_dom, pbmc_dom_b, internal = TRUE)

