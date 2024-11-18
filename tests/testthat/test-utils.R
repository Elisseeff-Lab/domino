test_that("bool conversion function works",{
  df <- data.frame(list(c1 = c("True", "False"),
                        c2 = c("False", "True"),
                        c3 = c(1, 2)),
                        c4 = c("a", "b"))

  c_df <- conv_py_bools(df)

  expect_equal(class(c_df$c1), "logical")
  expect_equal(class(c_df$c2), "logical")
  expect_equal(class(c_df$c3), "numeric")
  expect_equal(class(c_df$c4), "character")

})

test_that("read if char tries to read a file", {
  expect_error(read_if_char("./file_that_not_exists.csv",
                            "cannot open the connection"))
  expect_error(read_if_char(c('a', 'b')), "Length of obj must be one of: 1")
})

test_that("mandatory field absence yields error, presence does not", {
  expect_error(check_arg(arg = data.frame(a = c(1, 2), b = c(3, 4)),
                         allow_class = "data.frame",
                         need_vars = c("c")))
  expect_silent(check_arg(arg = data.frame(a = c(1, 2), b = c(3, 4)),
                          allow_class = "data.frame",
                          need_vars = c("a", "b")))
})

test_that("range checker works", {
  expect_error(check_arg(1, allow_class = "numeric", allow_range = c(2, 5)),
               "All values in 1 must be between 2 and 5")
  expect_silent(check_arg(3, allow_class = "numeric", allow_range = c(2, 5)))
})

test_that("check_arg works for class", {
  expect_error(check_arg(1, allow_class = "data.frame"),
               "Class of 1 must be one of: data.frame")
  expect_silent(check_arg(data.frame(a = c(1, 2), b = c(3, 4)),
                          allow_class = "data.frame"))
  expect_silent(check_arg(data.frame(a = c(1, 2), b = c(3, 4)),
                          allow_class = c("data.frame", "matrix")))
})

test_that("resolve_complexes maps complex names to their component genes, resolve_names maps ligand name to gene name", {
  data(CellPhoneDB)
  rl_map_tiny <- create_rl_map_cellphonedb(
    genes = CellPhoneDB$genes_tiny,
    proteins = CellPhoneDB$proteins_tiny,
    interactions = CellPhoneDB$interactions_tiny,
    complexes = CellPhoneDB$complexes_tiny)
  entry_for_resolve_names <- c(int_pair = "12oxoLeukotrieneB4_byPTGR1 & LTB4R", 
                               name_A = "12oxoLeukotrieneB4_byPTGR1", uniprot_A = "Q14914", gene_A = "PTGR1", type_A = "L", 
                               name_B = "LTB4R", uniprot_B = "Q15722", gene_B = "LTB4R", type_B = "R", 
                               annotation_strategy = "curated", 
                               source = "HMRbase;uniprot;reactome", 
                               database_name = "CellPhoneDB")
  rl_map_tiny <- rbind(rl_map_tiny, entry_for_resolve_names)
  
  data(SCENIC)
  regulon_list_tiny <- create_regulon_list_scenic(
    regulons = SCENIC$regulons_tiny)
  
  data(PBMC)
  pbmc_dom_tiny <- create_domino(
    rl_map = rl_map_tiny, features = SCENIC$auc_tiny,
    counts = PBMC$RNA_count_tiny, z_scores = PBMC$RNA_zscore_tiny,
    clusters = PBMC$clusters_tiny, tf_targets = regulon_list_tiny,
    use_clusters = TRUE, use_complexes = TRUE, remove_rec_dropout = FALSE)
  
  # Test resolve_names
  genes <- c("integrin_a6b4_complex", "TGFB3", "12oxoLeukotrieneB4_byPTGR1")
  ligand_names_resolved <- c("integrin_a6b4_complex" = "integrin_a6b4_complex", 
                             "TGFB3" = "TGFB3", 
                             "12oxoLeukotrieneB4_byPTGR1" = "PTGR1")
  expect_equal(resolve_names(pbmc_dom_tiny, genes), ligand_names_resolved)
  
  # Test resolve_complexes
  genes <- c("integrin_a6b4_complex", "IL7_receptor", "TGFBR3")
  complexes_resolved <- list("integrin_a6b4_complex" = c("ITGB4", "ITGA6"),
                             "IL7_receptor" = c("IL7R", "IL2RG"),
                             "TGFBR3" = "TGFBR3")
  expect_equal(resolve_complexes(pbmc_dom_tiny, genes), complexes_resolved)
})
