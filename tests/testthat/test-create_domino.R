## unit tests use internal data stored in R/create_domino.R

test_that(
  "read_rl_map_complexes: generate list of annotated complexes", {
    rl_parse_tiny <- data.frame(
      "R.gene" = c("GENEB1,GENEB2", "GENEC", "GENEC", "GENEB1,GENEB2"),
      "L.gene" = c("GENEA", "GENEA", "GENED1,GENED2", "GENED1,GENED2"),
      "R.name" = c("complexB", "GENEC", "GENEC", "complexB"),
      "L.name" = c("GENEA", "GENEA", "complexD", "complexD")
    )
    expected_complex_list <- list(
      "complexD" = c("GENED1", "GENED2"),
      "complexB" = c("GENEB1", "GENEB2")
    )
    expect_equal(
      read_rl_map_complexes(rl_parse_tiny, use_complexes = TRUE),
      expected_complex_list
    )
    expect_equal(
      read_rl_map_complexes(rl_parse_tiny, use_complexes = FALSE),
      NULL
    )
  }
)
test_that(
  "read_rl_map_complexes: warning that complexes with the same name have different components", {
    rl_parse_repeated <- data.frame(
      "R.gene" = c("GENEB1,GENEB2", "GENEB1,GENEB3"),
      "L.gene" = c("GENEA", "GENEA"),
      "R.name" = c("complexB", "complexB"),
      "L.name" = c("GENEA", "GENEA")
    )
    expect_warning(read_rl_map_complexes(rl_parse_repeated, use_complexes = TRUE))
  }
)
test_that(
  "read_rl_map_genes: Convert rl_map to columns stored in dom@misc$rl_map", {
    rl_map_tiny <- data.frame(
      gene_A = c("GENEA", "GENEA"),
      protein_A = c("simpleA", "simpleA"),
      type_A = c("L", "L"),
      name_A = c("GENEA", "GENEA"),
      gene_B = c("GENEB1,GENEB2", "GENEC"),
      protein_B = c("complexB1,complexB2", "simpleC"),
      type_B = c("R", "R"),
      name_B = c("complexB", "GENEC"),
      int_pair = c("GENEA & complexB", "GENEA & GENEC"),
      annotation_strategy = c("test", "test"),
      source = c("test", "test"),
      database_name = c("CellPhoneDB", "CellPhoneDB")
    )
    expected_rl_parse <- data.frame(
      "R.gene" = c("GENEB1,GENEB2", "GENEC"),
      "L.gene" = c("GENEA", "GENEA"),
      "R.name" = c("complexB", "GENEC"),
      "L.name" = c("GENEA", "GENEA")
    )
    expect_equal(
      read_rl_map_genes(rl_map_tiny),
      expected_rl_parse
    )
  }
)
