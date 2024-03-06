# unit tests for read_rl_map

test_that(
  "read_rl_map_genes: ", {
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
