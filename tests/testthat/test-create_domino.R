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
