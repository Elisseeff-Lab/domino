## unit tests use internal data stored in R/create_domino.R

test_that(
  "write_rec_lig_linkages: generate list of receptors capable of interacting with receptors", {
    rl_parse_tiny <- data.frame(
      "R.gene" = c("GENEB1,GENEB2", "GENEC", "GENEC", "GENEB1,GENEB2", "GENEE1,GENEE2"),
      "L.gene" = c("GENEA", "GENEA", "GENED1,GENED2", "GENED1,GENED2", "GENEF"),
      "R.name" = c("complexB", "GENEC", "GENEC", "complexB", "complexE"),
      "L.name" = c("GENEA", "GENEA", "complexD", "complexD", "GENEF")
    )
    
    expected_rec_lig <- list(
      "complexB" = c("GENEA", "complexD"),
      "GENEC" = c("GENEA", "complexD"),
      "complexE" = c("GENEF")
    )
    
    expect_equal(
      write_rec_lig_linkages(rl_parse_tiny),
      expected_rec_lig
    )
  }
)
