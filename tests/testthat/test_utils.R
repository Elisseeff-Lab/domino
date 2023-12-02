context("domino2")

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