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
})

test_that("mandatory field absence yields error, presence does not", {
  expect_error(check_arg(arg = data.frame(a = c(1, 2), b = c(3, 4)),
                         allow_class = "data.frame",
                         require_vars = c("c")))
  expect_silent(check_arg(arg = data.frame(a = c(1, 2), b = c(3, 4)),
                          allow_class = "data.frame",
                          require_vars = c("a", "b")))
})
