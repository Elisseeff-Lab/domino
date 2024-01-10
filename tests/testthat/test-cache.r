test_that("can add to a new cache", {
  bfc <- BiocFileCache::BiocFileCache(tempdir(), ask = FALSE)
  expect_equal(length(bfc), 0)

  #from Bioc example
  savepath <- BiocFileCache::bfcnew(bfc, "NewResource", ext = ".RData")
  df <- data.frame(list("a"=c(1, 2), "b"=c(3, 4)))
  save(df, file = savepath)

  expect_equal(length(bfc), 1)

  #cleanup
  BiocFileCache::removebfc(bfc, ask = FALSE)
})
