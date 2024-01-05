test_that("can add to a new cache", {
  tf <- tempfile(tmpdir = tempdir(check = TRUE))
  file.create(tf)
  bfc <- cache(path = tempdir())
  expect_equal(length(bfc), 0)
  BiocFileCache::cleanbfc(bfc, ask = FALSE)
  BiocFileCache::removebfc(bfc, ask = FALSE)
})
