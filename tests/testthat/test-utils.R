test_that("getmode returns the most frequent value", {
  expect_equal(getmode(c(1, 1, 2)), 1)
  expect_equal(getmode(c(5, 5, 5, 3, 3)), 5)
  expect_equal(getmode(c("a", "b", "b")), "b")
})

test_that("prettyMb formats genomic distances with the right unit", {
  expect_equal(prettyMb(1e7), "10Mb")
  expect_equal(prettyMb(5e5), "500Kb")
  expect_equal(prettyMb(1e3), "1Kb")
  expect_equal(prettyMb(2.5e8), "250Mb")
})
