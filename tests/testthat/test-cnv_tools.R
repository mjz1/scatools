test_that("calc_ratios divides each cell by its own centre", {
  sce <- make_toy_sce()
  res <- calc_ratios(sce, assay_name = "counts", fun = "mean")
  ratios <- SummarizedExperiment::assay(res, "counts_ratios")

  # cell1 counts c(1,2,3,4), mean 2.5 -> c(0.4,0.8,1.2,1.6)
  expect_equal(unname(ratios[, "cell1"]), c(0.4, 0.8, 1.2, 1.6))
  # cell3 is flat -> all ratios 1
  expect_equal(unname(ratios[, "cell3"]), c(1, 1, 1, 1))
})

test_that("logNorm applies the requested log transform", {
  sce <- make_toy_sce()
  SummarizedExperiment::assay(sce, "myratios") <- SummarizedExperiment::assay(sce, "counts")

  res <- logNorm(sce, transform = "log2", assay_name = "myratios", name = "mylog")
  out <- SummarizedExperiment::assay(res, "mylog")

  expect_equal(unname(out[, "cell1"]), round(log2(c(1, 2, 3, 4)), 2))
  expect_equal(unname(out[, "cell3"]), round(log2(c(10, 10, 10, 10)), 2))
})

test_that("logNorm replaces zeros before logging (no -Inf)", {
  sce <- make_toy_sce()
  m <- SummarizedExperiment::assay(sce, "counts")
  m[1, 1] <- 0
  SummarizedExperiment::assay(sce, "myratios") <- m

  res <- logNorm(sce, transform = "log2", assay_name = "myratios", name = "mylog")
  out <- SummarizedExperiment::assay(res, "mylog")

  expect_true(all(is.finite(out)))
  # 0 -> 1e-3 -> log2(1e-3)
  expect_equal(out[1, 1], round(log2(1e-3), 2))
})
