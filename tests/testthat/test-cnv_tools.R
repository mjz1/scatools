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

test_that("identify_normal warns rather than errors when n_normal_clusts is too large", {
  # Regression test: the warning string interpolated an unbalanced `{}`
  # expression, so cli failed to parse it and this branch raised
  # "Could not parse cli `{}` expression" instead of warning. See #9.
  sce <- make_toy_sce()
  sce$clusters <- c("c1", "c2", "c3")

  # cli_alert_warning signals a message condition, not a warning.
  expect_message(
    res <- identify_normal(
      sce,
      assay_name = "counts",
      group_by = "clusters",
      method = "min_sd",
      n_normal_clusts = 3,
      plot = FALSE
    ),
    "Setting n_normal_clusts to 2"
  )

  # n_normal_clusts is clamped to one fewer than the number of clusters.
  expect_equal(sum(!res$tumor_cell), 2)
  expect_equal(sum(res$tumor_cell), 1)
})

test_that("log_debug is silent unless scatools.debug is enabled", {
  withr::local_options(scatools.debug = NULL)
  expect_silent(log_debug("should not appear"))

  withr::local_options(scatools.debug = TRUE)
  expect_message(log_debug("should appear"), "should appear")
})
