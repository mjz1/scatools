test_that("do_qc computes bin/cell metrics and a summary", {
  sce <- make_toy_sce()
  res <- do_qc(sce, assay_name = "counts", plot = FALSE)

  expect_true(all(c("counts_per_cell", "median_percell_counts") %in% colnames(SummarizedExperiment::colData(res))))
  expect_true(all(c("counts_per_bin", "mean_perbin_counts", "frac_zero_perbin") %in%
    colnames(SummarizedExperiment::rowData(res))))

  qc <- S4Vectors::metadata(res)$qc_summary
  expect_equal(qc$n_cells, ncol(sce))
  expect_equal(qc$n_bins, nrow(sce))
  expect_equal(qc$assay, "counts")
  # frac_zero must be in [0, 1]
  expect_true(all(SummarizedExperiment::rowData(res)$frac_zero_perbin >= 0 &
    SummarizedExperiment::rowData(res)$frac_zero_perbin <= 1))
})

test_that("do_qc respects assay_name (not hardcoded to counts)", {
  sce <- make_toy_sce()
  SummarizedExperiment::assay(sce, "scaled") <- SummarizedExperiment::assay(sce, "counts") * 10
  res <- do_qc(sce, assay_name = "scaled", plot = FALSE)
  # counts_per_cell should reflect the 'scaled' assay (10x), not 'counts'
  expect_equal(
    unname(res$counts_per_cell),
    unname(Matrix::colSums(SummarizedExperiment::assay(sce, "scaled")))
  )
})
