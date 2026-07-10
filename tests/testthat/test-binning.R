test_that("bin_frags_chr handles cut-sites that fall outside all bins", {
  # Two bins on chr1 with a gap (101-200) between them — mimics arm-aware bins
  # that don't cover centromeres. A fragment landing in the gap must not break
  # the sparse-matrix construction.
  bins <- GenomicRanges::GRanges(
    "chr1",
    IRanges::IRanges(start = c(1, 201), end = c(100, 300))
  )
  frags <- GenomicRanges::GRanges(
    "chr1",
    IRanges::IRanges(start = c(10, 150, 250), end = c(50, 160, 260)),
    barcode = c("cellA", "cellB", "cellA")
  )
  cells <- c("cellA", "cellB")

  expect_no_error(mat <- bin_frags_chr(frags, bins = bins, cells = cells))
  expect_equal(dim(mat), c(2L, 2L))

  m <- as.matrix(mat)
  # frag1 (cellA) start+end both in bin1 -> 2; frag3 (cellA) both in bin2 -> 2;
  # frag2 (cellB) lands in the gap -> contributes nothing.
  expect_equal(unname(m["chr1_1_100", "cellA"]), 2)
  expect_equal(unname(m["chr1_201_300", "cellA"]), 2)
  expect_equal(sum(m[, "cellB"]), 0)
})
