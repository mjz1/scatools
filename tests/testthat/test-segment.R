test_that("segment_cnv tolerates NA bins without crashing", {
  np <- 25
  gr <- GenomicRanges::GRanges(
    rep(c("chr1", "chr2"), each = np),
    IRanges::IRanges(start = rep(seq(1, by = 1e6, length.out = np), 2), width = 1e6)
  )
  m <- matrix(1,
    nrow = 2 * np, ncol = 3,
    dimnames = list(paste0("b", seq_len(2 * np)), paste0("c", 1:3))
  )
  chr1_idx <- seq_len(np)
  chr2_idx <- (np + 1):(2 * np)
  m[chr2_idx, 2] <- 2 # cell c2: gain on chr2
  m[chr2_idx, 3] <- NA # cell c3: an entire chromosome is NA (the crashing pattern)

  sce <- SingleCellExperiment::SingleCellExperiment(
    assays = list(ratios = m), rowRanges = gr
  )
  rownames(sce) <- rownames(m)

  expect_no_error(
    sce <- segment_cnv(sce, assay_name = "ratios", arm_aware = FALSE, bpparam = BiocParallel::SerialParam())
  )
  seg <- SummarizedExperiment::assay(sce, "ratios_segment")

  # c3: the all-NA chromosome stays NA; the valid chromosome is segmented.
  expect_true(all(is.na(seg[chr2_idx, "c3"])))
  expect_true(all(!is.na(seg[chr1_idx, "c3"])))
  # c2's gain is recovered.
  expect_gt(mean(seg[chr2_idx, "c2"], na.rm = TRUE), mean(seg[chr1_idx, "c2"], na.rm = TRUE))
})

test_that("segment_cnv is arm-aware when arm annotation is present", {
  np <- 15
  # One chromosome, two arms (p, q); a step change exactly at the arm boundary.
  gr <- GenomicRanges::GRanges(
    "chr1",
    IRanges::IRanges(start = seq(1, by = 1e6, length.out = 2 * np), width = 1e6)
  )
  gr$arm <- rep(c("p", "q"), each = np)
  m <- matrix(c(rep(1, np), rep(3, np)),
    nrow = 2 * np, ncol = 1,
    dimnames = list(paste0("b", seq_len(2 * np)), "c1")
  )
  sce <- SingleCellExperiment::SingleCellExperiment(assays = list(ratios = m), rowRanges = gr)
  rownames(sce) <- rownames(m)

  sce <- segment_cnv(sce, assay_name = "ratios", arm_aware = TRUE, bpparam = BiocParallel::SerialParam())
  seg <- SummarizedExperiment::assay(sce, "ratios_segment")[, 1]
  # p arm segments to ~1, q arm to ~3 (the boundary is respected)
  expect_lt(mean(seg[1:np]), mean(seg[(np + 1):(2 * np)]))
  expect_equal(round(mean(seg[1:np]), 1), 1)
  expect_equal(round(mean(seg[(np + 1):(2 * np)]), 1), 3)
})
