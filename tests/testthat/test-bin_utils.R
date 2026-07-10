test_that("is_valid_bin flags bins on read count and N frequency", {
  valid <- is_valid_bin(
    counts = c(0, 5, 10),
    n_freq = c(0, 0, 0.1),
    min_reads = 1,
    max_N_freq = 0.05
  )
  expect_equal(valid, c(FALSE, TRUE, FALSE))
})

test_that("is_ideal_bin returns valid/ideal invariants", {
  counts <- c(0, 5, 100, 5, 6)
  gc <- c(0.4, 0.5, 0.5, 0.5, 0.5)
  n_freq <- c(0, 0, 0, 0.1, 0)

  res <- is_ideal_bin(counts, gc, n_freq, min_reads = 1, max_N_freq = 0.05)

  expect_s3_class(res, "data.frame")
  expect_named(res, c("ideal", "valid"))
  expect_equal(nrow(res), length(counts))

  # valid column must agree with is_valid_bin
  expect_equal(res$valid, is_valid_bin(counts, n_freq, min_reads = 1, max_N_freq = 0.05))

  # ideal bins are always a subset of valid bins
  expect_true(all(res$valid[res$ideal]))

  # a zero-count bin can never be valid or ideal
  expect_false(res$valid[1])
  expect_false(res$ideal[1])
})

test_that("bin id encode/decode round-trips", {
  gr <- GenomicRanges::GRanges(
    seqnames = c("chr1", "chr2"),
    ranges = IRanges::IRanges(start = c(1, 201), end = c(100, 300))
  )
  ids <- get_bin_ids(gr)
  expect_equal(ids, c("chr1_1_100", "chr2_201_300"))

  info <- get_bin_info(ids)
  expect_equal(as.character(info$chr), c("chr1", "chr2"))
  expect_equal(info$start, c(1, 201))
  expect_equal(info$end, c(100, 300))
})

test_that("length_normalize rescales only variable-width bins", {
  sce <- make_toy_sce()
  res <- length_normalize(sce, assay_name = "counts", assay_to = "counts_lenNorm")

  norm <- SummarizedExperiment::assay(res, "counts_lenNorm")
  orig <- SummarizedExperiment::assay(sce, "counts")

  # bins 1-3 are the modal width (100) -> unchanged
  expect_equal(norm[1:3, ], orig[1:3, ])
  # bin 4 is half width (50) -> counts doubled
  expect_equal(norm[4, ], orig[4, ] * 2)
})
