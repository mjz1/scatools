# The mappability bigWig is large (~1.3 GB) and lives outside the repo, so this
# test only runs when SCATOOLS_MAP_BW points at a local track.
test_that("add_map_freq computes per-bin mappability in [0, 1]", {
  bw <- Sys.getenv("SCATOOLS_MAP_BW", unset = "")
  skip_if_not(nzchar(bw) && file.exists(bw), "mappability bigWig not available")
  skip_if_not_installed("rtracklayer")

  data(bins_10mb, package = "scatools")
  bins <- add_map_freq(bins_10mb, bw)

  expect_true("map" %in% colnames(S4Vectors::mcols(bins)))
  m <- bins$map
  expect_true(all(m >= 0 & m <= 1, na.rm = TRUE))
  # most 10 Mb bins are highly mappable
  expect_gt(stats::median(m, na.rm = TRUE), 0.8)
})
