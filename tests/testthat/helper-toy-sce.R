# Shared fixtures for scatools unit tests.

# A tiny SingleCellExperiment with 4 bins on chr1 and 3 cells.
# The last bin is deliberately half-width to exercise length normalization.
make_toy_sce <- function() {
  m <- matrix(
    c(
      1, 2, 3, 4, # cell1
      2, 4, 6, 8, # cell2
      10, 10, 10, 10 # cell3
    ),
    nrow = 4, ncol = 3,
    dimnames = list(paste0("bin", 1:4), paste0("cell", 1:3))
  )
  gr <- GenomicRanges::GRanges(
    seqnames = "chr1",
    ranges = IRanges::IRanges(
      start = c(1, 101, 201, 301),
      end = c(100, 200, 300, 350) # widths 100, 100, 100, 50
    )
  )
  sce <- SingleCellExperiment::SingleCellExperiment(
    assays = list(counts = m),
    rowRanges = gr
  )
  rownames(sce) <- rownames(m)
  colnames(sce) <- colnames(m)
  sce
}
