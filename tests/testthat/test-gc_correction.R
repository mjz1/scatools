test_that("gc_cor_modal reduces GC bias and preserves length", {
  skip_if_not_installed("quantreg")
  skip_if_not_installed("polynom")

  set.seed(42)
  n <- 300
  gc <- runif(n, 0.35, 0.65)
  # Flat underlying copy number with a strong multiplicative GC bias.
  counts <- 100 * (1 + 1.5 * (gc - 0.5)) * exp(rnorm(n, 0, 0.05))
  counts <- pmax(counts, 1)

  corrected <- gc_cor_modal(counts = counts, gc = gc, results = "counts")

  expect_length(corrected, n)

  ok <- is.finite(corrected)
  expect_gt(sum(ok), n * 0.8) # most bins should be corrected

  raw_cor <- abs(stats::cor(counts, gc))
  corrected_cor <- abs(stats::cor(corrected[ok], gc[ok]))

  # Correction should substantially weaken the GC-count relationship.
  expect_lt(corrected_cor, raw_cor)
})

test_that("gc_cor_modal errors on mismatched input lengths", {
  expect_error(
    gc_cor_modal(counts = 1:10, gc = 1:5),
    "identical"
  )
})
