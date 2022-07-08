# TODO: Also code to regenerate mappability per bins from the hg38 hoffman bigwigs
#
# http://hgdownload.soe.ucsc.edu/gbdb/hg38/hoffmanMappability/
# https://bismap.hoffmanlab.org/
#

# Need function to define bins in cells as ideal or valid

#' Modal regression GC Correction for single cell data
#'
#' @param counts Vector of counts in bins
#' @param gc Vector of GC content in bins
#' @param bin_ids Bin ids
#' @param lowess_frac The smoother span. See [stats::lowess()] for more details
#' @param q A tuple of quantile bounds to compute quantile regressions over read counts
#' @param g GC quantiles over which to integrate regression curves
#' @param results Format of results to return. One of `full`, `default`, or `counts`.
#'
#' @return A dataframe with the modal quantile and GC corrected counts. Results are returned based on the `results` parameter as follows:
#' \describe{
#'   \item{full}{Full dataframe containing results from each quantile curve}
#'   \item{default}{Condensed dataframe from the selected modal quantile curve}
#'   \item{counts}{Return only a vector of corrected counts}
#' }
#'
#' @export
#'
gc_cor_modal <- function(counts, gc, bin_ids = seq_along(counts), lowess_frac = 0.2, q = c(0.1, 0.9), g = c(0.1, 0.9), results = c("default", "counts", "full")) {
  if (length(counts) != length(gc)) {
    stop("Length of counts and gc vectors are not identical")
  }

  if (length(q) != 2 | (q[1] > q[2])) {
    stop("Must provide lower and upper quantile")
  }

  results <- match.arg(results)

  # Into a dataframe pre-filtering
  df_regression_all <- data.frame(reads = counts, gc = gc, bin_ids = bin_ids)

  # Filter the data for non-zero gc and read counts
  df_regression <- df_regression_all[df_regression_all$gc > 0, ]

  df_regression <- df_regression_all[df_regression_all$reads > 0, ]

  # 2nd order polynomial quantile regression
  quantiles <- seq(q[1], q[2], by = 0.01)
  quantile_names <- paste0("q", quantiles * 100)

  if (nrow(df_regression) < 10) {
    message("Not enough data points")
    return(df_regression)
  }

  # Fit second order polynormial quantile regression
  poly2_quantile_model <- quantreg::rq(reads ~ gc + I(gc^2), tau = quantiles, data = df_regression)

  # Fit to our data
  poly_quantile_predict <- predict(object = poly2_quantile_model, newdata = df_regression)

  colnames(poly_quantile_predict) <- quantile_names # rename columns

  # Bind to our data
  df_regression <- cbind.data.frame(df_regression, poly_quantile_predict)

  poly2_quantile_params <- as.data.frame(poly2_quantile_model$coefficients)
  colnames(poly2_quantile_params) <- quantile_names # rename columns

  # integration and mode selection
  gc_min <- quantile(df_regression$gc, g[1])
  gc_max <- quantile(df_regression$gc, g[2])

  poly2_quantile_integration <- c(0, apply(X = poly2_quantile_params, MARGIN = 2, FUN = function(params) {
    poly2 <- polynom::polynomial(params)
    integ <- polynom::integral(poly2)
    integrand <- predict(integ, gc_max) - predict(integ, gc_min)
  }))

  # find the modal quantile
  distances <- diff(poly2_quantile_integration)

  df_dist <- data.frame(quantiles = quantiles, quantile_names = quantile_names, distances = distances)
  dist_max <- quantile(df_dist$distances, 0.95)
  df_dist_filter <- df_dist[df_dist$distances < dist_max, ]
  df_dist_filter$lowess <- lowess(y = df_dist_filter$distances, x = df_dist_filter$quantiles, f = lowess_frac, delta = 0)$y

  modal_quantile <- df_dist_filter[which.min(df_dist_filter$lowess), 2]

  # add values to table
  df_regression["modal_quantile"] <- modal_quantile
  df_regression["modal_curve"] <- df_regression[modal_quantile]
  df_regression["modal_corrected"] <- df_regression["reads"] / df_regression[modal_quantile]

  # Merge back the missing bins to ensure dimensions stay consistent
  df_regression <- dplyr::left_join(df_regression_all, df_regression, by = c("reads", "gc", "bin_ids"))

  # Do we want to reassign NAs as zeros?
  df_regression[is.na(df_regression$modal_corrected), "modal_corrected"] <- 0

  # What about below zero values...?
  df_regression[(df_regression$modal_corrected < 0), "modal_corrected"] <- 0

  if (results == "full") {
    return(df_regression)
  } else if (results == "default") {
    # Don't return each quantile curve
    return(df_regression[, which(!colnames(df_regression) %in% quantile_names)])
  } else if (results == "counts") {
    # only return the corrected counts
    return(as.vector(df_regression$modal_corrected))
  }
}




#' Perform GC Correction
#'
#' Performs GC correction using [BiocParallel::bplapply()] over a matrix of cell counts
#'
#' Note: If using `modal` must pass `results="counts"` as an argument.
#'
#' @param mat Count matrix for GC correction
#' @param gc GC corresponding to bins (rows) in the matrix
#' @param method Specifies the type of GC correction to perform. One of `'modal', 'copykit', or 'loess'`
#' @param ... Additional arguments to be passed to GC correction methods or [BiocParallel::bplapply()]
#'
#' @return Sparse matrix of corrected counts
#'
#' @export
perform_gc_cor <- function(mat, gc, method = c("modal", "copykit", "loess"), ...) {
  gc <- as.vector(gc)

  method <- match.arg(method)

  # If method is modal we must force result="counts" in the dot expansion
  # Not sure how to do this so just warn and exit
  if (method == "modal") {
    # Check for results="counts"
    dots <- list(...)

    if (dots["results"] != "counts") {
      stop("If using modal GC correction you must pass results='counts' as an argument")
    }
    # Ideally we would simply append this option to the dots
    # automatically but not sure how to dynamically alter the dots
    #   dots$result = "counts"
    #   (...) = dots
  }

  # Functions switch
  FUN <- switch(method,
    "modal" = get("gc_cor_modal"),
    "copykit" = get("gc_cor_copykit"),
    "loess" = get("gc_cor_loess")
  )

  # Parallel apply over the matrix
  counts_gc_list <- suppressMessages(BiocParallel::bplapply(X = seq_len(ncol(mat)), FUN = function(i) {
    FUN(gc, counts = as.vector(mat[, i]), ...)
  }))
  corrected <- as(do.call("cbind", counts_gc_list), "dgCMatrix")
  colnames(corrected) <- colnames(mat)
  rownames(corrected) <- rownames(mat)
  return(corrected)
}

# Copykit method of GC correction (Really https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3378858/) -- this ref also suggests a different model which corrects GC based on fragment length which may be better for atac data?
gc_cor_copykit <- function(counts, gc, span = 0.05) {
  counts <- as.vector(counts)
  gc_cor <- lowess(gc, log(counts + 1e-3), f = span)
  gc_cor_z <- approx(gc_cor$x, gc_cor$y, gc)
  counts_cor <- exp(log(counts) - gc_cor_z$y) * median(counts)
  return(counts_cor)
}

# Different way similar to HMMcopy
# Code adapated from `HMMcopy` package
gc_cor_hmm <- function(counts, gc, span1 = 0.03, span2 = 0.3) {
  counts <- as.vector(counts)
  rough <- loess(counts ~ gc, span = span1)
  idx <- seq(0, 1, by = 0.001)
  final <- loess(predict(rough, idx) ~ idx, span = span2)
  counts_cor <- counts / predict(final, gc)
  return(counts_cor)
}

# Modal quantile regression
# We use this function to wrap the modal_quantile_regression to just return the corrected counts. The original function returns a larger df in order to be as close to the original python function as possible
# This is now deprecated. Use gc_cor_modal
# @export
# gc_correct_modal <- function(counts, gc, ...) {
#   df <- modal_quantile_regression(counts = counts, gc = gc, ...)
#   counts_cor <- as.vector(df$modal_corrected)
#   return(counts_cor)
# }
