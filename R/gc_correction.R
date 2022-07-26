# TODO: Also code to regenerate mappability per bins from the hg38 hoffman bigwigs
#
# http://hgdownload.soe.ucsc.edu/gbdb/hg38/hoffmanMappability/
# https://bismap.hoffmanlab.org/
#


#' Add GC correction to SCE Object
#'
#' @param sce an SCE object
#' @param assay_name Name of assay on which to perform GC correction
#'
#' @inheritParams perform_gc_cor
#'
#' @return sce object with corrected GC count matrix in `assay(sce, 'counts_gc_[method]')`. See
#' @export
#'
add_gc_cor <- function(sce, gc = rowData(sce)$gc, assay_name = "counts", method = c("modal", "copykit", "loess"), verbose = FALSE, ncores = 1, ...) {

  method <- match.arg(method)

  gc_slot <- paste0("counts_gc_", method)

  # Check if valid bins exists and pass correctly
  if ("valid_bins" %in% names(assays(sce))) {
    if (verbose) {logger::log_info("Found valid bins in sce object")}
    valid_mat <- assay(sce, "valid_bins")
  } else {
    warning("No valid bins matrix in sce object. Defaulting to all valid.")
    valid_mat <- NULL
  }

  assay(sce, gc_slot) <- perform_gc_cor(
    mat = assay(sce, assay_name),
    gc = gc,
    valid_mat = valid_mat,
    method = method,
    ncores = ncores,
    verbose = verbose,
    ...
  )
  metadata(sce)$gc_cor_method <- method

  return(sce)
}


#' Perform GC Correction
#'
#' Performs GC correction using over a matrix of cell counts
#'
#' Note: If using `modal` must pass `results="counts"` as an argument. Uses [pbmcapply::pbmclapply()] for parallelization.
#'
#' @param mat Count matrix for GC correction
#' @param gc GC corresponding to bins (rows) in the matrix
#' @param valid_mat Matrix of TRUE/FALSE for valid bins. If none provided defaults to all TRUE
#' @param method Specifies the type of GC correction to perform. One of `'modal', 'copykit', or 'loess'`
#' @param ncores Number of cores to use if parallel backend is available
#' @param verbose Message verbosity (TRUE/FALSE)
#' @param ... Additional arguments to be passed to GC correction methods
#'
#' @return Sparse matrix of corrected counts
#'
#' @export
perform_gc_cor <- function(mat, gc, valid_mat = NULL, method = c("modal", "copykit", "loess"), ncores = 1, verbose = FALSE, ...) {
  if (verbose) {
    logger::log_info("Performing GC correction on {ncol(mat)} cells using {ncores} threads.")
    logger::log_info("GC correction method: {method}")
  }

  # Just pass all as true if not provided
  if (is.null(valid_mat)) {
    warning("No valid bin matrix provided. Defaulting to all bins valid. Could lead to erroneous results")
    valid_mat <- matrix(TRUE, nrow = nrow(mat), ncol = ncol(mat), dimnames = list(rownames(mat), colnames(mat)))
  }

  gc <- as.vector(gc)

  method <- match.arg(method)

  # Functions switch
  FUN <- switch(method,
    "modal" = get("gc_cor_modal"),
    "copykit" = get("gc_cor_copykit"),
    "loess" = get("gc_cor_loess")
  )

  # Parallel apply over the matrix
  if (requireNamespace("pbmcapply", quietly = TRUE)) {
    counts_gc_list <- pbmcapply::pbmclapply(X = seq_len(ncol(mat)), mc.cores = ncores, FUN = function(i) {
      FUN(gc = gc, counts = as.vector(mat[, i]), valid = as.vector(valid_mat[, i]), bin_ids = rownames(mat), ...)
    })
  } else {
    warning("No parallel backend used. GC correction may be slow.", call. = FALSE)
    counts_gc_list <- lapply(X = seq_len(ncol(mat)), FUN = function(i) {
      FUN(gc = gc, counts = as.vector(mat[, i]), valid = as.vector(valid_mat[, i]), bin_ids = rownames(mat), ...)
    })
  }

  if (verbose) {
    logger::log_success("GC correction completed!")
  }


  corrected <- matrix(unlist(counts_gc_list), ncol = length(counts_gc_list), byrow = FALSE)

  colnames(corrected) <- colnames(mat)
  rownames(corrected) <- rownames(mat)
  return(corrected)
}



#' Modal regression GC Correction for single cell data
#'
#' @param counts Vector of counts in bins
#' @param gc Vector of GC content in bins
#' @param valid Vector of values corresponding to valid bins. If none provided defaults to all TRUE
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
gc_cor_modal <- function(counts, gc, valid = rep(TRUE, length(counts)), bin_ids = names(counts), lowess_frac = 0.2, q = c(0.1, 0.9), g = c(0.1, 0.9), results = c("counts", "default", "full")) {
  if (length(counts) != length(gc)) {
    stop("Length of counts and gc vectors are not identical")
  }

  if (length(q) != 2 | (q[1] > q[2])) {
    stop("Must provide lower and upper quantile")
  }

  if (is.null(bin_ids)) {
    bin_ids = seq_len(length(counts))
  }

  results <- match.arg(results)

  # Into a dataframe pre-filtering
  df_regression_all <- data.frame(reads = counts, gc = gc, bin_ids = bin_ids, valid = valid)

  # Filter the data for non-zero gc and read counts
  df_regression <- df_regression_all[df_regression_all$valid == TRUE, ]

  # 2nd order polynomial quantile regression
  quantiles <- seq(q[1], q[2], by = 0.01)
  quantile_names <- paste0("q", quantiles * 100)

  if (nrow(df_regression) < 10) {
    df_regression <- .gc_warn_error(df_regression, quantile_names)
  } else {
    # Fit second order polynormial quantile regression
    poly2_quantile_model <- quantreg::rq(reads ~ gc + I(gc^2), tau = quantiles, data = df_regression)

    # Fit to our data
    poly_quantile_predict <- predict(object = poly2_quantile_model, newdata = df_regression)

    colnames(poly_quantile_predict) <- quantile_names # rename columns

    # Bind to our data
    df_regression <- dplyr::bind_cols(df_regression, poly_quantile_predict)

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

    # Error catch when df_dist_filter is empty
    if (nrow(df_dist_filter) == 0) {
      df_regression <- .gc_warn_error(df_regression, quantile_names)
    } else {
      df_dist_filter$lowess <- lowess(y = df_dist_filter$distances, x = df_dist_filter$quantiles, f = lowess_frac, delta = 0)$y

      modal_quantile <- df_dist_filter[which.min(df_dist_filter$lowess), 2]

      # add values to table
      df_regression["modal_quantile"] <- modal_quantile
      df_regression["modal_curve"] <- df_regression[modal_quantile]
      df_regression["modal_corrected"] <- df_regression["reads"] / df_regression[modal_quantile]
    }
  }



  # Below zeroes are NA? Not sure how to handle exactly
  df_regression[(df_regression$modal_corrected < 0 | is.na(df_regression$modal_corrected)), "modal_corrected"] <- NA

  # Merge back the missing bins to ensure dimensions stay consistent
  df_regression <- dplyr::left_join(df_regression_all, df_regression, by = c("reads", "gc", "bin_ids", "valid"))

  # Do we want to reassign NAs as zeros?
  # df_regression[is.na(df_regression$modal_corrected), "modal_corrected"] <- 0

  if (results == "full") {
    return(df_regression)
  } else if (results == "default") {
    # Only return selected quantile curve
    return(df_regression[, which(!colnames(df_regression) %in% quantile_names)])
  } else if (results == "counts") {
    # only return the corrected counts
    return(as.vector(df_regression$modal_corrected))
  }
}

# Warning for modal regression
.gc_warn_error <- function(df_regression, quantile_names) {
  warning("Not enough data points for modal GC regression", call. = FALSE)
  # Prepare for NA table
  df_regression[quantile_names] <- NA
  df_regression["modal_quantile"] <- NA
  df_regression["modal_curve"] <- NA
  df_regression["modal_corrected"] <- NA
  return(df_regression)
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
