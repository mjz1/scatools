# Need function to define bins in cells as ideal or valid

#' Modal regression GC Correction
#'
#' @param counts Vector of counts in bins
#' @param gc Vector of GC content in bins
#' @param bin_ids Bin ids
#' @param lowess_frac The smoother span. See [stats::lowess()] for more details
#' @param q A tuple of quantile bounds to compute quantile regressions over read counts
#' @param g GC quantiles over which to integrate regression curves
#' @param full_return Return each individual quantile curve in the return df (default: FALSE)
#'
#' @return A dataframe with the modal quantile and GC corrected counts
#' @export
#'
#' @examples
modal_quantile_regression <- function(counts, gc, bin_ids = seq_along(counts), lowess_frac = 0.2, q = c(0.1, 0.9), g = c(0.1, 0.9), full_return = FALSE) {

  if (length(counts) != length(gc)) {
    stop("Length of counts and gc vectors are not identical")
  }

  if (length(q) != 2 | (q[1] > q[2])) {
    stop("Must provide lower and upper quantile")
  }

  # Into a dataframe pre-filtering
  df_regression_all <- data.frame(reads = counts, gc = gc, bin_ids = bin_ids)

  # Filter the data for non-zero gc and read counts
  df_regression = df_regression_all[df_regression_all$gc > 0,]

  df_regression = df_regression_all[df_regression_all$reads > 0,]

  # 2nd order polynomial quantile regression
  quantiles = seq(q[1], q[2], by = 0.01)
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

  poly2_quantile_params = as.data.frame(poly2_quantile_model$coefficients)
  colnames(poly2_quantile_params) <- quantile_names # rename columns

  # integration and mode selection
  gc_min = quantile(df_regression$gc, g[1])
  gc_max = quantile(df_regression$gc, g[2])

  poly2_quantile_integration <- c(0, apply(X = poly2_quantile_params, MARGIN = 2, FUN = function(params) {
    poly2 = polynom::polynomial(params)
    integ = polynom::integral(poly2)
    integrand = predict(integ, gc_max) - predict(integ, gc_min)
  }))

  # find the modal quantile
  distances = diff(poly2_quantile_integration)

  df_dist <- data.frame(quantiles = quantiles, quantile_names = quantile_names, distances = distances)
  dist_max <- quantile(df_dist$distances, 0.95)
  df_dist_filter <- df_dist[df_dist$distances < dist_max,]
  df_dist_filter$lowess <- lowess(y = df_dist_filter$distances, x = df_dist_filter$quantiles, f = lowess_frac, delta = 0)$y

  modal_quantile <- df_dist_filter[which.min(df_dist_filter$lowess),2]

  # add values to table
  df_regression['modal_quantile'] = modal_quantile
  df_regression['modal_curve'] = df_regression[modal_quantile]
  df_regression['modal_corrected'] = df_regression['reads'] / df_regression[modal_quantile]

  # Merge back the missing bins to ensure dimensions stay consistent
  df_regression <- dplyr::left_join(df_regression_all, df_regression, by = c("reads", "gc", "bin_ids"))

  # Do we want to reassign NAs as zeros?
  df_regression[is.na(df_regression$modal_corrected),'modal_corrected'] <- 0

  # What about below zero values...?
  df_regression[(df_regression$modal_corrected < 0),'modal_corrected'] <- 0

  if (full_return) {
    return(df_regression)
  } else {
    # Don't return each quantile curve
    return(df_regression[,which(!colnames(df_regression) %in% quantile_names)])
  }

}
