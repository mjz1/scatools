#' Calculate allelic imbalance
#'
#' @param ref_counts vector of ref counts
#' @param alt_counts vector of alt counts
#'
#' @return measure of allelic imbalance
#' @export
#'
#' @examples
#' calc_ai(c(4, 2, 1), c(2, 3, 1))
calc_ai <- function(ref_counts, alt_counts) {

  # coerce to vectors
  ref_counts <- as.vector(ref_counts)
  alt_counts <- as.vector(alt_counts)

  x <- pmax(ref_counts, alt_counts)
  y <- pmin(ref_counts, alt_counts)
  ai <- (x - y) / x

  return(ai)
}


