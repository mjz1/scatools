state_cn_colors <- function() {
  cn_colors <- structure(
    c(
      "#3182BD", "#9ECAE1", "#CCCCCC", "#FDCC8A", "#FC8D59", "#E34A33",
      "#B30000", "#980043", "#DD1C77", "#DF65B0", "#C994C7", "#D4B9DA"
    ),
    names = c("0", "1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11+")
  )
  return(cn_colors)
}


#' @export
#' @noRd
logr_col_fun <- function(breaks = c(-2, 0, 2), colors = c("blue", "white", "red")) {
  circlize::colorRamp2(breaks = breaks, colors = colors)
}

#' @export
#' @noRd
counts_col_fun <- function(breaks = c(0, 2, 8), colors = c("blue", "white", "red")) {
  circlize::colorRamp2(breaks = breaks, colors = colors)
}

#' @export
#' @noRd
col_tumor_cells <- function() {
  c(`TRUE` = "#5E5E5E", `FALSE` = "#FF9A85")
}
