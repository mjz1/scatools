#' get_snp_counts
#'
#' Computes per feature counts across the dataset. Equivalent to pseudobulk summarization.
#'
#' @param sce `SingleCellExperiment` object with one or more assays
#' @param variables Factors upon which to split the counts. Defaults to `'all'` which counts across entire dataset
#' @param target_assays Assays in `sce` to count
#'
#' @return A `SingleCellExperiment` object with named counts for each assay in the `rowData(sce)` slot.
#' @export
#'
#' @examples
#' \dontrun{
#' sce <- get_snp_counts(sce, variables = c("all", "Sample", "Condition"), target_assays = c("ref", "alt"))
#' }
get_snp_counts <- function(sce, variables = "all", target_assays = c("ref", "alt")) {
  stopifnot(class(sce) %in% c("SingleCellExperiment"))

  if (!is.null(target_assays)) {
    # Check that all assays are present
    if (!all(target_assays %in% names(assays(sce)))) {
      missing_assays <- target_assays[which(!target_assays %in% names(assays(sce)))]
      warning("Missing assay(", missing_assays, ") from input sce.")
      target_assays <- target_assays[which(target_assays %in% names(assays(sce)))]
    }
  } else {
    # Set default assay if null
    target_assays <- 1
  }

  # Loop over variables and compute counts
  for (variable in variables) {
    # If variable is not 'all' and also missing, warn and skip
    if (!is.null(variable)) {
      if (is.null(sce[[variable]]) & (variable != "all")) {
        warning("Missing variable '", variable, "' from input colData.")
        next
      }
    }
    for (target in target_assays) {
      # If variable is "all", compute sums over all samples
      if (variable == "all") {
        dat <- data.frame("all" = as.integer(rowSums(assay(sce, target))))
      } else {
        dat <- sapply(levels(sce[[variable]]), FUN = function(X) {
          idx <- which(sce[[variable]] == X)
          as.integer(rowSums(assay(sce[, idx], target)))
        })
      }

      colnames(dat) <- paste(paste0(target, "Depth"), colnames(dat), sep = "_")

      # Write the data to the columns
      rowData(sce)[colnames(dat)] <- dat
    }
  }
  return(sce)
}


#' Get factor indices
#'
#' @param f a factor
#'
#' @return a named list of indices of each factor level
#' @export
#'
get_f_idx <- function(f) {
  f <- as.factor(f)
  l <- lapply(seq_along(levels(f)), FUN = function(i) {
    which(f == levels(f)[i])
  })
  names(l) <- levels(f)
  return(l)
}
