#' Calculate CNV Score
#'
#' @param sce SingleCellExperiment object
#' @param assay_name Name of assay
#'
#' @return SCE object with `cnv_score` appended
#' @export
#'
calc_cnv_score <- function(sce, assay_name = "counts") {

  dat <- assay(sce, assay_name)

  cnv_scores <- unlist(lapply(X = seq_len(ncol(dat)), FUN = function(i) {
    cnv_score <- mean(abs(dat[,i]))
  }))
  sce$cnv_score <- cnv_scores
  return(sce)
}

#' @export
scale_sub <- function(sce, assay_name = "counts", log2 = FALSE, scale = c("none", "cells", "bins", "both"), verbose = FALSE) {
  mat <- assay(sce, assay_name)

  scaled_mat <- scale_mat(mat, log2 = log2, scale = scale)

  new_assay <- paste0(assay_name, "_scaled_", scale)

  sce <- sce[rownames(scaled_mat), colnames(scaled_mat)]

  assay(sce, new_assay) <- scaled_mat
  if (verbose) {logger::log_info("Scaled assay: {new_assay}")}

  return(sce)
}


#' @export
scale_mat <- function(mat, log2 = FALSE, scale = c("none", "cells", "bins", "both")) {
  mat <- as.matrix(mat)

  scale <- match.arg(scale)

  # TODO: separate out cleaning of the matrix from this function

  logger::log_debug("Scaling: {scale}")

  if (scale == "none") {
    scale <- FALSE
  }

  # Remove fully NA or 0 columns
  keep_bins <- apply(mat, 1, FUN = function(x) !all(is.na(x)) & !all(x == 0))

  logger::log_debug("Keeping {sum(keep_bins)} of {nrow(mat)} bins")

  mat_names <- colnames(mat)

  mat <- as.matrix(mat[keep_bins, ])

  colnames(mat) <- mat_names

  # Replace remaining NAs with 0?
  mat[is.na(mat)] <- 0

  if (log2) {
    mat <- log2(mat + 1e-5)
  }

  if (scale == "both") {
    mat <- scale(t(scale(t(mat))))
  }

  if (scale == "cells") {
    mat <- scale(mat)
  }

  if (scale == "bins") {
    mat <- t(scale(t(mat)))
  }
  return(as.matrix(mat))
}


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
#' sce <- get_snp_counts(sce,
#'   variables = c("all", "Sample", "Condition"),
#'   target_assays = c("ref", "alt")
#' )
#' }
get_snp_counts <- function(sce, variables = "all", target_assays = c("ref", "alt")) {
  stopifnot(class(sce) %in% c("SingleCellExperiment"))

  if (!is.null(target_assays)) {
    # Check that all assays are present
    if (!all(target_assays %in% names(assays(sce)))) {
      missing_assays <- target_assays[which(!target_assays %in% names(assays(sce)))]
      logger::log_warn("Missing assay(", missing_assays, ") from input sce.")
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
        logger::log_warn("Missing variable '", variable, "' from input colData.")
        next
      }
    }
    for (target in target_assays) {
      # If variable is "all", compute sums over all samples
      if (variable == "all") {
        dat <- data.frame("all" = as.integer(Matrix::rowSums(assay(sce, target))))
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


.save_to <- function(object, save_to = "./", verbose = TRUE) {
  dir.create(dirname(save_to), recursive = TRUE, showWarnings = FALSE)
  if (verbose) {
    logger::log_info("Saving {deparse(substitute(object))} to \"{save_to}\"")
  }
  save(object, file = save_to)
}

#' @export
save_to <- .save_to

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

#' Bind sublists within lists
#'
#' Helper function for lists of lists, where the top level list is over multiple samples, and the sublists are results that are identical across samples that one wishes to bind
#'
#' @param toplist List of lists of `rbind`able data
#' @param sublist Name or index of sublist
#' @param what Specify either `cbind` or `rbind` (Currently only rbind implemented)
#' @param .add_id Add identifiers to each data entry prior to merging. Based on list names if available. (default=FALSE)
#' @param .id_name Name of the id column if added
#'
#' @export
#'
#' @examples
#' toplist <- list(
#'   sample_1 = list(
#'     result1 = data.frame(
#'       a = c(1, 2, 3),
#'       b = c("X", "Y", "Z")
#'     ),
#'     result2 = data.frame(
#'       height = 180,
#'       weight = 75
#'     )
#'   ),
#'   sample_2 = list(
#'     result1 = data.frame(
#'       a = c(6, 5, 4),
#'       b = c("A", "B", "C")
#'     ),
#'     result2 = data.frame(
#'       height = 155,
#'       weight = 60
#'     )
#'   )
#' )
#'
#' bind_sublist(toplist, sublist = 1, what = "rbind", .add_id = TRUE)
#'
#' bind_sublist(toplist, sublist = "result2", what = "rbind", .add_id = TRUE)
#'
#' bind_sublist(toplist, sublist = 2, what = "rbind", .add_id = FALSE)
#'
bind_sublist <- function(toplist, sublist, what = c("rbind"), .add_id = FALSE, .id_name = "id") {
  what <- match.arg(what)

  if (is.null(names(toplist))) {
    names(toplist) <- 1:length(toplist)
  }

  res <- do.call(what, lapply(names(toplist), FUN = function(name) {
    dat <- toplist[[name]][[sublist]]

    if (what == "rbind" & .add_id) {
      dat <- cbind(id = name, dat)
      colnames(dat)[1] <- .id_name
    }

    return(dat)
  }))

  return(res)
}

#' @export
getmode <- function(v) {
  uniqv <- unique(v)
  uniqv[which.max(tabulate(match(v, uniqv)))]
}

#' @export
prettyMb <- function(x, places = 3) {
  power <- pmin(6, floor(log(abs(x), 1000)))
  units <- c("B", "Kb", "Mb", "Gb", "Tb", "Pb", "Eb")[power + 1]
  x <- x / (1000^power)
  paste(prettyNum(signif(x, places)), units, sep = "")
}


chr_reorder <- function(chrs) {
  gtools::mixedsort(chrs)
}
