#' Integrate segments between assays
#'
#' Will merge bins of `y` into `x`, taking the weighted mean of the binned signal in `y` of `granges_signal_colname`
#'
#' @param x granges object to integrate on
#' @param y granges
#' @param granges_signal_colname column containing cna signal in granges
#' @param drop_na Logical. Drop bins with any NA values.
#'
#' @return `x` with averaged signal of `y` integrated over the bins
#' @export
#'
integrate_segments <- function(x, y, granges_signal_colname, drop_na = TRUE) {

  # Check and fix names if needed (in case the sample names are integer or invalid colnames)
  # This fixes the error that can occur during df conversion
  orig_colnames <- granges_signal_colname
  granges_signal_colname <- make.names(granges_signal_colname)

  # Ensure we use fixed column names from the outset
  colnames(mcols(y))[match(orig_colnames, colnames(mcols(y)))] <- granges_signal_colname

  # Make sure chrs match to the input
  df_x <- as.data.frame(x)
  df_y <- as.data.frame(y)

  if (all(grepl("chr", df_x$seqnames))) {
    no_chr <- FALSE
  } else {
    no_chr <- TRUE
  }

  # If x is no_chr and y is not
  if (no_chr & all(grepl("chr", df_y$seqnames))) {
    # Need to remove chr names
    df_y$seqnames <- gsub("chr", "", df_y$seqnames)
    df_y$seqnames <- gsub("X", "23", df_y$seqnames)
    df_y$seqnames <- gsub("Y", "24", df_y$seqnames)
  }
  # If x has chr and y does not
  if (!no_chr & !all(grepl("chr", df_y$seqnames))) {
    # Need to remove chr names
    df_y$seqnames <- gsub("23", "X", df_y$seqnames)
    df_y$seqnames <- gsub("24", "Y", df_y$seqnames)
    df_y$seqnames <- paste0("chr", df_y$seqnames)
  }


  gr.comb <- bind_rows(df_x, df_y) %>% makeGRangesFromDataFrame(keep.extra.columns = TRUE)

  x <- makeGRangesFromDataFrame(df_x, keep.extra.columns = TRUE)
  y <- makeGRangesFromDataFrame(df_y, keep.extra.columns = TRUE)

  # Disjoin all overlapping ranges
  gr.comb.disjoin <- disjoin(gr.comb)

  # Attach the CNA signal to respective bins
  olaps2 <- findOverlaps(gr.comb.disjoin, y)
  mcols(gr.comb.disjoin)[granges_signal_colname] <- NA
  mcols(gr.comb.disjoin[queryHits(olaps2)])[granges_signal_colname] <- mcols(y[subjectHits(olaps2)])[granges_signal_colname]

  gr.comb.disjoin$binwidth <- width(gr.comb.disjoin)

  # reaggregate the signal over set tile width of the input
  hits <- findOverlaps(x, gr.comb.disjoin)

  df <- as.data.frame(mcols(gr.comb.disjoin)[granges_signal_colname])

  res <- lapply(unique(queryHits(hits)), FUN = function(idx) {
    # Get the index in the matrix
    mat_idx <- subjectHits(hits)[which(queryHits(hits) == idx)]

    slice_binwidths <- gr.comb.disjoin$binwidth[mat_idx]

    # Apply the weighted mean calculation
    slice_means <- apply(as.matrix(df[mat_idx, ]), MARGIN = 2, FUN = function(x) {
      weighted.mean(x, w = slice_binwidths, na.rm = TRUE)
    })
  })

  scores_df <- data.frame(do.call("rbind", res))

  colnames(scores_df) <- orig_colnames

  mcols(x)[[orig_colnames]] <- scores_df[[orig_colnames]]

  # Drop bins that have any NA values
  if (drop_na) {
    na_idx <- unique(which(is.na(mcols(x)[orig_colnames]), arr.ind = TRUE)[, 1])
    x <- x[-na_idx, ]
  }

  return(x)

  # This version get's the maximum overlaps
  # # Remove ranges that are NA for either the modal total CN or the clone data
  # na_rows <- which(is.na(mcols(gr.comb.disjoin)[c(granges_signal_colname, group_columns)]), arr.ind = TRUE)[, 1]
  #
  # atac_comb <- gr.comb.disjoin[-na_rows]
  # atac_comb$binwidth <- width(atac_comb)
  #
  # # convert this to a singlecell experiment object for ease
  # atac_comb_sce <- SingleCellExperiment(assays = list(counts = as.matrix(mcols(atac_comb)[group_columns])), rowRanges = atac_comb)
  # colData(atac_comb_sce) <- colData(sce)

  # return(list(df = atac_comb, sce = atac_comb_sce))
}
