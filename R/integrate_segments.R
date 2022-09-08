#' Integrate segments between assays
#'
#' Will merge bins of `y` into `x`, taking the weighted mean of the binned signal in `y` of `granges_signal_colname`
#'
#' @param x granges object to integrate on
#' @param y granges
#' @param granges_signal_colname column containing cna signal in granges
#'
#' @return `x` with averaged signal of `y` integrated over the bins
#' @export
#'
integrate_segments <- function(x, y, granges_signal_colname) {
  # Make sure chrs match to the input
  df_x <- as.data.frame(x)
  df_y <- as.data.frame(y)

  if (all(grepl("chr", df_x$seqnames))) {
    no_chr = FALSE
  } else {
    no_chr = TRUE
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
  # EXPERIMENTAL: See if we can reaggregate the signal over set tile width of the input
  hits <- findOverlaps(x, gr.comb.disjoin)

  # Aggregate the mean score over the original bins
  # Take the binwidth weighted mean across the bins
  grange_agg_score <- aggregate(gr.comb.disjoin, hits,
                                score = weighted.mean(x = eval(as.symbol(granges_signal_colname)), w = binwidth, na.rm = TRUE))

  mcols(x)[granges_signal_colname] <- grange_agg_score$score

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
