#' Rebin `SingleCellExperiment` object
#'
#' Will rebin data into larger bins. Have not tested with smaller bins.
#'
#' @param sce `SingleCellExperiment` object
#' @param assays Names of assays to include in the rebinning process
#' @param new_bins GenomicRanges of new bins
#' @param cell_ids Cells to include in rebinning. If `NULL` defaults to all cells
#' @param ncores Number of cores to speed up rebinning
#'
#' @return Rebinned `SingleCellExperiment` object
#' @export
#'
rebin_sce <- function(sce, assays, new_bins, cell_ids = NULL, ncores = 1) {
  if (is.null(cell_ids)) {
    logger::log_info("Cell ids not provided. Rebinning all cells")
    cell_ids <- colnames(sce)
  }

  # Get the mapping between bins
  bin_map <- disjoint_bins_map(rowRanges(sce), new_bins)

  # New ranges
  new_ranges <- new_bins[unique(bin_map$y_bins)]

  new_assays <- vector(mode = "list")

  # External loop over assays
  for (assay_cur in assays) {
    message(glue("{assay_cur}"))
    # Loop over the new bins and map the weighted mean values on
    res <- pbmcapply::pbmclapply(unique(bin_map$y_bins), mc.cores = ncores, FUN = function(new_bin_idx) {
      # res <- lapply(unique(bin_map$y_bins), FUN = function(new_bin_idx) {
      # logger::log_debug("{new_bin_idx}")
      # Original bin indexes and disjoint widths
      orig_bin_idx <- bin_map[bin_map$y_bins == new_bin_idx, "x_bins"]
      orig_bin_widths <- bin_map[bin_map$y_bins == new_bin_idx, "disjoint_width"]
      # Now loop over cells and perform the weighted mean calculations
      # Catch single gene bins
      if (length(orig_bin_widths) == 1) {
        mat <- t(as.matrix(assay(sce, assay_cur)[orig_bin_idx, cell_ids]))
      } else {
        mat <- assay(sce, assay_cur)[orig_bin_idx, cell_ids]
      }
      slice_means <- apply(mat, MARGIN = 2, FUN = function(x) {
        weighted.mean(x, w = orig_bin_widths, na.rm = TRUE)
      })
    })

    new_assays[[assay_cur]] <- as.matrix(do.call("rbind", res))
  }

  # Create the new single cell experiment object and return
  sce_rebin <- SingleCellExperiment(new_assays, rowRanges = new_ranges, colData = colData(sce))
  rownames(sce_rebin) <- get_bin_ids(rowRanges(sce_rebin))

  return(sce_rebin)
}

#' Map between bins
#'
#' Will return a mapping between `x` and `y` bins using the `disjoin` function. The results can be used to map bins together taking weighted mean of the `disjoint_width`
#'
#' @param x x Granges
#' @param y y Granges
#'
#' @return data frame with bin mapping
#' @export
#'
disjoint_bins_map <- function(x, y) {
  # Split out just the bins to get the overlaps
  x_ranges <- GenomicRanges::granges(x)
  y_ranges <- GenomicRanges::granges(y)

  # Make sure chrs match to the input
  if (all(grepl("chr", seqnames(x_ranges)))) {
    no_chr <- FALSE
  } else {
    no_chr <- TRUE
  }

  # If x is no_chr and y is not
  if (no_chr & all(grepl("chr", seqlevels(y_ranges)))) {
    # Need to remove chr names
    seqlevels(y_ranges) <- gsub("chr", "", seqlevels(y_ranges))
    seqlevels(y_ranges) <- gsub("X", "23", seqlevels(y_ranges))
    seqlevels(y_ranges) <- gsub("Y", "24", seqlevels(y_ranges))
  }
  # If x has chr and y does not
  if (!no_chr & !all(grepl("chr", seqlevels(y_ranges)))) {
    # Need to remove chr names
    seqlevels(y_ranges) <- gsub("23", "X", seqlevels(y_ranges))
    seqlevels(y_ranges) <- gsub("24", "Y", seqlevels(y_ranges))
    seqlevels(y_ranges) <- paste0("chr", seqlevels(y_ranges))
  }

  # Combine the ranges
  gr.comb <- c(x_ranges, y_ranges)

  # Disjoin all overlapping ranges
  gr.comb.disjoin <- disjoin(gr.comb)
  gr.comb.disjoin$binwidth <- width(gr.comb.disjoin)

  # 'Attach' the CNA signal to respective bins
  # This overlap object will contain the matching y bins (in subjecthits) to the
  # query disjoint bins
  olaps2 <- data.frame(findOverlaps(gr.comb.disjoin, x_ranges)) %>%
    dplyr::rename(disjoint_bins = queryHits, x_bins = subjectHits)

  hits <- data.frame(findOverlaps(y_ranges, gr.comb.disjoin)) %>%
    dplyr::rename(y_bins = queryHits, disjoint_bins = subjectHits)

  df <- dplyr::full_join(hits, olaps2, by = "disjoint_bins", copy = TRUE, multiple = "all") %>%
    select(x_bins, y_bins, disjoint_bins)

  # Get the width of the disjoint bins and remove NA
  df$disjoint_width <- gr.comb.disjoin$binwidth[df$disjoint_bins]
  df <- na.omit(df)

  return(df)
}
