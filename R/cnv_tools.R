#' Smooth outlier copy number counts
#'
#' @param sce SingleCellExperiment object
#' @param assay_name Name of assay to smooth
#' @param ncores Number of cores to use
#' @param smooth_name Name of returned assay with smoothed counts
#' @param ...
#'
#' @return A SingleCellExperiument Object
#' @export
#'
smooth_counts <- function(sce, assay_name, ncores = 1, smooth_name = paste(assay_name, "smoothed", sep = "_"), ...) {
  chrs <- as.vector(seqnames(rowRanges(sce)))
  starts <- start(rowRanges(sce))
  sample_ids <- colnames(sce)
  logger::log_info("Smoothing {assay_name}")
  smoothed_counts <- pbmcapply::pbmclapply(1:ncol(sce), mc.cores = ncores, FUN = function(i) {
    x <- as.vector(assay(x = sce, assay_name)[, i])
    obj <- DNAcopy::CNA(genomdat = x, chrom = chrs, maploc = starts, data.type = "logratio", sampleid = sample_ids[i], presorted = T)
    res <- round(withr::with_seed(3, smoothed_CNA_counts <- DNAcopy::smooth.CNA(obj,
      smooth.region = 4,
      outlier.SD.scale = 4,
      smooth.SD.scale = 2,
      trim = 0.025
    ))[, 3], 2)
  })
  names(smoothed_counts) <- sample_ids
  smoothed_counts <- as.matrix(dplyr::bind_rows(smoothed_counts))
  rownames(smoothed_counts) <- rownames(sce)
  smoothed_counts[smoothed_counts <= 0] <- 1e-4 # negative or zero values post smoothed.

  assay(sce, smooth_name) <- smoothed_counts
  logger::log_success("Smoothing complete! Smoothed counts in assay '{smooth_name}'")
  return(sce)
}

#' scATAC CBS CNV segmentation
#'
#' @param sce A `SingleCellExperiment` object
#' @param assay_name Name of the assay to segment
#' @param new_assay Name of the new assay
#' @param verbose Verbosity
#' @param ncores Number of cores to use
#' @param ... Additional parameters to pass to [DNAcopy::segment]
#'
#' @inheritParams DNAcopy::segment
#'
#' @return A `SingleCellExperiment` object
#' @export
#'
segment_cnv <- function(sce, assay_name, new_assay = paste(assay_name, "segment", sep = "_"), alpha = 0.2, nperm = 10, min.width = 2, undo.splits = "none", verbose = 0, ncores = 1, ...) {
  chrs <- as.vector(seqnames(rowRanges(sce)))
  starts <- start(rowRanges(sce))
  sample_ids <- colnames(sce)
  # Perform segmentation
  # TODO: Make this function chromosome arm aware (ie segment within arms rather than chrs)
  logger::log_info("Segmenting CNVs")
  segmented_counts <- pbmcapply::pbmclapply(1:ncol(sce), mc.cores = ncores, FUN = function(i) {
    x <- as.vector(assay(sce, assay_name)[, i])
    obj <- DNAcopy::CNA(genomdat = x, chrom = chrs, maploc = starts, data.type = "logratio", sampleid = sample_ids[i], presorted = T)
    res <- withr::with_seed(3, smoothed_CNA_counts <- DNAcopy::segment(obj,
      alpha = alpha,
      nperm = nperm,
      min.width = min.width,
      undo.splits = undo.splits,
      verbose = verbose,
      ...
    ))

    # test0 <- res$segRows[[2]] + 1 - res$segRows[[1]]

    df <- data.frame(idx = 1:length(x), seg.mean = NA)

    for (j in 1:nrow(res$segRows)) {
      # Fails if no counts on final segments so we put try
      try(df[res$segRows[j, 1]:res$segRows[j, 2], "seg.mean"] <- res$output[j, "seg.mean"])
    }
    logger::log_success("Segmentation completed!")
    return(df$seg.mean)
  })

  names(segmented_counts) <- sample_ids
  segmented_counts <- as.matrix(dplyr::bind_rows(segmented_counts))
  rownames(segmented_counts) <- rownames(sce)

  assay(sce, new_assay) <- segmented_counts
  return(sce)
}


#' Merge segment levels
#'
#' Wrapper for [mergeLevels] to merge segments
#'
#' @inheritParams segment_cnv
#' @param smooth_assay name of assay with smoothed counts
#' @param segment_assay name of assay with segmented counts
#'
#' @return A `SingleCellExperiment` object
#' @export
#'
merge_segments <- function(sce, smooth_assay, segment_assay, new_assay = "segment_merged", ncores = 1) {
  smooth_counts <- log2(assay(sce, smooth_assay))

  segment_df <- as.data.frame(assay(sce, segment_assay))
  segment_df[segment_df == 0] <- 1e-4
  segment_df <- log2(segment_df)

  logger::log_info("Merging segments using {ncores} cores for {ncol(segment_df)} cells")

  seg_ml_list <- pbmcapply::pbmclapply(seq_along(segment_df), mc.cores = ncores, function(i) {
    cell_name <- names(segment_df)[i]
    seg_means_ml <- mergeLevels(
      vecObs = smooth_counts[, i],
      vecPred = segment_df[, i],
      verbose = 0,
      pv.thres = 1e-4
    )$vecMerged
  })

  names(seg_ml_list) <- names(segment_df)
  seg_ml_df <- dplyr::bind_cols(seg_ml_list)
  seg_ml_df <- round(2^seg_ml_df, 2)
  rownames(seg_ml_df) <- rownames(sce)
  assay(sce, new_assay) <- seg_ml_df

  # saving as segment ratios
  seg_ratios <- sweep(seg_ml_df, 2, apply(seg_ml_df, 2, mean, na.rm = T), "/")
  rownames(seg_ratios) <- rownames(sce)
  assay(sce, paste(new_assay, "ratios", sep = "_")) <- as.matrix(round(seg_ratios, 2))

  sce <- logNorm(sce, assay = paste(new_assay, "ratios", sep = "_"), name = paste(new_assay, "logratios", sep = "_"))

  logger::log_info("Merged segments in: {new_assay}")
  logger::log_info("Merged segments ratios in: {paste(new_assay, 'ratios', sep = '_')}")
  logger::log_info("Merged segments logratios in: {paste(new_assay, 'logratios', sep = '_')}")

  return(sce)
}


#' Identify normal cells from scATAC data
#'
#' Uses the standard deviation of the difference between each bin to estimate tumor and normal cell clusters. Using method `gmm` will calculate the per cluster median of the sd, and then fit a two-component GMM to determine tumor cells. If method `min_sd` is specified (or if only two clusters are provided), simply uses the `n_normal_clusts` with the lowest median sd.
#'
#'
#' @param sce SingleCellExperiment Object
#' @param assay_name Name of assay from which to calculate metrics from. It is highly recommended that these are segmented and merged data.
#' @param group_by Name of column containing the grouping information
#' @param method One of `gmm` or `min_sd`
#' @param n_normal_clusts Number of expected normal clusters (only for method `min_sd`)
#' @param plot Plot cluster scores and tumor/normal identifications
#' @param use_cnv_score Also use CNV score (absolute mean of the assay)
#'
#' @return SingleCellExperiment object with column `tumor_cell`
#' @export
#'
identify_normal <- function(sce, assay_name, group_by = "clusters", method = c("gmm", "min_sd"), n_normal_clusts = NULL, plot = TRUE, use_cnv_score = TRUE) {
  # TODO: Need a fallback when all clusters are very close or not sure if tumor or normal. Perhaps spike in normal distribution assumed based on reference population
  # TODO: See if we can apply this without pre-clustering

  method <- match.arg(method, choices = c("min_sd", "gmm"))

  if (length(unique(sce[[group_by]])) == 1) {
    logger::log_warn("Only one group detected. Cannot identify normal cells")
    sce$tumor_cell <- NA
    return(sce)
  }

  if (length(unique(sce[[group_by]])) == 2) {
    logger::log_warn("Only two groups detected. Defaulting to method = 'min_sd'")
    method <- "min_sd"
    n_normal_clusts <- 1
  }

  sce$seg_sd <- colSdDiffs(assay(sce, assay_name), na.rm = TRUE)


  # Take the per cluster medians
  s <- split(sce[["seg_sd"]], sce[[group_by]])

  mus <- lapply(s, median) %>% unlist()

  if (method == "gmm") {
    if (use_cnv_score) {
      sce$cnv_score <- apply(assay(sce, assay_name), MARGIN = 2, FUN = function(X) abs(mean(X, na.rm = T)))
      s2 <- split(sce[["cnv_score"]], sce[[group_by]])
      mus2 <- lapply(s2, median) %>% unlist()
      mus <- cbind(mus, mus2)
      mod <- mclust::densityMclust(mus, G = 2, plot = FALSE, verbose = FALSE)
    } else {
      mod <- mclust::densityMclust(mus, G = 2, plot = FALSE, verbose = FALSE)
    }

    # Pick the classifications with lower mus
    mus <- as.data.frame(mus)

    m1 <- mean(mus[which(mod$classification == 1), 1])
    m2 <- mean(mus[which(mod$classification == 2), 1])

    normal_class <- ifelse(test = (m1 < m2), yes = 1, no = 2)
    normal_clust <- names(which(mod$classification == normal_class))
  }

  # Simply identify grou with lowest median sequential segmental difference
  if (method == "min_sd") {
    if (is.null(n_normal_clusts)) {
      n_normal_clusts <- 1
    }
    if (n_normal_clusts >= length(unique(sce[[group_by]]))) {
      logger::log_warn("Provided n_normal_clusts = {n_normal_clusts} with {length(unique(sce[[group_by]])} clusters. Setting n_normal_clusts to {length(unique(sce[[group_by]])) - 1}.")
      n_normal_clusts <- length(unique(sce[[group_by]])) - 1
    }
    normal_clust <- names(sort(mus)[1:n_normal_clusts])
  }

  sce$tumor_cell <- FALSE

  sce$tumor_cell[which(!sce[[group_by]] %in% normal_clust)] <- TRUE

  logger::log_info("{table(sce$tumor_cell)[[1]]} normal cells identified in {length(normal_clust)} clusters using {method} method. Clusters = {paste(normal_clust, collapse =', ')}")

  if (plot) {
    p1 <- suppressWarnings(qplot(x = sce[[group_by]], y = sce[["seg_sd"]], geom = "boxplot", fill = sce[["tumor_cell"]]) + scale_fill_manual(values = col_tumor_cells()) + labs(x = paste0(group_by), y = paste0(assay_name, " cell sd"), fill = "Tumor cell"))
    if (use_cnv_score) {
      p2 <- colData(sce) %>%
        as.data.frame() %>%
        ggplot(aes(x = seg_sd, y = cnv_score, color = tumor_cell)) +
        geom_density_2d() +
        labs(x = "Cell sd", y = "CNV score", color = "Tumor cell")

      print(p1 + p2)
    } else {
      print(p1)
    }
  }
  return(sce)
}

#' @export
calc_ratios <- function(sce, assay_name, fun = c("mean", "median"), new_assay = paste(assay_name, "ratios", sep = "_")) {
  fun <- match.arg(fun)

  ratios <- sweep(assay(sce, assay_name), 2, apply(assay(sce, assay_name), 2, fun, na.rm = T), "/")
  assay(sce, new_assay) <- round(ratios, 2)
  return(sce)
}


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Merge Levels
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# mergeLevels and combine.func are present in https://doi.org/10.1038/ng.3641
# They can also be found in the package aCGH from Peter Dimitrov
# https://www.bioconductor.org/packages/release/bioc/html/aCGH.html
# under GPL-2 licence

# This is taken from copykit to avoid additional dependencies

#' @export
#' @importFrom stats ansari.test wilcox.test
mergeLevels <- function(vecObs, vecPred, pv.thres = 1e-04, ansari.sign = 0.05,
                        thresMin = 0.05, thresMax = 0.5, verbose = 1, scale = TRUE) {
  if (thresMin > thresMax) {
    cat("Error, thresMax should be equal to or larger than thresMin\n")
    return()
  }
  thresAbs <- thresMin
  sq <- numeric()
  j <- 0
  ansari <- numeric()
  lv <- numeric()
  flag <- 0
  if (thresMin == thresMax) {
    flag <- 2
  } else {
    l.step <- signif((thresMax - thresMin) / 10, 1)
    s.step <- signif((thresMax - thresMin) / 200, 1)
  }
  while (1) {
    if (verbose >= 1) {
      cat("\nCurrent thresAbs: ", thresAbs, "\n")
    }
    j <- j + 1
    sq[j] <- thresAbs
    vecPredNow <- vecPred
    mnNow <- unique(vecPred)
    mnNow <- mnNow[!is.na(mnNow)]
    cont <- 0
    while (cont == 0 & length(mnNow) > 1) {
      mnNow <- sort(mnNow)
      n <- length(mnNow)
      if (verbose >= 2) {
        cat("\r", n, ":", length(unique(vecPred)), "\t")
      }
      if (scale) {
        d <- (2 * 2^mnNow)[-n] - (2 * 2^mnNow)[-1]
      } else {
        d <- (mnNow)[-n] - (mnNow)[-1]
      }
      dst <- cbind(
        abs(d)[order(abs(d))], (2:n)[order(abs(d))],
        (1:(n - 1))[order(abs(d))]
      )
      for (i in 1:nrow(dst)) {
        cont <- 1
        out <- combine.func(
          diff = dst[i, 1], vecObs,
          vecPredNow, mnNow, mn1 = mnNow[dst[i, 2]],
          mn2 = mnNow[dst[i, 3]], pv.thres = pv.thres,
          thresAbs = if (scale) {
            2 * 2^thresAbs - 2
          } else {
            thresAbs
          }
        )
        if (out$pv > pv.thres) {
          cont <- 0
          vecPredNow <- out$vecPredNow
          mnNow <- out$mnNow
          break
        }
      }
    }
    ansari[j] <- ansari.test(sort(vecObs - vecPredNow), sort(vecObs -
      vecPred))$p.value
    if (is.na(ansari[j])) {
      ansari[j] <- 0
    }
    lv[j] <- length(mnNow)
    if (flag == 2) {
      break
    }
    if (ansari[j] < ansari.sign) {
      flag <- 1
    }
    if (flag) {
      if (ansari[j] > ansari.sign | thresAbs == thresMin) {
        break
      } else {
        thresAbs <- signif(thresAbs - s.step, 3)
        if (thresAbs <= thresMin) {
          thresAbs <- thresMin
        }
      }
    } else {
      thresAbs <- thresAbs + l.step
    }
    if (thresAbs >= thresMax) {
      thresAbs <- thresMax
      flag <- 2
    }
  }
  return(list(
    vecMerged = vecPredNow, mnNow = mnNow, sq = sq,
    ansari = ansari
  ))
}


#' @export
combine.func <- function(diff, vecObs, vecPredNow, mnNow, mn1, mn2, pv.thres = 1e-04,
                         thresAbs = 0) {
  vec1 <- vecObs[which(vecPredNow == mn1)]
  vec2 <- vecObs[which(vecPredNow == mn2)]
  if (diff <= thresAbs) {
    pv <- 1
  } else {
    if ((length(vec1) > 10 & length(vec2) > 10) | sum(
      length(vec1),
      length(vec2)
    ) > 100) {
      pv <- wilcox.test(vec1, vec2)$p.value
    } else {
      pv <- wilcox.test(vec1, vec2, exact = TRUE)$p.value
    }
    if (length(vec1) <= 3 | length(vec2) <= 3) {
      pv <- 0
    }
  }
  index.merged <- numeric()
  if (pv > pv.thres) {
    vec <- c(vec1, vec2)
    index.merged <- which((vecPredNow == mn1) | (vecPredNow ==
      mn2))
    vecPredNow[index.merged] <- median(vec, na.rm = TRUE)
    mnNow[which((mnNow == mn1) | (mnNow == mn2))] <- median(vec,
      na.rm = TRUE
    )
    mnNow <- unique(mnNow)
  }
  list(mnNow = mnNow, vecPredNow = vecPredNow, pv = pv)
}


#' logNorm()
#'
#' Computes a log transformation of the selected assay. Function taken from the `copyKit` package.
#'
#' @param scCNA scCNA object.
#' @param transform String specifying the transformation to apply to the selected
#' assay.
#' @param assay String with the name of the assay to pull data from to run the
#' segmentation.
#' @param name String with the name for the target slot for the resulting
#' transformed counts.
#'
#' @return A data frame with log transformed counts inside the `assay` slot.
#'
#' @export
logNorm <- function(scCNA,
                    transform = c("log", "log2", "log10", "log1p"),
                    assay = "segment_ratios",
                    name = "logr") {
  transform <- match.arg(transform)

  # obtaining data
  seg_ratios <- assay(scCNA, assay)

  # saving logr
  seg_ratios[seg_ratios == 0] <- 1e-3

  if (transform == "log") {
    seg_ratios_logr <- log(seg_ratios)
  } else if (transform == "log2") {
    seg_ratios_logr <- log2(seg_ratios)
  } else if (transform == "log1p") {
    seg_ratios_logr <- log1p(seg_ratios)
  } else if (transform == "log10") {
    seg_ratios_logr <- log10(seg_ratios)
  }

  assay(scCNA, name) <- round(seg_ratios_logr, 2)

  return(scCNA)
}
