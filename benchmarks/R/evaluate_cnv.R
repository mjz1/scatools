# Dataset-agnostic concordance metrics: scatools CN vs ground truth.
#
# The truth is provided as a GRanges carrying a copy-number signal column; this
# harness rebins it onto the scatools query bins (reusing scatools::integrate_segments)
# and computes correlation + gain/loss classification concordance.

#' Pseudobulk (or clone-level) per-bin profile from a scatools SCE.
#'
#' @param sce scatools SingleCellExperiment
#' @param assay_name assay holding relative CN / logratios (e.g. "logr_modal")
#' @param group_by optional colData column; if given, returns one column per group
#' @return GRanges (rowRanges of sce) with a `signal` column (or one per group)
pseudobulk_profile <- function(sce, assay_name, group_by = NULL) {
  gr <- SummarizedExperiment::rowRanges(sce)
  m <- SummarizedExperiment::assay(sce, assay_name)
  if (is.null(group_by)) {
    S4Vectors::mcols(gr)$signal <- Matrix::rowMeans(m, na.rm = TRUE)
  } else {
    ids <- sce[[group_by]]
    for (g in unique(ids)) {
      S4Vectors::mcols(gr)[[paste0("signal_", g)]] <-
        Matrix::rowMeans(m[, ids == g, drop = FALSE], na.rm = TRUE)
    }
  }
  gr
}

#' Rebin a truth CN GRanges onto query bins (weighted mean over overlaps).
#'
#' Dogfoods scatools::integrate_segments.
#' @param query_bins GRanges of scatools bins
#' @param truth_gr GRanges with a numeric CN/signal column
#' @param signal_col name of the signal column in truth_gr
#' @return query_bins with an integrated `truth` column
rebin_truth <- function(query_bins, truth_gr, signal_col) {
  out <- scatools::integrate_segments(
    x = query_bins,
    y = truth_gr,
    granges_signal_colname = signal_col,
    drop_na = FALSE
  )
  out
}

#' Classify a CN signal into loss / neutral / gain.
#'
#' @param x numeric signal
#' @param mode "logratio" (thresholds around 0) or "integer" (around ploidy)
#' @param t threshold (logratio) — |x| < t is neutral
#' @param ploidy baseline for integer mode
classify_cn <- function(x, mode = c("logratio", "integer"), t = 0.15, ploidy = 2) {
  mode <- match.arg(mode)
  center <- if (mode == "logratio") 0 else ploidy
  tol <- if (mode == "logratio") t else 0.5
  out <- rep("neutral", length(x))
  out[x > center + tol] <- "gain"
  out[x < center - tol] <- "loss"
  factor(out, levels = c("loss", "neutral", "gain"))
}

#' Evaluate scatools CN against a truth profile on shared bins.
#'
#' @param query_gr GRanges from pseudobulk_profile() with a `signal` column (scatools)
#' @param truth_gr GRanges with a truth CN column
#' @param query_col signal column in query_gr
#' @param truth_col signal column in truth_gr
#' @param query_mode,truth_mode "logratio" or "integer" for classification
#' @return one-row data.frame of concordance metrics
evaluate_cnv <- function(query_gr, truth_gr,
                         query_col = "signal", truth_col = "cn",
                         query_mode = "logratio", truth_mode = "integer") {
  # Put truth onto the query bins.
  merged <- rebin_truth(query_gr, truth_gr, signal_col = truth_col)

  q <- S4Vectors::mcols(merged)[[query_col]]
  tr <- S4Vectors::mcols(merged)[[truth_col]]

  ok <- is.finite(q) & is.finite(tr)
  q <- q[ok]; tr <- tr[ok]
  if (length(q) < 10) stop("Too few overlapping bins to evaluate (", length(q), ")")

  # Continuous concordance
  pearson <- suppressWarnings(stats::cor(q, tr, method = "pearson"))
  spearman <- suppressWarnings(stats::cor(q, tr, method = "spearman"))

  # Categorical concordance (gain/neutral/loss)
  qc <- classify_cn(q, mode = query_mode)
  tc <- classify_cn(tr, mode = truth_mode)
  tab <- table(truth = tc, query = qc)
  accuracy <- sum(diag(tab)) / sum(tab)
  kappa <- cohen_kappa(tab)

  data.frame(
    n_bins = length(q),
    pearson = round(pearson, 3),
    spearman = round(spearman, 3),
    class_accuracy = round(accuracy, 3),
    kappa = round(kappa, 3)
  )
}

#' Cohen's kappa from a square confusion matrix.
cohen_kappa <- function(tab) {
  n <- sum(tab)
  po <- sum(diag(tab)) / n
  pe <- sum(rowSums(tab) * colSums(tab)) / n^2
  if (isTRUE(all.equal(pe, 1))) return(NA_real_)
  (po - pe) / (1 - pe)
}
