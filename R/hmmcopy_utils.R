#' Single Cell HMMcopy
#'
#' Runs `HMMCopy` on single cell binned counts.
#'
#' @param chr Vector of chromosomes
#' @param start Vector of bin start positions
#' @param end Vector of bin end positions
#' @param counts Vector of bin corrected counts. These should ideally be GC corrected counts from [modal_quantile_regression()]
#' @param reads A vector of raw read counts per bin
#' @param param A matrix of parameter values in columns for each state in rows. See [HMMcopy::HMMsegment()] for more information.
#' @param cell_id Cell id
#' @param multiplier Ploidy multiplier
#' @param verbose Print verbose
#' @param maxiter The maximum number of iterations allows for the Maximum-Expectation algorithm, reduce to decrease running time at the expense of robustness.
#' @param n_cutoff Cutoff for the number of bins in a given state for calculating the `true_multiplier` value. Defaults to 5% of bins.
#'
#' @return A list with the following objects:
#' \describe{
#'   \item{bincounts}{Data frame of bin counts with copy state}
#'   \item{modal_seg}{Data frame of genome segments}
#'   \item{mstats}{A single row of cell summary statistics}
#'   \item{df_params}{Data frame of parameters used for each iteration and state}
#' }
#'
#' @export
#'
hmmcopy_singlecell <- function(chr, start, end, counts, reads, param = sc_hmm_params(), cell_id, multiplier = 1, verbose = FALSE, maxiter = 200, n_cutoff = NULL) {

  # In the original function, copy, cor_gc, and modal_corrected columns are all the same

  # Format the bincount data into a table
  bincounts <- data.frame(chr, start, end, reads, counts, multiplier)

  # Adjust the counts by the ploidy multiplier
  bincounts$copy <- bincounts$counts * multiplier

  # Run HMMcopy
  segmented <- HMMcopy::HMMsegment(bincounts, param = param, verbose = verbose, maxiter = maxiter)

  bincounts$state <- segmented$state - 1 # Why we subtract 1?

  # TWEAK
  # Trying to understand what this is doing -- we do an initial segmentation above, ... and adjust the copy values
  # Get median values per state
  meds <- bincounts %>%
    dplyr::group_by(state) %>%
    dplyr::summarize(median = median(copy, na.rm = TRUE), n = n()) %>%
    dplyr::mutate(fix = state / median)

  # Do we need to do this ordering
  # meds <- meds[order(meds$n, decreasing = TRUE), ]

  # Original cutoff was 200 for 500kb bins (of which there are 6073 which is about 3.3% of bins)
  if (is.null(n_cutoff)) {
    n_cutoff <- round(0.05 * sum(meds$n))
  }

  true_multiplier <- multiplier * mean(dplyr::filter(meds, n > n_cutoff)$fix, na.rm = TRUE)

  # Store the original copy value and states
  bincounts$copy_orig <- bincounts$copy
  bincounts$state_orig <- bincounts$state

  # Adjust the counts
  bincounts$copy <- bincounts$copy_orig * true_multiplier

  segmented_2 <- HMMcopy::HMMsegment(bincounts, param, verbose = verbose, maxiter = maxiter)

  # BASED 0 STATE
  bincounts$state <- segmented_2$state - 1

  modal_seg <- segmented_2$segs
  modal_seg$multiplier <- multiplier
  modal_seg$state <- as.numeric(as.character(modal_seg$state)) - 1

  # Sample stats
  stats <- modal_seg %>%
    dplyr::group_by(multiplier) %>%
    dplyr::summarise(MSRSI_non_integerness = median(abs(median - state), na.rm = TRUE))

  # add the seg means to the individual bins
  # TODO: Check for correct ordering
  rleseg <- rle(paste0(bincounts$chr, ":", bincounts$state))
  bincounts$median <- rep(modal_seg$median, rleseg$lengths)

  bincounts$halfiness <- -log2(abs(pmin(abs(bincounts$median - bincounts$state), 0.499) - 0.5)) - 1

  stats2 <- bincounts %>%
    dplyr::group_by(multiplier) %>%
    dplyr::summarise(
      MBRSI_dispersion_non_integerness = median(abs(copy - state), na.rm = TRUE),
      MBRSM_dispersion = median(abs(copy - median), na.rm = TRUE),
      autocorrelation_hmmcopy = tail(acf(copy_orig, 1, na.action = na.pass, type = "correlation", plot = FALSE)$acf, 1),
      cv_hmmcopy = sd(copy_orig, na.rm = TRUE) / mean(copy_orig, na.rm = TRUE),
      empty_bins_hmmcopy = sum(reads == 0, na.rm = TRUE),
      mad_hmmcopy = mad(copy_orig, constant = 1, na.rm = TRUE),
      mean_hmmcopy_reads_per_bin = mean(reads, na.rm = TRUE),
      median_hmmcopy_reads_per_bin = median(reads, na.rm = TRUE),
      std_hmmcopy_reads_per_bin = sd(reads, na.rm = TRUE),
      total_mapped_reads_hmmcopy = sum(reads, na.rm = TRUE),
      total_halfiness = sum(halfiness, na.rm = TRUE),
      scaled_halfiness = sum(halfiness / (state + 1), na.rm = TRUE)
    )

  stats3 <- bincounts %>%
    dplyr::group_by(state, multiplier) %>%
    dplyr::summarise(
      state_mads = mad(copy_orig, constant = 1, na.rm = TRUE),
      state_vars = var(copy, na.rm = TRUE)
    )

  stats4 <- stats3 %>%
    dplyr::group_by(multiplier) %>%
    dplyr::summarise(
      mean_state_mads = mean(state_mads, na.rm = TRUE),
      mean_state_vars = mean(state_vars, na.rm = TRUE)
    )

  mstats <- merge(merge(stats, stats2), stats4)
  neumad <- subset(stats3, state == 2)$state_mads
  mstats$mad_neutral_state <- ifelse(length(neumad) == 1, neumad, NA)

  mstats$breakpoints <- nrow(modal_seg) - length(unique(modal_seg$chr))
  mstats$mean_copy <- mean(bincounts$copy, na.rm = TRUE)
  mstats$state_mode <- as.numeric(names(tail(sort(table(bincounts$state)), 1)))
  mstats$log_likelihood <- tail(segmented_2$loglik, 1)
  mstats$true_multiplier <- true_multiplier

  # HAPLOID POISON
  ones <- bincounts$state == 1
  if (sum(ones) / length(ones) > 0.7) {
    mstats$scaled_halfiness <- Inf
  }

  df.params <- format_parameter_table(segmented_2, param)

  # add cellid
  df.params$cell_id <- cell_id
  bincounts$cell_id <- cell_id
  modal_seg$cell_id <- cell_id
  mstats$cell_id <- cell_id

  hmmresult <- list(bincounts = bincounts, modal_seg = modal_seg, mstats = mstats, df_params = df.params)

  return(hmmresult)
}



#' Single Cell HMMcopy
#'
#' A conveniece wrapper function for [hmmcopy_singlecell()] to test multiple candidate `multiplier` (aka ploidy) values and return either the best result, or a list of all results for downstream analysis.
#'
#' @inherit hmmcopy_singlecell
#'
#' @param multipliers Positive integer list of ploidy multipliers to test
#' @param mult_res a character. One of `best` or `all` to either return the result for the best ploidy only, or a list of results for all ploidies
#'
#' @export
#'
run_sc_hmmcopy <- function(chr, start, end, counts, reads, param = sc_hmm_params(), cell_id, multipliers = 1:6, verbose = FALSE, maxiter = 200, n_cutoff = NULL, mult_res = c("best", "all")) {
  # check integer multipliers
  if (!all(multipliers %% 1 == 0) | any(multipliers < 0)) {
    stop("Multipliers must be positive integers")
  }

  mult_res <- match.arg(mult_res)

  # Run the HMM for each multiplier state and compile the results into a list
  hmm_results <- lapply(multipliers, FUN = function(multiplier) {
    if (verbose) {
      message("Running single cell HMMcopy for multiplier: ", multiplier)
    }
    res <- hmmcopy_singlecell(
      chr = chr,
      start = start,
      end = end,
      counts = counts,
      reads = reads,
      param = param,
      cell_id = cell_id,
      multiplier = multiplier,
      verbose = verbose,
      maxiter = maxiter,
      n_cutoff = n_cutoff
    )
    return(res)
  })
  names(hmm_results) <- multipliers

  # Gather the summary statistics from each multiplier
  seg.best <- do.call("rbind", lapply(multipliers, FUN = function(multiplier) {
    hmm_results[[multiplier]]$mstats
  }))

  # scaledpenalty from original code is scaled_halfiness
  seg.best$red <- FALSE
  seg.best$red[which(seg.best$scaled_halfiness == min(seg.best$scaled_halfiness))] <- TRUE

  pick <- dplyr::filter(seg.best, red)$multiplier
  if (length(pick) > 1) {
    pick <- pick[1]
  }


  hmm_results["best"] <- pick

  message("Best ploidy: ", pick)

  if (mult_res == "best") {
    res <- list(
      bincounts = hmm_results[[pick]]$bincounts,
      modal_seg = hmm_results[[pick]]$modal_seg,
      mstats = hmm_results[[pick]]$mstats,
      df_params = hmm_results[[pick]]$df_params
    )
    return(res)
  }

  if (mult_res == "all") {
    return(hmm_results)
  }

  # auto_ploidy.reads <- hmm_results[[pick]]$bincounts
  # write.table(auto_ploidy.reads, sep = ",", quote = FALSE, row.names = FALSE, file = file.path(auto_output, "reads.csv"))

  # auto_ploidy.segs <- hmm_results[[pick]]$modal_seg
  # write.table(auto_ploidy.segs, sep = ",", quote = FALSE, row.names = FALSE, file = file.path(auto_output, "segs.csv"))

  # auto_ploidy.metrics <- hmm_results[[pick]]$mstats
  # write.table(auto_ploidy.metrics, sep = ",", quote = FALSE, row.names = FALSE, file = file.path(auto_output, "metrics.csv"))

  # auto_ploidy.params <- hmm_results[[pick]]$df_params
  # write.table(auto_ploidy.params, sep = ",", quote = FALSE, row.names = FALSE, file = file.path(auto_output, "params.csv"))
}


stack_params <- function(data, paramname) {
  data <- data.frame(data)
  colnames(data) <- 1:length(data) - 1
  data$state <- as.numeric(row.names(data)) - 1
  # TODO: Replace reshape2::melt with tidyr::pivot_longer here
  data <- reshape2::melt(data, id.vars = "state", value.name = "value", variable.name = "iteration")
  data$parameter <- paramname
  return(data)
}

format_parameter_table <- function(samp.segmented, new.params) {
  # mus - state medians
  # lambdas - state precision (inverse variance)
  # pi - state distribution
  # loglik  - likelihood values of each EM iteration

  num_iter <- ncol(samp.segmented$mus)

  loglik <- stack_params(t(samp.segmented$loglik), "loglik")
  loglik$state <- NaN

  nus <- stack_params(new.params$nu, "nus")
  nus$iteration <- NaN

  df.params <- rbind(
    stack_params(samp.segmented$mus, "mus"),
    stack_params(samp.segmented$lambdas, "lambdas"),
    stack_params(samp.segmented$pi, "pi"),
    loglik, nus
  )

  return(df.params)
}

#' Single Cell HMMcopy parameters
#'
#' Function to generate the parameter matrix for [HMMcopy::HMMsegment()]. Set-up with default values for single-cell analysis.
#'
#' @param e Probability of extending a segment, increase to lengthen segments, decrase to shorten segments. Range: (0, 1)
#' @param strength Strength of initial `e` suggestion, reducing allows e to change, increasing makes e undefinable. Range: [0, Inf)
#' @param mu Suggested median for copy numbers in state, change to readjust classification of states. Range: (-Inf, Inf)
#' @param lambda Suggested precision (inversed variance) for copy numbers in state, increase to reduce overlap between states. Range: [0, Inf)
#' @param nu Suggested degree of freedom between states, increase to reduce overlap between states. Range: [0, Inf)
#' @param kappa Suggested distribution of states. Should sum to 1. Must be same length as `mu`.
#' @param m Optimal value for mu, difference from corresponding mu value determines elasticity of the mu value. i.e. Set to identical value as mu if you don't want mu to move much.
#' @param eta Mobility of mu, increase to allow more movement. Range: [0, Inf)
#' @param gamma Prior shape on lambda, gamma distribution. Effects flexibility of lambda.
#' @param S Prior scale on lambda, gamma distribution. Effects flexibility of lambda.
#'
#' @inherit HMMcopy::HMMsegment details
#'
#' @return A data frame of parameters for [HMMcopy::HMMsegment()]
#' @export
#'
#' @examples
#' param <- sc_hmm_params()
sc_hmm_params <- function(e = (1 - 1e-6),
                          strength = 1000,
                          mu = 0:11,
                          lambda = 20,
                          nu = 2.1,
                          kappa = c(100, 100, 700, 100, 25, 25, 25, 25, 25, 25, 25, 25),
                          m = mu,
                          eta = 50000,
                          gamma = 3,
                          S = 1
                          ) {

  if (length(mu) != length(kappa)) {
    stop("kappa and mu must be the same length")
  }

  param <- data.frame(
    strength = strength,
    e = e,
    mu = mu,
    lambda = lambda,
    nu = nu,
    kappa = kappa,
    m = m,
    eta = eta,
    gamma = gamma,
    S = S
  )
  return(param)
}
