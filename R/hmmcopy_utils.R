
#' Add HMMCopy results to sce object
#'
#' @param sce an SCE object
#' @param verbose Message verbosity
#' @param ncores Number of cores
#' @param save_raw_hmm Path to save raw hmm data in an `rda` file
#' @param slot_suffix Suffix to add to newly created `copy` and `state` assay slots.
#' @param assay_name Name of the assay with counts to input into HMMcopy. Ideally these are GC corrected.
#'
#' @return an sce object with hmm copy metadata added to coldata, and new slots `copy` and `state`
#' @export
#'
#' @examples
#' data(test_sce)
#' test_sce <- add_hmmcopy(test_sce)
add_hmmcopy <- function(sce,
                        verbose = FALSE,
                        ncores = 1,
                        assay_name = paste0("counts_gc_", sce@metadata$gc_cor_method),
                        save_raw_hmm = NULL,
                        slot_suffix = NULL) {
  if (verbose) {
    logger::log_info("Running HMMcopy on ", ncol(sce), " cells. Using ", ncores, " threads")
  }


  if (verbose) {
    logger::log_info("Input assay: {assay_name}")
  }


  # Get the matrices we need
  count_mat <- SummarizedExperiment::assay(sce, assay_name)
  reads_mat <- SummarizedExperiment::assay(sce, "raw_counts")
  ideal_mat <- SummarizedExperiment::assay(sce, "ideal_bins")

  chr <- as.factor(seqnames(rowRanges(sce)))
  start <- BiocGenerics::start(rowRanges(sce))
  end <- BiocGenerics::end(rowRanges(sce))

  hmm_results <- pbmcapply::pbmclapply(X = seq_len(ncol(count_mat)), mc.cores = ncores, FUN = function(i) {
    res <- suppressMessages(run_sc_hmmcopy(
      chr = chr,
      start = start,
      end = end,
      counts = count_mat[, i],
      reads = reads_mat[, i],
      ideal = ideal_mat[, i],
      cell_id = colnames(count_mat)[i],
      return = "all"
    ))
    return(res)
  })


  names(hmm_results) <- colnames(sce)

  if (verbose) {
    logger::log_success("HMMcopy completed!")
  }

  if (!is.null(save_raw_hmm)) {
    save_to(hmm_results, save_to = save_raw_hmm)
  }

  if (verbose) {
    logger::log_info("Grabbing best HMMcopy results")
  }
  hmm_results_best <- grab_hmm_res(hmm_results, ncores = 1)

  if (verbose) {
    logger::log_info("Adding HMMcopy metadata to sce")
  }
  hmm_metadata <- bind_sublist(hmm_results_best, "mstats")
  if (!all(hmm_metadata$cell_id == rownames(colData(sce)))) {
    logger::log_warn("Cell ids in HMMcopy metadata do not match the original sce object. Merge may be incomplete")
  }

  # Merge the metadata
  colData(sce) <- cbind(colData(sce), hmm_metadata[match(hmm_metadata$cell_id, rownames(colData(sce))), ])

  if (verbose) {
    logger::log_info("Adding HMMcopy data to sce")
  }

  copy_mat <- do.call("cbind", lapply(names(hmm_results_best), FUN = function(name) {
    dat <- hmm_results_best[[name]][["bincounts"]]$copy
  }))
  rownames(copy_mat) <- rownames(sce)
  colnames(copy_mat) <- names(hmm_results_best)

  state_mat <- do.call("cbind", lapply(names(hmm_results_best), FUN = function(name) {
    dat <- hmm_results_best[[name]][["bincounts"]]$state
  }))
  rownames(state_mat) <- rownames(sce)
  colnames(state_mat) <- names(hmm_results_best)

  if (!is.null(slot_suffix)) {
    copy_slot <- paste("copy", slot_suffix, sep = "_")
    state_slot <- paste("state", slot_suffix, sep = "_")
  } else {
    copy_slot <- "copy"
    state_slot <- "state"
  }

  if (verbose) {logger::log_info("Adding copy and state data as assays: {copy_slot} AND {state_slot}")}

  assay(sce, copy_slot) <- copy_mat
  assay(sce, state_slot) <- state_mat

  if (verbose) {
    logger::log_success("HMMcopy data added!")
  }

  # Store the modal segments as unstructured data

  return(sce)
}

#' Single Cell HMMcopy
#'
#' Runs `HMMCopy` on single cell binned counts.
#'
#' @param chr Vector of chromosomes
#' @param start Vector of bin start positions
#' @param end Vector of bin end positions
#' @param counts Vector of bin corrected counts. These should ideally be GC corrected counts from [gc_cor_modal()]
#' @param reads A vector of raw read counts per bin
#' @param ideal A logical vector indicating which bins are ideal for analysis. See [is_ideal_bin()]
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
hmmcopy_singlecell <- function(chr, start, end, counts, reads, ideal = rep(TRUE, length(counts)), param = params_sc_hmm(), cell_id, multiplier = 1, verbose = FALSE, maxiter = 200, n_cutoff = NULL) {

  # Format the bincount data into a table
  bincounts <- data.frame(
    chr = as.factor(chr),
    start = as.vector(start),
    end = as.vector(end),
    reads = as.vector(reads),
    counts = as.vector(counts),
    ideal = as.vector(ideal),
    multiplier
  )

  # Adjust the counts by the ploidy multiplier
  bincounts$copy <- bincounts$counts * multiplier

  # For the initial fit we do not include the ideal bins (this is to parallel the original codes intention)
  bincounts$copy[!bincounts$ideal] <- NA

  if (all(is.na(bincounts$copy))) {
    warning("No count data. Unable to segment.")
    hmmresult <- list(bincounts = bincounts, modal_seg = NA, mstats = .dummy_mstats(cell_id), df_params = NA)
    # TODO: See how this is handled in later merging
  } else {
    # Run HMMcopy
    segmented <- HMMcopy::HMMsegment(bincounts, param = param, verbose = verbose, maxiter = maxiter)

    bincounts$state <- segmented$state - 1 # Why we subtract 1?

    # TWEAK
    # Trying to understand what this is doing -- we do an initial segmentation above, ... and adjust the copy values
    # Get median values per state
    meds <- bincounts %>%
      dplyr::filter(ideal == TRUE) %>%
      dplyr::group_by(state) %>%
      dplyr::summarize(median = median(copy, na.rm = TRUE), n = dplyr::n()) %>%
      dplyr::mutate(fix = state / median)

    # Do we need to do this ordering
    # meds <- meds[order(meds$n, decreasing = TRUE), ]

    # Original cutoff was 200 for 500kb bins (of which there are 6073 which is about 3.3% of bins)
    if (is.null(n_cutoff)) {
      n_cutoff <- round(0.05 * sum(meds$n))
    }

    true_multiplier <- multiplier * mean(dplyr::filter(meds, n > n_cutoff)$fix, na.rm = TRUE)

    # Store the original copy value and states
    bincounts$copy_orig <- bincounts$counts # as with the original code, now reinclude the non-ideal bins in the refit
    bincounts$state_orig <- bincounts$state

    # Adjust the counts
    bincounts$copy <- bincounts$copy_orig * true_multiplier

    resegmented <- HMMcopy::HMMsegment(bincounts, param, verbose = verbose, maxiter = maxiter)

    # BASED 0 STATE
    bincounts$state <- resegmented$state - 1

    modal_seg <- resegmented$segs
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
      dplyr::filter(ideal == TRUE) %>%
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
      dplyr::filter(ideal == TRUE) %>%
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
    mstats$mean_copy <- mean(dplyr::filter(bincounts, ideal == TRUE)$copy, na.rm = TRUE)
    mstats$state_mode <- getmode(dplyr::filter(bincounts, ideal == TRUE)$state)
    mstats$log_likelihood <- tail(resegmented$loglik, 1)
    mstats$true_multiplier <- true_multiplier

    # HAPLOID POISON
    ones <- dplyr::filter(bincounts, ideal == TRUE)$state == 1
    if (sum(ones) / length(ones) > 0.7) {
      mstats$scaled_halfiness <- Inf
    }

    df.params <- format_parameter_table(resegmented, param)

    # add cellid
    df.params <- cbind(cell_id = as.factor(cell_id), df.params)
    bincounts <- cbind(cell_id = as.factor(cell_id), bincounts)
    modal_seg <- cbind(cell_id = as.factor(cell_id), modal_seg)
    mstats <- cbind(cell_id = as.factor(cell_id), mstats)

    # Memory savings here
    bincounts$chr <- as.factor(bincounts$chr)
    bincounts$state <- as.integer(bincounts$state)
    bincounts$state_orig <- as.integer(bincounts$state_orig)
    modal_seg$chr <- as.factor(modal_seg$chr)
    modal_seg$state <- as.integer(modal_seg$state)
    df.params$state <- as.integer(df.params$state)
    df.params$parameter <- as.factor(df.params$parameter)


    hmmresult <- list(bincounts = bincounts, modal_seg = modal_seg, mstats = mstats, df_params = df.params)
  }

  return(hmmresult)
}



#' Single Cell HMMcopy
#'
#' A conveniece wrapper function for [hmmcopy_singlecell()] to test multiple candidate `multiplier` (aka ploidy) values and return either the best result, or a list of all results for downstream analysis.
#'
#' @inherit hmmcopy_singlecell
#'
#' @param multipliers Positive integer list of ploidy multipliers to test
#' @param return a character. One of `best` or `all` to either return the result for the best ploidy only, or a list of results for all ploidies
#'
#' @export
#'
run_sc_hmmcopy <- function(chr, start, end, counts, reads, ideal = rep(TRUE, length(counts)), param = params_sc_hmm(), cell_id, multipliers = 1:6, verbose = FALSE, maxiter = 200, n_cutoff = NULL, return = c("best", "all")) {
  # check integer multipliers
  if (!all(multipliers %% 1 == 0) | any(multipliers < 0)) {
    stop("Multipliers must be positive integers")
  }

  return <- match.arg(return)

  # Run the HMM for each multiplier state and compile the results into a list
  hmm_results <- lapply(multipliers, FUN = function(multiplier) {
    if (verbose) {
      logger::log_info("Running single cell HMMcopy for multiplier: ", multiplier)
    }
    res <- hmmcopy_singlecell(
      chr = chr,
      start = start,
      end = end,
      counts = counts,
      reads = reads,
      ideal = ideal,
      param = param,
      cell_id = cell_id,
      multiplier = multiplier,
      verbose = verbose,
      maxiter = maxiter,
      n_cutoff = n_cutoff
    )
    return(res)
  })
  names(hmm_results) <- paste0("m", multipliers)

  # Gather the summary statistics from each multiplier
  seg.best <- do.call("rbind", lapply(names(hmm_results), FUN = function(multiplier) {
    hmm_results[[multiplier]]$mstats
  }))

  seg.best <- bind_sublist(hmm_results, sublist = "mstats")

  # Catch error cases when not enough ideal bins to segment
  if (all(is.na(seg.best$multiplier))) {
    # This will ensure mstats keeps track of failed cells
    pick <- "fail"
    hmm_results <- list("fail" = hmm_results[[1]]) # Just to reduce on space usage only need one
    hmm_results["best"] <- pick
    logger::log_info("Best ploidy: ", pick)
    pick_m <- pick
  } else {
    # scaledpenalty from original code is scaled_halfiness
    seg.best$red <- FALSE
    seg.best$red[which(seg.best$scaled_halfiness == min(seg.best$scaled_halfiness))] <- TRUE

    pick <- dplyr::filter(seg.best, red)$multiplier
    if (length(pick) > 1) {
      pick <- pick[1]
    }

    hmm_results["best"] <- pick

    logger::log_info("Best ploidy: ", pick)

    pick_m <- paste0("m", pick)
  }



  if (return == "best") {
    res <- list(
      bincounts = hmm_results[[pick_m]]$bincounts,
      modal_seg = hmm_results[[pick_m]]$modal_seg,
      mstats = hmm_results[[pick_m]]$mstats,
      df_params = hmm_results[[pick_m]]$df_params
    )
    return(res)
  }

  if (return == "all") {
    return(hmm_results)
  }
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

  df.params$parameter <- as.factor(df.params$parameter)

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
#' param <- params_sc_hmm()
params_sc_hmm <- function(e = (1 - 1e-6),
                          strength = 1000,
                          mu = 0:11,
                          lambda = 20,
                          nu = 2.1,
                          kappa = c(100, 100, 700, 100, 25, 25, 25, 25, 25, 25, 25, 25),
                          m = mu,
                          eta = 50000,
                          gamma = 3,
                          S = 1) {
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


#' Grab HMM results from returned list
#'
#' Helper function to parse large list of HMM results
#'
#' Note: Parallel currently disabled as it seems slower than single core
#'
#' @param hmm_results hmm results in multilist format
#' @param grab What to grab ("best")
#' @param ncores number of cores for parallelization
#'
#' @return A list of each cells best result\
#' @export
#'
grab_hmm_res <- function(hmm_results, grab = "best", ncores = 1) {
  if (requireNamespace("pbmcapply", quietly = TRUE) & TRUE == 1) {
    pbmcapply::pbmclapply(hmm_results, mc.cores = ncores, function(cell_hmm) {
      pick <- cell_hmm[[grab]]
      res <- cell_hmm[[pick]]
    })
  } else {
    lapply(hmm_results, function(cell_hmm) {
      pick <- cell_hmm[[grab]]
      res <- cell_hmm[[pick]]
    })
  }
}

# Should autogenerate based on test data within package on build...
.mstats_cols <- function() {
  c(
    "cell_id", "multiplier", "MSRSI_non_integerness", "MBRSI_dispersion_non_integerness",
    "MBRSM_dispersion", "autocorrelation_hmmcopy", "cv_hmmcopy",
    "empty_bins_hmmcopy", "mad_hmmcopy", "mean_hmmcopy_reads_per_bin",
    "median_hmmcopy_reads_per_bin", "std_hmmcopy_reads_per_bin",
    "total_mapped_reads_hmmcopy", "total_halfiness", "scaled_halfiness",
    "mean_state_mads", "mean_state_vars", "mad_neutral_state", "breakpoints",
    "mean_copy", "state_mode", "log_likelihood", "true_multiplier"
  )
}

.dummy_mstats <- function(cell_id) {
  x <- data.frame(t(rep(NA, length(.mstats_cols()))))
  colnames(x) <- .mstats_cols()
  x$cell_id <- cell_id
  return(x)
}
