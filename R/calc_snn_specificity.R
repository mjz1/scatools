#' Calculate SNN Specificity
#'
#' Calculates the fraction of nearest neighbours divided by the expected fraction of nearest neighbours in the patient subgraph.
#'
#' Also computes the neighbor purity estimation for each identity.
#'
#' @param snn_graph A SNN matrix
#' @param cell_idents1 Vector of cell identities for which the specificity will be calculated against. Ideally sample or patient.
#' @param cell_idents2 Vector of cell identities for which define cell types or states. The variable for which we want to see if the factor is specific to cell_idents1.
#' @param weighted Logical. Whether to weight the identities by their inverse size. Only relevent for the neighbor purity estimated
#' @param sampling_rate Proportion of cells to subsample when performing calculations.
#' @param ncores Number of cores to use
#'
#'
#' @return A dataframe
#' @export
#'
calc_snn_specificity <- function(snn_graph,
                                 cell_idents1,
                                 cell_idents2,
                                 sampling_rate = 0.1,
                                 weighted = T,
                                 ncores = 1) {
  # Refactored from: /work/shah/users/uhlitzf/projects/spectrum_explore/140_SPECTRUM_freeze_patient_mixing.Rmd

  ### For testing ###
  # snn_graph <- base::readRDS("/work/shah/users/uhlitzf/data/SPECTRUM/freeze/v5/cohort_snn_matrix.rds")
  #
  # meta_tbl <- readxl::read_excel("/work/shah/users/uhlitzf/projects/spectrum_explore/_data/small/MSK SPECTRUM - Single cell RNA-seq.xlsx", sheet = 3) %>%
  #   dplyr::mutate(pid1 = stringr::str_remove_all(patient_id, "SPECTRUM-OV-")) %>%
  #   dplyr::rename(sample = isabl_id)
  #
  # seu_tbl <- readr::read_tsv("/work/shah/isabl_data_lake/analyses/16/52/1652/cells.tsv") %>%
  #   dplyr::left_join(meta_tbl, by = "sample") %>%
  #   dplyr::filter(therapy == "pre-Rx") %>%
  #   dplyr::mutate(cell_type = stringr::str_replace_all(cell_type, "\\.", " "))
  # snn_graph <- snn_graph[seu_tbl$cell_id, seu_tbl$cell_id]
  # cell_idents1 <- stringr::str_sub(rownames(snn_graph), 13, 15)
  # cell_idents2 <- seu_tbl$cell_type
  #
  # ncores = 8
  # sampling_rate = 0.01
  ###############

  cell_idents1 <- as.factor(cell_idents1)
  cell_idents2 <- as.factor(cell_idents2)

  if (weighted == TRUE) {
    # Weight inversely by size of each identity
    c1_w <- 1 / table(cell_idents1)
    c1_w <- as.numeric(c1_w[as.character(cell_idents1)])

    c2_w <- 1 / table(cell_idents2)
    c2_w <- as.numeric(c2_w[as.character(cell_idents2)])
  } else {
    c1_w <- c2_w <- rep(1, length(cell_idents1))
  }


  # Compute the fraction of cell_idents2 from each cell_idents1
  cell_idents1_freq <- as.data.frame.matrix(prop.table(table(cell_idents1, cell_idents2), margin = 2)) %>%
    tibble::rownames_to_column(var = "cell_idents1") %>%
    tidyr::pivot_longer(cols = !cell_idents1, names_to = "cell_idents2", values_to = "nrel")

  # Subsample the graph according to the sampling rate
  withr::with_seed(3, code = {
    sample_idx <- sample.int(n = length(cell_idents1), size = floor(length(cell_idents1) * sampling_rate))
  })

  snn_graph_test <- snn_graph[, sample_idx]
  cell_idents1_test <- cell_idents1[sample_idx]
  cell_idents2_test <- cell_idents2[sample_idx]
  cell_ids_test <- colnames(snn_graph_test)

  # Loop over cells and compute the fraction of snn belonging to each cell_idents2
  # library(furrr)
  # prog_list <- list(
  #   type = "iterator",
  #   format = "Calculating {cli::pb_current}/{cli::pb_total} {cli::pb_bar} {cli::pb_percent} ETA: {cli::pb_eta}",
  #   clear = TRUE)

  # res <- furrr::future_map(cell_ids_test, .progress = T, .f = (cell_id) {
  res <- pbmcapply::pbmclapply(X = cell_ids_test, mc.preschedule = T, mc.cores = ncores, FUN = function(cell_id) {
    # Get neighbors of the current cell
    snn_idx <- which(snn_graph_test[, cell_id] != 0)

    snn <- snn_graph_test[snn_idx, cell_id]

    # Get the weights when summing the SNN
    c1_snn_w <- c1_w[snn_idx]
    c2_snn_w <- c2_w[snn_idx]

    # Identities of the neighbors for the cell being tested
    nb_cid1 <- cell_idents1[snn_idx]
    nb_cid2 <- cell_idents2[snn_idx]

    return(data.frame(snn, c1_snn_w, c2_snn_w, nb_cid1, nb_cid2))
  })

  res_df <- purrr::list_rbind(res, names_to = "idx") %>%
    dplyr::mutate(
      cell_id = cell_ids_test[idx],
      cid1 = cell_idents1_test[idx],
      cid2 = cell_idents2_test[idx], .after = "idx"
    ) %>%
    dplyr::left_join(cell_idents1_freq, by = c("cid1" = "cell_idents1", "cid2" = "cell_idents2")) %>%
    dplyr::mutate(
      cid1_match = (cid1 == nb_cid1),
      cid2_match = (cid2 == nb_cid2)
    ) %>%
    dplyr::group_by(cell_id) %>%
    dplyr::mutate(
      nexp = nrel * length(cell_id),
      nobs = sum(cid1 == nb_cid1),
      score = log2(nobs / (nexp + 1)),
      score2 = log2(nobs / nexp),
      cid1_w_total = sum(c1_snn_w),
      cid2_w_total = sum(c2_snn_w),
      cid1_w_target = sum(c1_snn_w * cid1_match),
      cid2_w_target = sum(c2_snn_w * cid2_match)
    ) %>%
    dplyr::distinct(
      cell_id,
      cid1,
      cid2,
      nexp,
      nobs,
      nrel,
      score,
      score2,
      cid1_w_total,
      cid1_w_target,
      cid2_w_total,
      cid2_w_target
    ) %>%
    dplyr::mutate(
      cid1_purity = cid1_w_target / cid1_w_total,
      cid2_purity = cid2_w_target / cid2_w_total
    )
  return(res_df)
}
