#' Get Gene Copies
#'
#' Calculate per group copy number for a list of genes.
#'
#' @param sce Single cell experiment object
#' @param assay_name Name of the assay from which to pull copy number
#' @param group_var Grouping variable (defaults to all)
#' @param gain_cutoff Relative copy increase to be considered a copy gain
#' @param loss_cutoff Relative copy decrease to be considered a copy loss
#' @param clonal_thr Proportion of clones harboring a gain or loss to be considered clonal
#' @param gene_names List of genes to get copy number
#' @param cancer_genes Logical. Whether to include results for cancer genes. Default=TRUE.
#'
#' @return A list with three elements:
#' \describe{
#'   \item{copy}{Wide form data frame with one row per  gene, and columns reflecting the average per clone copy, the relative copy increase or decrease versus clone average copy, clone average copy}
#'   \item{copy_long}{A melted version of `copy`}
#'   \item{clonal}{A filtered version of `copy_long` that classifies each mutation as private, shared, or clonal, and whether the type of gain or loss matches the gene's oncogene or tumor suppressor type.}
#' }
#' @export
#'
get_gene_copy <- function(sce,
                          assay_name,
                          group_var = "all",
                          gain_cutoff = 0.75,
                          loss_cutoff = -gain_cutoff,
                          clonal_thr = 0.75,
                          gene_names = NULL,
                          cancer_genes = TRUE) {
  # Check for the overlaps and create if possible
  if (is.null(sce@metadata$gene_overlap)) {
    # Attempt to perform the overlaps on the fly
    if (requireNamespace("EnsDb.Hsapiens.v86")) {
      logger::log_warn("No gene overlaps detected in SCE input. Performing overlaps now.")
      sce <- overlap_genes(sce = sce, ensDb = EnsDb.Hsapiens.v86::EnsDb.Hsapiens.v86, gene_biotype = "protein_coding")
    } else {
      logger::log_error("No gene overlaps detected in SCE input. Please run 'overlap_genes' prior to labelling genes.")
    }
  }


  bulk_all <- pseudobulk_sce(sce = sce, assay_name = assay_name, group_var = "all")

  bulk <- pseudobulk_sce(sce = sce, assay_name = assay_name, group_var = group_var)

  # Calculate average copy per sample and cluster so that we can highlight gains a losses from
  # uninteresting events
  mean_copies <- apply(assay(bulk_all), 2, mean)
  mean_copies_group <- apply(assay(bulk), 2, mean)
  # median_copies <- apply(assay(bulk), 2, median)


  if (!is.null(gene_names)) {
    gene_idx <- which(sce@metadata$gene_overlap$gene_name %in% gene_names &
      !is.na(sce@metadata$gene_overlap$bin_id))

    missing <- data.frame(cbind(gene_names, gene_names %in% sce@metadata$gene_overlap$gene_name)) %>%
      dplyr::filter(V2 == FALSE) %>%
      dplyr::pull(gene_names)

    if (length(missing) >= 1) {
      logger::log_warn("Unable to find genes: {paste(missing, collapse = '; ')}")
    } else {
      logger::log_info("Found all genes in 'gene_names'")
    }
  }

  if (cancer_genes == TRUE) {
    oncokb_idx <- which(sce@metadata$gene_overlap$onco_kb_annotated == "Yes" &
      !is.na(sce@metadata$gene_overlap$bin_id))

    gene_idx <- sort(union(gene_idx, oncokb_idx))
  }

  gene_ranges <- sce@metadata$gene_overlap[gene_idx]

  bulk_idx <- match(gene_ranges$bin_id, get_bin_ids(rowRanges(bulk)))

  if (group_var != "all") {
    new_colnames <- paste(assay_name, colnames(bulk), sep = "_")
    rel_colnames <- paste("relative_copies", colnames(bulk), sep = "_")
  } else {
    new_colnames <- paste(assay_name, group_var, sep = "_")
    rel_colnames <- paste("relative_copies", group_var, sep = "_")
  }

  mcols(gene_ranges)[new_colnames] <- assay(bulk[bulk_idx, ], assay_name)

  # Use sweep to subtract from each column
  mcols(gene_ranges)[rel_colnames] <- sweep(assay(bulk[bulk_idx, ], assay_name), 2, mean_copies_group, `-`)
  # mcols(onco_ranges)["relative_copies_all"] <- assay(bulk_all[bulk_idx, ], assay_name) - mean_copies

  df <- as.data.frame(gene_ranges)

  # df$sample_ploidy <- mean_copies

  # df[paste0("mean_ploidy_", names(mean_copies_group))] <- mean_copies_group


  df_long <- df %>%
    pivot_longer(
      cols = c(starts_with(assay_name), starts_with("relative_copies")),
      names_to = c(".value", "clone"),
      names_pattern = "(.+)_(.+)"
    ) %>%
    mutate(
      cn_cat = factor(
        ifelse(relative_copies >= gain_cutoff, "gain",
          ifelse(relative_copies <= loss_cutoff, "loss", "neutral")
        ),
        levels = c("gain", "neutral", "loss")
      )
    )

  n_clones <- length(unique(df_long$clone))


  clonal <- df_long %>%
    dplyr::filter(cn_cat != "neutral") %>%
    dplyr::group_by(gene_name, cn_cat) %>%
    dplyr::add_tally(name = "n") %>%
    dplyr::select(
      seqnames, bin_id, gene_name, gene_biotype, is_oncogene,
      is_tumor_suppressor_gene,
      clone, cn_cat, n, relative_copies, as.symbol(assay_name)
    ) %>%
    dplyr::mutate(
      clonal_prop = as.numeric(n / n_clones),
      clonality = factor(
        ifelse(clonal_prop >= clonal_thr, "clonal",
          ifelse(n == 1, "private", "shared")
        ),
        levels = c("private", "shared", "clonal")
      )
    )

  # Label instances of gains of oncogenes, and losses of tumor suppressors
  clonal$match_og <-
    with(clonal, factor(ifelse((cn_cat == "gain" & is_oncogene == "Yes") |
      (cn_cat == "loss" & is_tumor_suppressor_gene == "Yes"), "yes", "no")))

  return(list(copy = df, copy_long = df_long, clonal = clonal))
}
