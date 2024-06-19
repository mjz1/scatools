# devtools::load_all("~/work/repos/ascat/ASCAT/")


atac_ascpf <- function(...) {
  sce_bulk <- sce_clone
  # Parameters for eventual function
  normal_clust <- "N"
  logr_assay <- "counts_gc_modal_smoothed_ratios"
  baf_assay <- "mhf2"
  gender <- "XX"
  genomeVersion <- "hg38"


  ascat_listnames <- function() {
    c(
      "Tumor_LogR", "Tumor_BAF", "Tumor_LogR_segmented", "Tumor_BAF_segmented",
      "Germline_LogR", "Germline_BAF", "SNPpos", "ch", "chr", "chrs",
      "samples", "gender", "sexchromosomes", "X_nonPAR", "isTargetedSeq",
      "failedarrays"
    )
  }

  ascat_obj <- vector(mode = "list", length = length(ascat_listnames()))
  names(ascat_obj) <- ascat_listnames()

  tum_idx <- which(colnames(sce_bulk) != normal_clust)


  ### With normal as normal
  ascat_obj$Tumor_LogR <- as.data.frame(assay(sce_bulk[, tum_idx], logr_assay))
  ascat_obj$Tumor_BAF <- as.data.frame(assay(sce_bulk[, tum_idx], baf_assay))

  ascat_obj$Germline_LogR <- as.data.frame(assay(sce_bulk[, -tum_idx], logr_assay))
  ascat_obj$Germline_BAF <- as.data.frame(assay(sce_bulk[, -tum_idx], baf_assay))

  ## Or treat normal as tumor
  # ascat_obj$Tumor_LogR <- as.data.frame(assay(sce_bulk[,], logr_assay))
  # ascat_obj$Tumor_BAF <- as.data.frame(assay(sce_bulk[,], baf_assay))
  # ascat_obj$Germline_LogR <- as.data.frame(assay(sce_bulk[, -tum_idx], logr_assay))
  # ascat_obj$Germline_BAF <- as.data.frame(assay(sce_bulk[, -tum_idx], baf_assay))


  ascat_obj$SNPpos <- start(rowRanges(sce_bulk))

  # Split indices of each chromosome
  ascat_obj$ch <- split(seq_along(as.vector(seqnames(rowRanges(sce)))), as.vector(seqnames(rowRanges(sce))))
  ascat_obj$ch <- ascat_obj$ch[gtools::mixedsort(names(ascat_obj$ch))]

  # Split indices of each segmentable region (in this case chr arms)
  arms <- paste0(as.vector(seqnames(rowRanges(sce))), rowRanges(sce)$arm)
  ascat_obj$chr <- split(seq_along(arms), arms)
  ascat_obj$chr <- ascat_obj$chr[gtools::mixedsort(names(ascat_obj$chr))]

  # Vector of chromosome names
  ascat_obj$chrs <- unique(as.vector(seqnames(rowRanges(sce_bulk))))

  ascat_obj$gender <- gender

  if (!is.null(genomeVersion)) {
    if (genomeVersion == "hg19") {
      X_nonPAR <- c(2699521, 154931043)
    } else if (genomeVersion == "hg38") {
      X_nonPAR <- c(2781480, 155701382)
    } else {
      stop("genomeVersion must be either 'hg19' or 'hg38'.")
    }
  } else {
    X_nonPAR <- NULL
  }

  ascat_obj$X_nonPAR <- X_nonPAR

  ascat_obj$sexchromosomes <- c("X", "Y")

  ascat_obj$isTargetedSeq <- FALSE # Not sure what this one does yet

  ascat_obj$samples <- colnames(ascat_obj$Tumor_LogR)


  test <- ascat.asmultipcf(ascat_obj, penalty = 3, out.dir = NA, flip_baf = FALSE, winsorize = FALSE, min_bins = 3)

  # Pull the values out
  tum_logr_seg <- test$Tumor_LogR_segmented

  # tum_baf_seg <- do.call(cbind, test$Tumor_BAF_segmented)
  tum_baf_seg <- test$Tumor_BAF_segmented

  missing <- rownames(sce_bulk)[which(!rownames(sce_bulk) %in% rownames(tum_baf_seg))]

  sce_bulk_new <- sce_bulk[intersect(rownames(sce_bulk), rownames(tum_baf_seg)), colnames(test$Tumor_BAF_segmented)]

  assay(sce_bulk_new, "mhf_jointseg") <- assay(sce_bulk_new, "mhf2")

  assay(sce_bulk_new, "mhf_jointseg")[, tum_idx] <- tum_baf_seg
  # assay(sce_bulk_new, "mhf_jointseg")[, ] <- tum_baf_seg


  assay(sce_bulk_new, "ratios_jointseg")[, tum_idx] <- tum_logr_seg[rownames(sce_bulk_new), ]


  sce_bulk_new <- logNorm(sce_bulk_new, transform = "log2", assay_name = "ratios_jointseg", name = "logratios_jointseg")
  # assay(sce_bulk_new, "ratios_jointseg")[, ] <- tum_logr_seg[rownames(sce_bulk_new), ]
}
