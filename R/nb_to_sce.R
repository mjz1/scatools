#' @export
get_counts1 <- function(df_allele, column) {
  # This is less memory efficient but possible faster original version
  counts <- df_allele %>%
    tidyr::pivot_wider(id_cols = c("cell"), names_from = "snp_id", values_from = column) %>%
    tibble::column_to_rownames("cell") %>%
    as.matrix()
  counts[is.na(counts)] <- 0

  counts <- as(t(counts), "sparseMatrix")
}

#' @export
get_counts <- function(df_allele, column, ncores = 1) {
  requireNamespace("data.table")

  # Much more memory efficient
  dt_split <- split(data.table::as.data.table(df_allele), df_allele$CHROM)

  if (requireNamespace("pbmcapply", quietly = TRUE) & ncores > 1) {
    message("Faster snps")
    res <- pbmcapply::pbmclapply(seq_along(dt_split), mc.cores = ncores, FUN = function(i) {
      message(i)
      dt_chr <- dt_split[[i]]
      dt_chr_cast <- data.table::dcast(dt_chr, snp_id ~ cell, value.var = column, fill = 0, drop = FALSE)

      dt_res <- mltools::sparsify(dt_chr_cast, sparsifyNAs = TRUE)
      rownames(dt_res) <- dt_chr_cast$snp_id

      rm(dt_chr_cast)
      return(dt_res)
    })
  } else {
    message("slower snps")

    res <- lapply(seq_along(dt_split), FUN = function(i) {
      message(i)
      dt_chr <- dt_split[[i]]
      dt_chr_cast <- dcast(dt_chr, snp_id ~ cell, value.var = column, fill = 0, drop = FALSE)

      dt_res <- mltools::sparsify(dt_chr_cast, sparsifyNAs = TRUE)
      rownames(dt_res) <- dt_chr_cast$snp_id

      rm(dt_chr_cast)

      return(dt_res)
    })
  }

  counts <- merge.sparse(l = res)

  return(counts)
}

#' @export
numbat_to_sce <- function(df_allele, ncores = 1, min_snp_cov = 1) {
  alt <- get_counts(df_allele, column = "AD", ncores = ncores)
  tot <- get_counts(df_allele, column = "DP", ncores = ncores)
  ref <- tot - alt

  # keep_snps <- names(which(rowSums(tot) >= min_snp_cov))

  snp_ranges <- as.data.frame(str_split(rownames(alt), pattern = "_", simplify = TRUE))
  colnames(snp_ranges) <- c("chr", "start", "ref", "alt")
  snp_ranges$snp_id <- rownames(alt)

  # snp_metadata <- distinct(df_allele[, c("CHROM", "POS", "REF", "ALT", "GT", "gene", "sample", "group")]) %>%
  # mutate(snp_id = paste(CHROM, POS, REF, ALT, sep = "_"))

  cell_metadata <- distinct(df_allele[, c("cell", "sample", "group")])
  rownames(cell_metadata) <- cell_metadata$cell

  cell_metadata <- cell_metadata[colnames(alt), ]

  snp_ranges <- GenomicRanges::makeGRangesFromDataFrame(snp_ranges, end.field = "start", keep.extra.columns = TRUE) %>%
    sort()

  alt <- alt[snp_ranges$snp_id, ]
  ref <- ref[snp_ranges$snp_id, ]

  sce <- SingleCellExperiment::SingleCellExperiment(list(ref = ref, alt = alt), rowRanges = snp_ranges, colData = cell_metadata)

  return(sce)
}


#' Merge sparse matrices
#'
#' Can merge when rows are non matching
#'
#' @param l List of sparse matrices
#'
#' @return A merged matrix
#' @export
#'
merge.sparse <- function(l) {
  # https://stackoverflow.com/questions/43117608/r-binding-sparse-matrices-of-different-sizes-on-rows

  require("Matrix")

  cnnew <- character()
  rnnew <- character()
  x <- vector()
  i <- numeric()
  j <- numeric()

  for (M in l) {
    cnold <- colnames(M)
    rnold <- rownames(M)

    cnnew <- union(cnnew, cnold)
    rnnew <- union(rnnew, rnold)

    cindnew <- match(cnold, cnnew)
    rindnew <- match(rnold, rnnew)
    ind <- Matrix::summary(M)
    i <- c(i, rindnew[ind[, 1]])
    j <- c(j, cindnew[ind[, 2]])
    x <- c(x, ind[, 3])
  }

  m <- Matrix::sparseMatrix(i = i, j = j, x = x, dims = c(length(rnnew), length(cnnew)), dimnames = list(rnnew, cnnew))
  return(m)
}
