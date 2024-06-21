#' Bin atac fragments
#'
#' Convenience wrapper to bin fragments of a given size and save them as `.mtx` files
#'
#' @inherit run_scatools
#'
#' @param return_mat Logical. Return the binned depth matrix (default = FALSE)
#'
#' @return If `return_mat=TRUE`, returns a sparse binned depth matrix. Otherwise returns `NULL`, but will save the results to a 10X style market matrix `outdir`
#'
#' @export
bin_atac_frags <- function(sample_id,
                           fragment_file,
                           cells = NULL,
                           bins,
                           bin_name = NULL,
                           blacklist = NULL,
                           outdir,
                           ncores = 1,
                           bpparam = BiocParallel::SerialParam(),
                           overwrite = FALSE,
                           return_mat = FALSE) {

  stopifnot(class(bins) %in% "GRanges")

  if (is.null(bin_name)) {
    bin_name = prettyMb(getmode(width(bins)))
  }

  # Compute fragments per bins and combine
  bin_dir <- file.path(outdir, paste0(bin_name, "_counts"))
  dir.create(bin_dir, showWarnings = FALSE, recursive = TRUE)

  # Only bin frags if not done already
  if (!file.exists(file.path(bin_dir, "matrix.mtx.gz")) | overwrite) {

    # Convert to GenomicRanges
    logger::log_info("{sample_id} -- Loading fragments...")
    fragments <- data.table::fread(fragment_file, header = FALSE)
    logger::log_success("Fragments loaded!")
    colnames(fragments)[1:5] <- c("chr", "start", "end", "barcode", "umi")

    # If no cells provided use all barcodes and warn
    if (is.null(cells)) {
      cells = unique(fragments$barcode)
      if (length(cells) >= 1e5) {
       logger::log_warn("No cell barcodes specified. {length(cells)} unique cells in fragments file..." )
      }
    }

    fragments <- fragments[fragments$barcode %in% cells & fragments$chr %in% GenomeInfoDb::seqlevelsInUse(bins),]
    fragments <- GenomicRanges::makeGRangesFromDataFrame(fragments, keep.extra.columns = T)
    fragments <- GenomicRanges::split(fragments, GenomeInfoDb::seqnames(fragments))

    # Only keep fragments in called cells
    logger::log_info("{sample_id} -- Computing fragments in {bin_name} bins using {ncores} cores")
    count_mat <- BiocParallel::bplapply(
      X = fragments,
      FUN = bin_frags_chr,
      bins = bins,
      blacklist = blacklist,
      cells = cells,
      BPPARAM = BiocParallel::SerialParam()
    ) %>%
      do.call("rbind", .)



    if (!is.null(bin_dir)) {
      cat("Writing to", bin_dir, "\n")
      dir.create(bin_dir, showWarnings = FALSE, recursive = TRUE)

      # Write output using dropletutils
      DropletUtils::write10xCounts(
        path = bin_dir,
        x = count_mat,
        barcodes = colnames(count_mat),
        gene.id = rownames(count_mat),
        version = "3",
        overwrite = TRUE,
        gene.type = "Bin Counts"
      )
    }

    if (return_mat) {
      return(count_mat)
    }
  } else {
    logger::log_info("Binned counts file already found for {sample_id}")
  }
}


#' Bin scATAC fragments
#'
#' `bin_frags_chr` computes the fragments across bins in a single chromosome from an ArchR ArrowFile
#'
#' @param fragments_chr A single chromosome fragments
#'
#' @inherit run_scatools
#'
#' @return Sparse matrix of binned fragment counts
bin_frags_chr <- function(fragments_chr,
                          bins,
                          cells,
                          blacklist = NULL) {
  stopifnot(dplyr::n_distinct(as.vector(GenomeInfoDb::seqnames(fragments_chr))) == 1)
  stopifnot(class(bins) %in% "GRanges")
  # stopifnot(file.exists(fragment_file))

  chrom <- GenomeInfoDb::seqlevelsInUse(fragments_chr)

  logger::log_info("Binning {chrom}")

  # Remove fragments overlapping with blacklist regions
  if (!is.null(blacklist)) {
    if (!class(blacklist) %in% "GRanges") {
      logger::log_warn("Blacklist must be of class 'GRanges'. Provided: {class(blacklist)}")
    } else {
      blacklisted <- GenomicRanges::findOverlaps(subject = fragments_chr, query = blacklist)

      fragments_chr$blacklist <- FALSE
      mcols(fragments_chr)[S4Vectors::to(blacklisted), "blacklist"] <- TRUE

      fragments_chr <- fragments_chr[!fragments_chr$blacklist]
    }
  }


  # Chromosome bins
  bins_chr <- bins[as.vector(BSgenome::seqnames(bins)) %in% chrom]

  # Get overlapping indices
  # Note: Each fragment represents two transposition events
  # Therefore we count both the start and end site of each independently
  start_hits <- GenomicRanges::findOverlaps(
    subject = bins_chr,
    query = GRanges(
      seqnames = chrom,
      IRanges::IRanges(
        start = GenomicRanges::start(fragments_chr),
        end = GenomicRanges::start(fragments_chr)
      )
    )
  )
  end_hits <- GenomicRanges::findOverlaps(
    subject = bins_chr,
    query = GenomicRanges::GRanges(
      seqnames = chrom,
      IRanges::IRanges(
        start = GenomicRanges::end(fragments_chr),
        end = GenomicRanges::end(fragments_chr)
      )
    )
  )

  # Match Cells
  matchID <- S4Vectors::match(mcols(fragments_chr)$barcode, cells)

  # Create Sparse Matrix
  mat <- Matrix::sparseMatrix(
    i = c(S4Vectors::to(start_hits), S4Vectors::to(end_hits)),
    j = as.vector(c(matchID, matchID)),
    x = rep(1, 2 * length(fragments_chr)),
    dims = c(length(bins_chr), length(cells))
  )

  # Name matrix
  colnames(mat) <- cells
  rownames(mat) <- paste(GenomicRanges::seqnames(bins_chr),
                         GenomicRanges::start(bins_chr),
                         GenomicRanges::end(bins_chr),
                         sep = "_"
  )

  return(mat)
}

#' Get chromosome arm bins
#'
#' @param genome Genome version ('hg38', 'hg19')
#' @param calc_gc Logical: Whether or not to calculate GC content per bin
#' @param bs_genome BSgenome object. Must be passed if `calc_gc` is set to `TRUE`
#'
#' @return A GRanges object of chromosome arm bins
#' @export
#'
#' @examples
#' bins <- get_chr_arm_bins("hg38")
get_chr_arm_bins <- function(genome = "hg38", calc_gc = FALSE, bs_genome = NULL) {
  bins <- get_cytobands(genome = genome) %>%
    dplyr::group_by(CHROM, arm, genome) %>%
    dplyr::summarise(
      "start" = min(start),
      "end" = max(end)
    ) %>%
    dplyr::mutate("bin_id" = paste(CHROM, start, end, sep = "_")) %>%
    GenomicRanges::makeGRangesFromDataFrame(keep.extra.columns = T)
  bins$binwidth <- IRanges::width(bins)

  if (calc_gc) {
    if (is.null(bs_genome)) {
      stop("To calculate GC content you must pass a BSgenome object")
    }
    stopifnot(class(bs_genome) %in% "BSgenome")

    bins <- add_gc_freq(bins = bins, bs_genome = bs_genome)
  }
  bins <- sort(bins)

  bins$chr_arm <- paste0(GenomeInfoDb::seqnames(bins), bins$arm)

  return(bins)
}

#' Get tiled bins
#'
#' @param bs_genome A BSgenome object
#' @param tilewidth Bin size
#' @param select_chrs Vector of chromosomes to include
#' @param respect_chr_arms logical If `TRUE`, bins will be created with respect to chromosome arms (ie. not crossing arm boundries)
#'
#' @return A GRanges object of bins
#' @export
#'
#' @examples
#' \dontrun{
#' bins <- get_tiled_bins(BSgenome.Hsapiens.UCSC.hg38::BSgenome.Hsapiens.UCSC.hg38, tilewidth = 500000)
#' }
get_tiled_bins <- function(bs_genome = NULL,
                           tilewidth = 500000,
                           select_chrs = NULL,
                           respect_chr_arms = TRUE) {
  if (is.null(bs_genome)) {
    logger::log_error("Must provide 'bs_genome'")
    stop()
  }

  stopifnot(class(bs_genome) %in% "BSgenome")

  if (is.null(select_chrs)) {
    select_chrs <- paste("chr", c(1:22, "X"), sep = "")
  }
  logger::log_info("Creating chromosome bins for {unique(genome(bs_genome))[1]}")
  logger::log_info("respect_chr_arms = {respect_chr_arms}")
  logger::log_info("Chromosomes: {paste(select_chrs, collapse = ';')}")
  logger::log_info("Binwidth = {prettyMb(tilewidth)}")

  if (respect_chr_arms) {
    if (is.null(bs_genome)) {
      logger::log_error("Option: 'respect_chr_arms' requires specification of 'bs_genome'")
      stop()
    }
    # Get arm bins
    arm_bins <- get_chr_arm_bins(genome = unique(GenomeInfoDb::genome(bs_genome))[1])
    bins <- lapply(X = seq_along(arm_bins), FUN = function(i) {
      r <- arm_bins[i]

      starts <- seq(GenomicRanges::start(r), GenomicRanges::end(r), by = tilewidth)
      ends <- c(starts[2:(length(starts))] - 1, GenomicRanges::end(r))

      r_new <- GenomicRanges::GRanges(seqnames = seqnames(r), ranges = IRanges(start = starts, end = ends), arm = r$arm)
    }) %>%
      GenomicRanges::GRangesList() %>%
      unlist()
  } else {
    bins <- GenomicRanges::tileGenome(BSgenome::seqinfo(bs_genome)[select_chrs],
                                      tilewidth = tilewidth,
                                      cut.last.tile.in.chrom = TRUE
    )
  }

  bins$binwidth <- IRanges::width(bins)

  bins$bin_id <- paste(GenomeInfoDb::seqnames(bins), GenomicRanges::start(bins), GenomicRanges::end(bins), sep = "_")

  if (!is.null(select_chrs)) {
    bins <- GenomeInfoDb::keepSeqlevels(x = bins, value = select_chrs, pruning.mode = "coarse")
  }

  bins <- add_gc_freq(bs_genome, bins)
  bins <- sort(bins)
  idx <- which(as.vector(GenomeInfoDb::seqnames(bins)) %in% select_chrs)
  bins <- bins[idx,]

  return(bins)
}


#' Get genome cytobands
#'
#' @param genome Genome version (hg38, hg19, mm10)
#'
#' @return Dataframe of genome cytobands
#'
get_cytobands <- function(genome = "hg38") {
  cyto_url <- paste0("http://hgdownload.cse.ucsc.edu/goldenpath/", genome, "/database/cytoBand.txt.gz")
  cyto <- readr::read_delim(file = cyto_url, col_names = c("CHROM", "start", "end", "cytoband", "unsure"), show_col_types = FALSE) %>%
    dplyr::filter(!is.na(cytoband)) %>%
    dplyr::mutate(dplyr::across(where(is.character), as.factor),
                  "start" = start + 1,
                  "arm" = factor(substr(cytoband, 0, 1)),
                  "genome" = genome
    )
  return(cyto)
}


#' Add GC frequency
#'
#' @param bs_genome BSGenome object
#' @param bins GRanges bins object
#'
#' @return GRanges bin object with GC and N frequency per bin
add_gc_freq <- function(bs_genome, bins) {
  stopifnot(class(bs_genome) %in% "BSgenome")
  logger::log_info("Computing GC content for {prettyMb(getmode(width(bins)))} size bins")
  freqs <- BSgenome::alphabetFrequency(BSgenome::getSeq(bs_genome, bins))
  bins$gc <- (freqs[, "C"] + freqs[, "G"]) / rowSums(freqs)

  # Add N frequency
  bins$n_freq <- (freqs[, "N"]) / rowSums(freqs)

  return(bins)
}


#' Get ideal bin matrix
#'
#' Given a matrix of bin counts, bin gc and N frequency, and filtering parameters, return a boolean matrix flagging ideal bins
#'
#' @param mat,sce A count matrix or SCE object depending on the function
#' @param ncores number of cores for parallel evaluation (requires `pbmcapply` package)
#' @param assay_name Name of assay
#' @param verbose message verbosity
#'
#' @inherit is_ideal_bin
#' @return Boolean matrices of ideal and valid bins
get_ideal_mat <- function(mat, gc, n_freq, map, min_reads = 1, max_N_freq = 0.05, reads_outlier = 0.01, gc_outlier = 0.001, min_map = 0.9, ncores = 1, verbose = FALSE) {
  if (verbose) {
    logger::log_info("Computing ideal bins in {ncol(mat)} cells using {ncores} threads")
  }
  if (requireNamespace("pbmcapply")) {
    res <- do.call(
      "cbind",
      pbmcapply::pbmclapply(X = seq_len(ncol(mat)), mc.cores = ncores, FUN = function(i) {
        is_ideal_bin(
          counts = mat[, i],
          gc = gc,
          n_freq = n_freq,
          map = map,
          min_reads = min_reads,
          max_N_freq = max_N_freq,
          reads_outlier = reads_outlier,
          gc_outlier = gc_outlier,
          min_map = min_map
        )
      })
    )
  } else {
    logger::log_warn("No parallel backend detected. Ideal mat computation may be slow", call. = FALSE)
    res <- do.call(
      "cbind",
      lapply(X = seq_len(ncol(mat)), FUN = function(i) {
        is_ideal_bin(
          counts = mat[, i],
          gc = gc,
          n_freq = n_freq,
          map = map,
          min_reads = min_reads,
          max_N_freq = max_N_freq,
          reads_outlier = reads_outlier,
          gc_outlier = gc_outlier,
          min_map = min_map
        )
      })
    )
  }

  if (verbose) {
    logger::log_success("Computing ideal bins completed!")
  }

  # Sort of silly but works for now to return both matrices
  ideal_mat <- res[, grep("ideal", colnames(res))]
  colnames(ideal_mat) <- colnames(mat)
  rownames(ideal_mat) <- rownames(mat)
  valid_mat <- res[, grep("valid", colnames(res))]
  colnames(valid_mat) <- colnames(mat)
  rownames(valid_mat) <- rownames(mat)

  return(list(ideal = as.matrix(ideal_mat), valid = as.matrix(valid_mat)))
}


#' Length normalize counts
#'
#' By default this procedure will only impact counts in the bins that are variable length (for example tail ends of chromosomes).
#'
#' @param sce SCE object
#' @param assay_name Name of assay to normalize
#' @param assay_to Name of assay to save
#' @param binwidths Vector of binwidths
#' @param by_factor Multiplication factor for counts
#' @param verbose Print verbose (TRUE/FALSE)
#'
#' @return An sce object with counts length normalized in `assay(sce, 'counts_permb')`
length_normalize <- function(sce,
                             assay_name = "counts",
                             assay_to = "counts_lenNorm",
                             binwidths = width(rowRanges(sce)),
                             by_factor = getmode(binwidths),
                             verbose = FALSE) {
  if (verbose) {
    logger::log_info("Performing bin-length normalization. Storing as assay(sce, '", assay_to, "')")
  }

  # Bin length normalize upfront
  # Mainly necessary for our chr arm analysis where bins are variable in size
  # This is effectively an reads per Megabase (RPBMb) calculation. For reminder: https://www.rna-seqblog.com/rpkm-fpkm-and-tpm-clearly-explained/

  assay(sce, assay_to) <- assay(sce, assay_name) / binwidths * by_factor

  return(sce)
}



#' @rdname get_ideal_mat
#' @return SCE object with ideal and valid boolean matrices
add_ideal_mat <- function(sce,
                          assay_name = "counts",
                          gc = rowData(sce)$gc,
                          n_freq = rowData(sce)$n_freq,
                          map = rowData(sce)$map,
                          min_reads = 1,
                          max_N_freq = 0.05,
                          reads_outlier = 0.01,
                          gc_outlier = 0.001,
                          min_map = 0.9,
                          ncores = 1,
                          verbose = FALSE) {
  id_val_mats <- get_ideal_mat(
    mat = assay(sce, assay_name),
    gc = gc,
    n_freq = n_freq,
    map = map,
    min_reads = min_reads,
    max_N_freq = max_N_freq,
    reads_outlier = reads_outlier,
    gc_outlier = gc_outlier,
    min_map = min_map,
    ncores = ncores,
    verbose = verbose
  )

  assay(sce, "ideal_bins") <- id_val_mats$ideal
  assay(sce, "valid_bins") <- id_val_mats$valid

  return(sce)
}


#' Flag ideal bins
#'
#' `is_ideal_bin` will apply a set of bin-wise filters, based on high count outliers, high or low gc outliers, minimum read counts, minimum mappability, or maximum allowable frequency of N bases per bin.
#'
#' @param counts Vector of bin counts for single cell
#' @param gc Vector of gc content
#' @param n_freq Vector of bin N frequency (proportion of N bases in a bin)
#' @param map Vector of bin mappability
#' @param min_reads Minimum number of reads to consider a bin
#' @param max_N_freq Maximum allowable frequency of N bases to consider a bin. Range (0, 1)
#' @param reads_outlier Flag bins with reads in the top quantile given by this value. Range (0, 1)
#' @param gc_outlier Flag bins with GC content in the top and bottom quantule given by this value. Range (0, 1)
#' @param min_map Minimum allowable mappability score for a bin. Range (0, 1)
#'
#' @return A dataframe of two columns meet the `valid` or `ideal` criteria
is_ideal_bin <- function(counts, gc, n_freq, map = NULL, min_reads = 0, max_N_freq = 0.05, reads_outlier = 0.01, gc_outlier = 0.001, min_map = 0.9) {
  counts <- as.vector(counts)
  gc <- as.vector(gc)
  n_freq <- as.vector(n_freq)

  # Currently we use a placeholder for mappability. It doesn't do antyhing in the function yet
  if (is.null(map)) {
    map <- as.vector(rep(1, length(counts)))
  } else {
    map <- as.vector(map)
  }

  # Check lengths are equal and not zero
  if (var(unlist(lapply(list(counts, gc, n_freq, map), length))) != 0) {
    stop("counts, gc, n_freq, and map must be the same length. Check input data!")
  }

  if (length(counts) == 0) {
    stop("No data...")
  }

  # First identify valid bins
  valid <- is_valid_bin(counts = counts, n_freq = n_freq, min_reads = min_reads, max_N_freq = max_N_freq)

  # Subset for the valid counts/bins in computing quantiles
  # Remove high read outliers
  read_range <- stats::quantile(counts[valid], probs = c(0, 1 - reads_outlier))

  # Remove outlier GC bins on both sides
  gc_range <- stats::quantile(gc[valid], probs = c(gc_outlier, 1 - gc_outlier))

  # Is ideal if meeting all the following criteria
  ideal <- valid &
    (map > min_map) &
    (counts < read_range[2]) &
    (counts >= read_range[1]) &
    (gc < gc_range[2]) &
    (gc > gc_range[1])

  return(data.frame(ideal = ideal, valid = valid))
}


is_valid_bin <- function(counts, n_freq, min_reads = 0, max_N_freq = 0.05) {
  # Only valid if having both min reads and min n_freq
  counts > min_reads & n_freq <= max_N_freq
}

get_bin_ids <- function(granges) {
  # get the bin_ids
  bin_ids <- as.data.frame(granges) %>%
    dplyr::select(seqnames, start, end) %>%
    tidyr::unite("bin_id") %>%
    dplyr::pull()
  return(bin_ids)
}

#' @noRd
get_bin_info <- function(bin_ids) {
  # get the bin_ids
  bin_info <- data.frame(stringr::str_split_fixed(bin_ids, pattern = "_", n = 3))
  colnames(bin_info) <- c("chr", "start", "end")
  bin_info$chr <- factor(bin_info$chr, levels = gtools::mixedsort(unique(bin_info$chr)))
  bin_info$start <- as.numeric(bin_info$start)
  bin_info$end <- as.numeric(bin_info$end)
  return(bin_info)
}

#' Overlap genes with bins
#'
#' Given an annotation object, this function will overlap genes with bins and place the results in the metadata slot `gene_overlaps`. This function will also attempt to annotate cancer genes using OncoKB.
#'
#' @param sce SingleCellExperiment object
#' @param ensDb EnsemblDb object such as [EnsDb.Hsapiens.v86::EnsDb.Hsapiens.v86]
#' @param gene_biotype Gene biotypes to overlap. Defaults to "all". Use "protein_coding" to keep only protein coding genes
#'
#' @return SingleCellExperiment object with gene overlaps in `sce@metadata$gene_overlaps`
#' @export
#'
overlap_genes <- function(sce, ensDb, gene_biotype = "all") {
  # TODO: Allow selection of individual genes
  bin_ranges <- SummarizedExperiment::rowRanges(sce)

  # Pull genes
  g <- GenomicFeatures::genes(ensDb)
  g <- GenomeInfoDb::keepStandardChromosomes(g, pruning.mode = "coarse")

  # Keep the chromosome naming styles consistent prior to overlapping
  GenomeInfoDb::seqlevelsStyle(g) <- "UCSC"
  GenomeInfoDb::seqlevelsStyle(bin_ranges) <- "UCSC"

  # Subset g for chr levels in bin_ranges
  GenomeInfoDb::seqlevels(bin_ranges)
  g <- GenomeInfoDb::keepSeqlevels(x = g, value = GenomeInfoDb::seqlevels(bin_ranges), pruning.mode = "coarse")

  if (gene_biotype != "all") {
    logger::log_info("Filtering for gene_biotypes: {paste(gene_biotype, collapse = '; ')}")
    g <- g[which(mcols(g)$gene_biotype %in% gene_biotype)]
  }

  hits <- IRanges::findOverlaps(bin_ranges, g)

  # Want to attach this as a metadata granges to the sce where the new range for
  # each gene is the bin in which it resides

  mcols(g)["bin_id"] <- NA

  mcols(g[S4Vectors::subjectHits(hits)])["bin_id"] <- get_bin_ids(bin_ranges[S4Vectors::queryHits(hits)])

  # Label oncogenes and TS
  oncokb_df <- get_oncokb_genelist()

  # Match using both symbol and entrez for maximum overlap
  match_idx1 <- match(oncokb_df$hugo_symbol, g$gene_name)
  match_idx2 <- match(oncokb_df$entrez_gene_id, g$entrezid)
  # take from idx2
  match_idx1[is.na(match_idx1)] <- match_idx2[is.na(match_idx1)]

  oncokb_df$match_idx <- match_idx1
  # Log oncogenes which are missing from the match
  missing_g <- oncokb_df[is.na(oncokb_df$match_idx), "hugo_symbol"]
  if (length(missing_g > 0)) {
    logger::log_warn("Cancer genes missing from overlap: {paste(missing_g, collapse = '; ')}")
  }

  # Filter down to non-na to facilate merging
  oncokb_df <- oncokb_df[!is.na(oncokb_df$match_idx), ]

  # Perform the merge
  mcols(g)[oncokb_df$match_idx, colnames(oncokb_df)] <- oncokb_df

  sce@metadata$gene_overlap <- g

  return(sce)
}

#' @noRd
get_oncokb_genelist <- function(link = "https://www.oncokb.org/api/v1/utils/cancerGeneList.txt") {
  # Column names are problematic. Load these first
  colnames_cleaned <- read.table(
    file = link, header = F, sep = "\t",
    check.names = TRUE, nrows = 1, comment.char = ""
  )
  colnames_cleaned <- janitor::make_clean_names(colnames_cleaned[1, ])

  df <- read.table(
    file = link, header = FALSE, sep = "\t",
    col.names = colnames_cleaned, fill = F, skip = 1
  )

  return(df)
}

remove_zero_bins <- function(sce, assay_name = "counts", threshold = 0.85) {
  zeroes <- colSums(apply(assay(sce, assay_name), 1, function(x) x == 0))

  zeroes_prop <- zeroes / ncol(sce)

  keep <- zeroes_prop < threshold

  return(sce[names(which(keep)), ])
}

#' Pseudobulk groups
#'
#' @param sce sce obj
#' @param assay_name One or more assays to pseudobulk
#' @param group_var grouping variable
#' @param FUN funtion to use when pseudobulking
#' @param na.rm logical -- remove NAs in `FUN`
#' @param ... additional parameters to pass to `FUN`
#'
#' @return A pseudobulked sce obj
#' @export
#'
pseudo_groups <- function(sce, assay_name, group_var = NULL, FUN = mean, na.rm = TRUE, ...) {
  if (is.null(group_var)) {
    group_var <- "all"
    sce$all <- "all"
  }

  ids <- sce[[group_var]]

  by.group <- split(seq_along(ids), ids, drop = TRUE)

  res_list <- vector(mode = "list")

  for (a in assay_name) {
    res <- lapply(X = by.group, FUN = function(x) {
      apply(X = as.matrix(assay(sce, a)[, x]), MARGIN = 1, FUN = FUN, na.rm = na.rm, ...)
    })
    res_list[[a]] <- do.call("cbind", res)
  }

  res <- SingleCellExperiment(res_list)

  rowRanges(res) <- rowRanges(sce)
  rownames(res) <- rownames(sce)
  res$ncells <- unlist(lapply(by.group, length))
  res$ids <- colnames(res)
  res[[group_var]] <- colnames(res)

  return(res)
}
