library(argparse)

# Parse command line arguments
parser <- ArgumentParser()

parser$add_argument("--sample_name", help="Sample name")
parser$add_argument("--fragment_file", help="Path to fragment file")
parser$add_argument("--genome_ver", help="Genome version (hg19, hg38) [default %(default)s]", default = "hg38")
parser$add_argument("--arrowdir", help="output directory for arrow file")
parser$add_argument("--binsize", help="Bin size for counts matrix [default %(default)s]", default = 1e7)
parser$add_argument("--bindepthdir", help="output directory for binned depth file")
parser$add_argument("--scatoolsdir", help="output directory for scatools processed files")
parser$add_argument("--reformat_frags", help="Logical: Reformat fragments before arrow [default %(default)s]", action="store_true")
parser$add_argument("--verbose", help="Message verbosity [default %(default)s]", action="store_true")
parser$add_argument("--overwrite", help="Logical: Overwrite existing files [default %(default)s]", action="store_true")


ncores = parallelly::availableCores()[[1]]

args <- parser$parse_args()

logger::log_info("Command line arguments:")
print(args)

sample_name = args$sample_name
fragment_file = args$fragment_file
genome_ver = args$genome_ver
arrowdir = args$arrowdir
binsize = args$binsize
bindepthdir = args$bindepthdir
scatoolsdir = args$scatoolsdir
reformat_frags = args$reformat_frags
verbose = args$verbose
overwrite = args$overwrite

logger::log_info("THREADS: {ncores}")


# devtools::load_all("~/work/repos/ArchR/")
library(ArchR)
devtools::load_all("~/repos/scatools")

# TESTING
# fragment_file = "/home/zatzmanm/work/repos/scatac_fetal_atlas/data/fragments/sample_59_bonemarrow.fragments.hg38.liftover.bed"
# sample_name = "sample_59_bonemarrow_TEST"
# arrowdir = "~/dump/create_arrows_test/ArrowFiles"
# bindepthdir = "~/dump/create_arrows_test/binned_depth"
# scatoolsdir = "~/dump/create_arrows_test/scatools_analysis"

# function from archr has bug to fix here
reformatFragmentFiles <- function(
    fragmentFiles = NULL,
    checkChrPrefix = getArchRChrPrefix()
){

  ArchR:::.validInput(input = fragmentFiles, name = "fragmentFiles", valid = c("character"))
  ArchR:::.validInput(input = checkChrPrefix, name = "checkChrPrefix", valid = c("boolean"))

  options(scipen = 999)
  ArchR:::.requirePackage("data.table")
  ArchR:::.requirePackage("Rsamtools")
  new_names <- vector(mode = "list")
  for(i in seq_along(fragmentFiles)){
    logger::log_info("{i} of {length(fragmentFiles)}: {fragmentFiles[i]}")
    fileNew <- gsub(".tsv.bgz|.tsv.gz|.bed.gz", "-Reformat.tsv", fragmentFiles[i])
    if (file.exists(paste0(fileNew, ".bgz"))) {
      logger::log_info("SKIPPING -- Reformatted fragments file already found: ", paste0(fileNew, ".bgz"))
      new_name <- paste0(fileNew, ".bgz")
      return(new_name)
      next
    }
    dt <- data.table::fread(fragmentFiles[i])
    dt <- dt[order(dt$V1,dt$V2,dt$V3), ]
    if(checkChrPrefix){
      idxRemove1 <- which(substr(dt$V1,1,3) != "chr")
    }else{
      idxRemove1 <- c()
    }
    idxRemove2 <- which(dt$V2 != as.integer(dt$V2))
    idxRemove3 <- which(dt$V3 != as.integer(dt$V3))
    #get all
    idxRemove <- unique(c(idxRemove1, idxRemove2, idxRemove3))
    if(length(idxRemove) > 0){
      dt <- dt[-idxRemove,]
    }
    if(nrow(dt) == 0){
      if(checkChrPrefix){
        stop("No fragments found after checking for integers and chrPrefix!")
      }else{
        stop("No fragments found after checking for integers!")
      }
    }
    #Make sure no spaces or #
    dt$V4 <- gsub(" |#", ".", dt$V4)
    # fileNew <- gsub(".tsv.bgz|.tsv.gz", "-Reformat.tsv", fragmentFiles[i])
    #### This is the update that fixes things
    ####
    data.table::fwrite(dt, fileNew, sep = "\t", col.names = FALSE)
    Rsamtools::bgzip(fileNew, overwrite = TRUE)
    file.remove(fileNew)
    new_name <- paste0(fileNew, ".bgz")
    # .fileRename(paste0(fileNew, ".bgz"), paste0(fileNew, ".gz"))
    new_names[[i]] <- new_name
  }
  return(new_names)
}

addArchRThreads(ncores)
addArchRGenome("hg38")

names(fragment_file) <- sample_name

if (reformat_frags == TRUE) {
  logger::log_info("REFORMATTING FRAGMENTS")
    new_names <- reformatFragmentFiles(
        fragmentFiles = fragment_file,
        checkChrPrefix = getArchRChrPrefix()
    )

    fragment_file <- unlist(new_names)
    names(fragment_file) <- sample_name
}

dir.create(arrowdir, recursive = TRUE, showWarnings = FALSE)
setwd(arrowdir)

logger::log_info("CREATING ARROW FILES")

ArrowFiles <- createArrowFiles(
  inputFiles = fragment_file,
  sampleNames = names(fragment_file),
  minTSS = 4, #Dont set this too high because you can always increase later
  minFrags = 1000,
  addTileMat = TRUE,
  addGeneScoreMat = TRUE,
  force = overwrite
)

logger::log_info("CALCULATING DOUBLETS")
doubScores <- addDoubletScores(
  input = ArrowFiles,
  k = 10, #Refers to how many cells near a "pseudo-doublet" to count.
  knnMethod = "LSI", #Refers to the embedding to use for nearest neighbor search.
  LSIMethod = 1,
  force = overwrite
)

proj <- ArchRProject(
  ArrowFiles = ArrowFiles,
  outputDirectory = arrowdir,
  copyArrows = FALSE #This is recommended so that you maintain an unaltered copy for later usage.
)

# Scatools processing
# Binned depth file output location
bins_out <- bindepthdir

dir.create(bins_out, recursive = TRUE, showWarnings = FALSE)

# Generate bins
bins <- get_tiled_bins(bs_genome = BSgenome.Hsapiens.UCSC.hg38::BSgenome.Hsapiens.UCSC.hg38, tilewidth = 1e7)

bin_name <- prettyMb(getmode(width(bins)))

# Bin the fragments
bin_atac_frags(ArrowFiles = getArrowFiles(proj), bins = bins, outdir = bins_out, ncores = ncores, overwrite = overwrite, return_mat = FALSE)


samples <- file.path(bins_out, bin_name, sample_name)

for (i in seq_along(samples)) {
  samp_dir <- samples[i]
  samp_name <- sample_name[i]

  samp_outdir <- file.path(scatoolsdir, bin_name, samp_name)

  raw_out <- file.path(samp_outdir, "sce", "01_raw.sce")
  final_out <- file.path(samp_outdir, "sce", "02_hmm.sce")
  hmm_out <- file.path(samp_outdir, "hmm", "hmm_results.rda")

  logger::log_info("Processing sample {i} of {length(samples)}: {samp_name}")
  
  if (file.exists(final_out) & !overwrite) {
    logger::log_info("Final output exists! Skipping to next sample...")
    next
  }

  if (file.exists(raw_out) & !overwrite) {
    logger::log_info("Raw sce object found -- Loading...")
    sce <- get(load(raw_out))
  } else {
    sce <- load_atac_bins(
      samples = samp_dir,
      sample.names = samp_name,
      ArchR_Proj = proj,
      bins = bins,
      BPPARAM = BiocParallel::bpparam(),
      save_to = raw_out,
      verbose = verbose
    )
  }
  
  sce <- sce %>%
    add_ideal_mat(ncores = ncores, verbose = verbose) %>%
    add_gc_cor(method = "modal", verbose = verbose, ncores = ncores) %>%
    add_hmmcopy(verbose = verbose, ncores = ncores, save_raw_hmm = hmm_out)
  
  save_to(object = sce, save_to = final_out)
}
