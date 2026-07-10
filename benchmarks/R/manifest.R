# Read the benchmark dataset registry and resolve on-disk paths.
# Data lives outside git under SCATOOLS_BENCH_DATA.

#' Root directory for benchmark data (outside git).
bench_data_root <- function() {
  root <- Sys.getenv("SCATOOLS_BENCH_DATA", unset = "")
  if (!nzchar(root)) {
    root <- file.path(path.expand("~"), "work", "scatools_benchmark_data")
  }
  root
}

#' Read datasets.yaml into a list of dataset entries.
read_manifest <- function(path = file.path(dirname(dirname(this_dir())), "datasets.yaml")) {
  if (!requireNamespace("yaml", quietly = TRUE)) {
    stop("Package 'yaml' is required to read the benchmark manifest.")
  }
  yaml::read_yaml(path)$datasets
}

#' Resolve the on-disk directories for a dataset id.
dataset_paths <- function(id, create = FALSE) {
  base <- file.path(bench_data_root(), id)
  paths <- list(
    base = base,
    atac = file.path(base, "atac"),
    truth = file.path(base, "truth"),
    results = file.path(base, "results")
  )
  if (create) {
    lapply(paths, dir.create, showWarnings = FALSE, recursive = TRUE)
  }
  paths
}

#' Look up a single dataset entry by id.
get_dataset <- function(id, manifest = read_manifest()) {
  ids <- vapply(manifest, `[[`, character(1), "id")
  hit <- manifest[[which(ids == id)]]
  if (is.null(hit)) stop(sprintf("Dataset '%s' not found in manifest", id))
  hit
}

#' Tabular overview of the registry: id, aneuploidy class, truth type, status.
list_datasets <- function(manifest = read_manifest()) {
  do.call(rbind, lapply(manifest, function(d) {
    data.frame(
      id = d$id %||% NA,
      aneuploidy = d$aneuploidy %||% NA,
      assay = d$assay %||% NA,
      truth = (d$truth$type %||% NA),
      status = d$status %||% NA,
      stringsAsFactors = FALSE
    )
  }))
}

`%||%` <- function(a, b) if (is.null(a)) b else a

# Best-effort location of this script's directory (for default manifest path).
this_dir <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  f <- sub("^--file=", "", args[grep("^--file=", args)])
  if (length(f)) return(dirname(normalizePath(f)))
  # Fall back to the conventional repo location.
  "benchmarks/R"
}
