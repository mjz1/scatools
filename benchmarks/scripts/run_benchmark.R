#!/usr/bin/env Rscript
# Run scatools on a benchmark dataset and evaluate against ground truth.
#
# Usage: Rscript benchmarks/scripts/run_benchmark.R <dataset_id> [assay_name]
#
# Expects the dataset's atac/ and truth/ dirs to be populated under
# SCATOOLS_BENCH_DATA (see benchmarks/scripts/download_<id>.sh).

suppressPackageStartupMessages({
  library(scatools)
})
source(file.path("benchmarks", "R", "manifest.R"))
source(file.path("benchmarks", "R", "evaluate_cnv.R"))

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) stop("Usage: run_benchmark.R <dataset_id> [assay_name]")
id <- args[[1]]
assay_name <- if (length(args) >= 2) args[[2]] else "logr_modal"

ds <- get_dataset(id)
paths <- dataset_paths(id, create = TRUE)
cli::cli_alert_info("Benchmarking dataset '{id}' ({ds$aneuploidy} aneuploidy)")

# 1) Run scatools on the ATAC fragments -----------------------------------
#    (Assumes atac/ holds a fragments.tsv.gz; adapt per dataset format.)
frag <- list.files(paths$atac, pattern = "fragments.*tsv.gz$", full.names = TRUE)
if (length(frag) != 1) {
  cli::cli_abort("Expected exactly one fragments file in {paths$atac}; found {length(frag)}")
}

bins <- get(data(bins_10mb))
sce <- run_scatools(
  sample_id = id,
  fragment_file = frag,
  bins = bins,
  outdir = paths$results,
  segment = TRUE,
  save_h5ad = FALSE
)

# 2) scatools pseudobulk / tumor profile ----------------------------------
query_gr <- pseudobulk_profile(sce[, sce$tumor_cell %in% TRUE], assay_name = assay_name)

# 3) Load ground truth (adapter per truth type) ---------------------------
#    TODO: implement load_truth_cn(ds, paths$truth) dispatching on ds$truth$type:
#      scDNA    -> per-cell integer CN -> clone consensus GRanges
#      bulk_wgs -> segment GRanges
#      multiome_gex -> paired-GEX-derived profile
truth_gr <- NULL
if (is.null(truth_gr)) {
  cli::cli_alert_warning("No truth adapter wired for '{ds$truth$type}' yet — skipping evaluation.")
  quit(save = "no")
}

# 4) Evaluate -------------------------------------------------------------
metrics <- evaluate_cnv(query_gr, truth_gr,
  query_col = "signal", truth_col = "cn",
  query_mode = "logratio", truth_mode = "integer"
)
metrics$dataset <- id
print(metrics)
saveRDS(metrics, file.path(paths$results, "metrics.rds"))
cli::cli_alert_success("Benchmark complete for '{id}'")
