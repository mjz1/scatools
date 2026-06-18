# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## What this is

`scatools` is an R package (Bioconductor-style) for **copy number analysis in single-cell ATAC (scATAC) data**. It bins scATAC fragments across the genome, corrects technical biases, and infers copy number alterations (CNAs) per cell, with downstream clustering, segmentation, allele-specific analysis, and visualization.

The package is distributed via GitHub (`devtools::install_github("mjz1/scatools")`) and documented at https://mjz1.github.io/scatools/.

## Development commands

This is a standard R package using `devtools`/`roxygen2`. Run from an R session at the repo root:

```r
devtools::load_all()      # Load all functions for interactive development (Cmd/Ctrl+Shift+L in RStudio)
devtools::document()      # Regenerate man/*.Rd and NAMESPACE from roxygen2 comments — REQUIRED after changing any roxygen block
devtools::check()         # Full R CMD check (build + tests + examples + vignettes)
devtools::install()       # Install the package locally
devtools::build_vignettes()
pkgdown::build_site()     # Build the documentation website locally
```

Command line equivalents:

```bash
R CMD build .             # Build the source tarball
R CMD check scatools_*.tar.gz --no-manual
```

### Tests

There is **no `tests/` directory yet** — the package currently has no automated unit tests despite declaring `Config/testthat/edition: 3`. CI relies on `R CMD check` (examples + vignettes) for validation. If adding tests, create `tests/testthat/` and use `devtools::test()` (run a single file with `testthat::test_file("tests/testthat/test-foo.R")`).

### Dependencies / environment

- The project uses **`renv`** (`renv.lock`, `renv/`). `.Rprofile` sources `renv/activate.R`, so dependencies are restored into a project-local library. Run `renv::restore()` to install the locked versions; `renv::status()` to check drift.
- Many dependencies are **Bioconductor** packages (`SingleCellExperiment`, `GenomicRanges`, `DNAcopy`, `HMMcopy`, `ComplexHeatmap`, `rtracklayer`, etc.) plus large annotation packages (`EnsDb.Hsapiens.v86`, optional `BSgenome.Hsapiens.UCSC.hg38`). Installs are heavy.
- `Dockerfile` builds on `zatzmanm/rstudio` and installs the package via `devtools::install_github`.

## Architecture

### Central data structure

Everything revolves around a **`SingleCellExperiment` (SCE)** object:
- **rows = genomic bins** (a `GenomicRanges` `rowRanges`, e.g. 10Mb tiles or chromosome arms), **columns = cells** (barcodes).
- Per-bin metadata (GC content, blacklist overlap, `ideal` bin flags) lives in `rowData`; per-cell metadata (clusters, QC, `tumor_cell`/normal calls) lives in `colData`.
- The pipeline progressively **adds named assays** rather than mutating in place. The assay name encodes its processing lineage — this naming convention is load-bearing because functions take an `assay_name` input and write to a derived `assay_to`/`name` output:

  ```
  counts
    → counts_lenNorm           (length_normalize)
    → counts_gc_modal          (add_gc_cor, method = "modal")
    → counts_gc_modal_smoothed (smooth_counts)
    → ..._ratios               (calc_ratios)
    → logr_modal               (logNorm)
    → segment_merged_logratios (segment_cnv → merge_segments)
  ```

  When adding or wiring functions, preserve this `assay_name` in / `assay_to` out contract so steps remain composable via `%>%`.

### The main pipeline

`run_scatools()` (`R/run_scatools.R`) is the top-level convenience wrapper and the best map of the intended workflow. Given a `sample_id` and a `fragment_file`, it:

1. **`bin_atac_frags()`** — bins fragments into per-cell counts over `bins`, writing intermediate counts to `outdir` (uses `data.table` for efficient fragment loading).
2. **`load_atac_bins()`** — assembles the binned counts into an SCE.
3. **Filtering & correction chain** (piped): remove small leftover bins → `length_normalize` → `filter_sce` (GC range, min bin/cell counts & proportions) → `add_ideal_mat` → `add_gc_cor` → `smooth_counts` → `calc_ratios` → `logNorm` → `cluster_seurat`.
4. **Optional segmentation** (`segment = FALSE` by default): `segment_cnv` → `merge_segments` → `identify_normal`.
5. Saves `.sce` files and optionally `.h5ad` (requires `zellkonverter` + `anndata`).

Output files are written to `outdir` and named by bin size, e.g. `10Mb_raw.sce`, `10Mb_processed.sce`. If `<bin_name>_processed.sce` exists and `overwrite = FALSE`, the run short-circuits and reloads it.

### Code organization (`R/`)

Files are grouped by function, not by object. Key modules:

- **Binning / IO**: `bin_utils.R` (largest; binning, normalization, `get_tiled_bins`, `get_chr_arm_bins`, pseudobulk helpers), `load_atac_bins.R`, `rebin_sce.R`, `get_blacklist.R`.
- **CNV core**: `cnv_tools.R`, `gc_correction.R` (GC bias correction, modal method), `hmmcopy_utils.R` (HMMcopy segmentation), `integrate_segments.R`, `summarise_chr_arm.R`, `correct_atac_bias.R`.
- **Allele-specific**: `phase_snps.R`, `bin_snp_data.R`, `read_vartrix.R`, `aspcf.R`, `ascn.R`, `calc_allelic_imbalance.R`.
- **Clustering / QC**: `clustering.R` (`cluster_seurat`, `leiden_wrapper`), `quality_control.R`, `calc_snn_specificity.R`, `calc_clonal_diversity.R`, `identify_normal` logic.
- **Plotting**: `plotting.R` (largest plotting module; `cnaHeatmap`, `plot_cell_cna`, etc.), `colors.R`.
- **Interop**: `nb_to_sce.R` (numbat → SCE), `get_gene_copy.R`.
- **Infra**: `scatools-package.R` (`@import`/`@importFrom` roxygen for the whole package), `zzz.R` (`.onLoad`: sets `logger` color layout and ggplot `theme_classic`), `utils.R`, `utils-pipe.R` (`%>%` re-export), `data.R` (dataset docs).

The public API is whatever is `@export`ed (see `NAMESPACE`); roughly 60 of ~100 functions are exported. Internal helpers (e.g. `length_normalize`, `add_ideal_mat`, `calc_ratios`, `save_to`, `getmode`, `prettyMb`) are not exported but are used throughout the pipeline.

### Bundled data

- `data/bins_10mb.rda` — example hg38 10Mb `GenomicRanges` bins (regenerated by `data-raw/bins_10mb.R`).
- `data/test_sce.rda` — a small processed SCE (67 cells, 255 bins; from `data-raw/test_sce.R`) used in examples/vignettes.
- `inst/extdata/fragments.bed.gz` — example scATAC fragments (100 normal mammary cells) reachable via `system.file("extdata", "fragments.bed.gz", package = "scatools")`.

Scripts in `data-raw/` regenerate bundled datasets; rerun and `devtools::document()` after changing them.

## Conventions

- **Documentation is generated.** Edit roxygen2 `#'` comments above functions, then run `devtools::document()` — never hand-edit `man/*.Rd` or `NAMESPACE`. `README.md` is generated from `README.Rmd`; edit the `.Rmd`. The pkgdown site and `vignettes/TCGA.Rmd` (note the `.Rmd.orig` precomputation pattern) follow the same generated-artifact rule.
- **Logging**: use the `logger` package (`logger::log_info/log_success/log_error/log_warn`) with `glue`-style `{}` interpolation, matching existing code — not bare `message()`/`cat()`.
- **Parallelism**: prefer `BiocParallel` (`bpparam` argument, `BiocParallel::SerialParam()` default) over `mclapply`. The codebase is mid-migration toward this — newer functions (`add_gc_cor`, `segment_cnv`, `merge_segments`) take a `bpparam`; some still take `ncores`. Match the surrounding function's parameter.
- **Style**: 2-space indentation, `<-` for assignment, snake_case for functions and variables. Tidyverse (`dplyr`/`purrr`/`magrittr` pipe) and `data.table` are both used depending on the module.
- **Versioning**: bump `Version` in `DESCRIPTION` and add an entry to `NEWS.md` for user-facing changes. The `.9000` suffix denotes a development version.

## CI

GitHub Actions in `.github/workflows/`:
- `R-CMD-check.yaml` — `R CMD check` on macOS/Windows/Ubuntu across R release/oldrel (triggers on push/PR to `main`).
- `pkgdown.yaml` — builds and deploys the docs site.
- `pr-commands.yaml` — responds to `/document` and `/style` slash commands on PRs.
- `build_docker.yaml`, `rhub.yaml` — Docker image build and R-hub checks.
