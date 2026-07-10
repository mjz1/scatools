# scatools benchmarks

Reproducible copy-number benchmarking for `scatools` against public and
in-house tumor scATAC datasets with orthogonal ground-truth copy number.

This directory holds **only** the registry, scripts, and evaluation harness.
The actual data (fragments, matrices, truth CN, results) lives **outside git**
under a data root, because the files are large.

## Design goals

- **Span the aneuploidy spectrum.** Datasets are tagged `normal` → `quiet` →
  `moderate` → `high` so we test both *specificity* (don't hallucinate CNVs in
  near-diploid genomes) and *sensitivity* (resolve complex aneuploidy).
- **Multiple ground-truth modalities.** Concordance is computed against whatever
  truth a dataset has: matched single-cell DNA (**DLP+**, the strongest — per-cell
  integer CN), bulk WGS/WES segments, multiome GEX-derived CN, SNP array, or a
  known karyotype.
- **One registry, many datasets.** Every dataset is a row in
  [`datasets.yaml`](datasets.yaml). Adding a dataset = adding an entry + a prep
  script; the harness is dataset-agnostic.

## Layout

```
benchmarks/
  README.md
  datasets.yaml          # the dataset registry (source of truth)
  R/
    manifest.R           # read the registry, resolve on-disk paths
    evaluate_cnv.R       # truth adapters + concordance metrics
  scripts/
    download_<id>.sh     # per-dataset fetch/prep (one per dataset)
    run_benchmark.R      # run scatools on a dataset -> results
```

Data root (NOT in git), set via `SCATOOLS_BENCH_DATA`
(default `~/work/scatools_benchmark_data`):

```
$SCATOOLS_BENCH_DATA/<dataset_id>/
  atac/        # fragments.tsv.gz / peak or bin matrix
  truth/       # ground-truth CN (DLP, WGS segments, ...)
  results/     # scatools output + evaluation
```

## Ground-truth modalities

| type          | source                              | resolution      |
|---------------|-------------------------------------|-----------------|
| `scDNA_DLP`   | Shah lab DLP+ (SIGNALS/HMMcopy)     | per-cell integer CN |
| `scWGS`       | single-cell WGS (non-DLP)           | per-cell / consensus |
| `bulk_wgs`    | matched bulk WGS/WES segments       | segment-level   |
| `multiome_gex`| paired GEX (inferCNV/Numbat)        | in-cell, coarse |
| `snp_array`   | SNP array CN                        | segment-level   |
| `karyotype`   | known/published karyotype           | arm-level       |

## Adding a dataset

1. Add an entry to `datasets.yaml` (id, aneuploidy class, atac source/format,
   truth type + accession).
2. Write `scripts/download_<id>.sh` that populates
   `$SCATOOLS_BENCH_DATA/<id>/atac/` and `.../truth/`.
3. Run the benchmark (`scripts/run_benchmark.R <id>`), which runs scatools and
   calls the evaluation harness.

## Status

Bootstrapping. See `datasets.yaml` for the current spectrum and per-dataset
`status` (`pending` → `downloaded` → `processed` → `benchmarked`).
