# scatools roadmap

Living plan for `scatools` — copy-number calling from single-cell ATAC-seq.
High-level here; discrete tasks are tracked as
[GitHub issues](https://github.com/mjz1/scatools/issues); datasets in
[`benchmarks/datasets.yaml`](benchmarks/datasets.yaml).

## Vision & positioning

Ship a **state-of-the-art, publishable, cleanly-coded** scATAC CNV method. The
defensible novelty claim (see competitive landscape below):

> The first **unified, reference-free** scATAC framework to jointly infer
> **absolute (ploidy-aware) total copy number and allele-specific copy number
> at single-cell resolution** — without matched bulk/scDNA sequencing —
> natively on `SingleCellExperiment`.

Load-bearing qualifiers: **reference-free** (vs Alleloscope, which needs matched
bulk DNA), **single-cell** (vs TeaCNV, subclone-level), **unified absolute +
allele-specific** (no competitor spans both).

## Status (done)

Green/fast CI, real test suite, ~13 fewer hard deps (Seurat dropped), runs
end-to-end on real normal (PBMC) and tumor (SNU601) scATAC. See merged PRs
#22–#31. Depth-CNV pipeline (`run_scatools`): bin → GC-correct (modal) → smooth
→ ratio → cluster (`cluster_sce`) → segment → merge → identify normal.

## Priorities (roughly ordered)

1. **Reference/population normalization** *(foundational)* — normal cells are
   not flat; per-cell depth ratios retain the systematic chromatin-accessibility
   profile, which masquerades as CNV. Needs normalization against a population
   median / normal reference before CNV calling is clean. Gates everything below.
2. **Validate against ground truth** — score depth calls vs matched scDNA
   (SNU601 `SNU601_dna.rds`) and vs published competitor results, in-notebook.
3. **Wire in absolute CN** — integrate the implemented-but-orphaned single-cell
   HMMcopy engine (`add_hmmcopy`, ploidy/multiplier search) into `run_scatools`
   as an absolute-CN mode. *Headline feature 1.*
4. **Finish the reference-free ASCN path** — `read_vartrix` → `phase_snps`
   (EM/population phasing) → BAF/LOH → allele-specific CN. Rewrite the `aspcf`
   scratch stub. *Headline feature 2.* (GH #8)
5. **Real mappability tracks** — currently a no-op placeholder. (GH #4)
6. **chr-arm-aware segmentation** (GH #6), **multipcf** (GH #7).
7. **Broaden benchmarks** — matched in-house DLP+ / scATAC (gold-standard
   truth) + more tumor datasets. NOTE (2026-07): turnkey `fragments.tsv.gz` for
   tumor scATAC is **rare** — the competitor papers (AtaCNV, epiAneufinder,
   Alleloscope) distribute **peak/bin matrices** (GEO) or raw BAM/FASTQ
   (SRA/dbGaP), not fragments. See `benchmarks/datasets.yaml` for the catalog
   (BCC GSE129785, pGBM GSE163655, PDAC GSE147726, OC/EC GSE173682). Turnkey
   fragments are basically only 10x public data (+ normal-only
   Human-scATAC-Corpus).
8. **Pre-binned / peak-matrix entry point** *(unlocks priority 7)* — add a
   `load_bin_matrix()`/`load_peak_matrix()` entry point (+ export
   `add_ideal_mat`/`calc_ratios`) so we can benchmark on the **same processed
   matrices the competitors used**, without needing fragments. Demonstrated
   ad-hoc for SNU601; make it first-class.
9. **Release** — CRAN vs Bioconductor decision (deferred); either way: minimize
   deps (done-ish), `--as-cran` clean, package size, no network in examples.

## Method-development findings / open questions

- **Normal ≠ flat (PBMC).** Depth ratios carry cell-type accessibility; needs
  normalization (priority 1).
- **Real data catches bugs.** Benchmark runs found 3 crashes hidden by the
  pre-subset test fixture: arm-gap binning, GC bin-boundary overshoot, and NA
  segmentation. Keep running real data early.
- **SNU601 not yet validated.** Pipeline runs (3001 cells → 4 subclones) but the
  segmented signal was modest; accuracy vs scDNA truth is unmeasured.

## Competitive landscape (who to beat)

- **epiAneufinder** (Nat Commun 2023) — relative 3-state; popular baseline.
- **CopyscAT**, **AtaCNV** (Cell Rep Methods 2025), **RIDDLER** — relative CN.
- **TeaCNV** (Brief Bioinform 2026) — clonal **absolute** CN (subclone-level).
- **Alleloscope** (Nat Biotechnol 2021), **Numbat-multiome** (2025) —
  **allele-specific** (Alleloscope needs matched bulk DNA).

## Known limitations

Depth normalization (above); mappability no-op; ASCN path unfinished; absolute
CN engine not wired in; `run_scatools` starts from fragments only (no pre-binned
entry point yet).

## Dev workflow

R via Singularity (`Rscript_ -v 4.3`); analysis/benchmarks in Quarto `.qmd`
notebooks under `benchmarks/notebooks/`; fast loop via `make_mini.sh`. See
[CLAUDE.md](CLAUDE.md).
