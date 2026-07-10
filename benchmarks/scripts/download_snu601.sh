#!/usr/bin/env bash
# Fetch SNU601 gastric scATAC processed matrices + matched scDNA truth from the
# Alleloscope repository (github.com/seasoncloud/Alleloscope, samples/SNU601).
# Processed data -> no cellranger needed. Raw FASTQ fallback: SRA PRJNA674903.
set -euo pipefail

ID=snu601
ROOT="${SCATOOLS_BENCH_DATA:-$HOME/work/scatools_benchmark_data}"
BASE="$ROOT/$ID"
ATAC="$BASE/atac"; TRUTH="$BASE/truth"
mkdir -p "$ATAC" "$TRUTH"

REPO="$BASE/Alleloscope"
if [[ ! -d "$REPO/.git" ]]; then
  echo "[$ID] sparse-cloning Alleloscope samples/SNU601 ..."
  git clone --depth 1 --filter=blob:none --sparse \
    https://github.com/seasoncloud/Alleloscope.git "$REPO"
  git -C "$REPO" sparse-checkout set samples/SNU601
else
  echo "[$ID] Alleloscope repo already present: $REPO"
fi

SRC="$REPO/samples/SNU601"
if [[ -d "$SRC" ]]; then
  # scATAC: bin-by-cell fragment counts + SNP ref/alt matrices, barcodes, VCF
  cp -rn "$SRC/scATAC/." "$ATAC/" 2>/dev/null || true
  # scDNA truth (matched single-cell genome; HMM CN segments)
  cp -rn "$SRC/scDNA/." "$TRUTH/" 2>/dev/null || true
  echo "[$ID] copied processed matrices:"
  ls -1 "$ATAC" | sed 's/^/    atac: /'
  ls -1 "$TRUTH" | sed 's/^/    truth: /'
else
  echo "[$ID] WARNING: $SRC not found — verify the sparse-checkout path in the repo."
  echo "[$ID] Raw FASTQ fallback: SRA PRJNA674903 (requires cellranger-atac)."
fi

echo "[$ID] done. Truth = matched scDNA (Andor 2020, PRJNA498809)."
