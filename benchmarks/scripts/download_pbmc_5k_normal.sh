#!/usr/bin/env bash
# Fetch the 10x PBMC 5k scATAC fragments (euploid negative control).
# ~961 MB. No cellranger needed.
set -euo pipefail

ID=pbmc_5k_normal
ROOT="${SCATOOLS_BENCH_DATA:-$HOME/work/scatools_benchmark_data}"
ATAC="$ROOT/$ID/atac"
mkdir -p "$ATAC"

BASE=https://cf.10xgenomics.com/samples/cell-atac/1.2.0/atac_pbmc_5k_nextgem
FRAG="$ATAC/fragments.tsv.gz"

if [[ -s "$FRAG" ]]; then
  echo "[$ID] fragments already present: $FRAG"
else
  echo "[$ID] downloading fragments (~961 MB) -> $FRAG"
  wget -q --show-progress -O "$FRAG" "$BASE/atac_pbmc_5k_nextgem_fragments.tsv.gz"
fi

# Tabix index (optional; scatools reads the fragments directly)
if [[ ! -s "$FRAG.tbi" ]]; then
  wget -q -O "$FRAG.tbi" "$BASE/atac_pbmc_5k_nextgem_fragments.tsv.gz.tbi" || \
    echo "[$ID] note: .tbi not fetched (not required)"
fi

echo "[$ID] done. Truth = diploid (no CNVs expected)."
