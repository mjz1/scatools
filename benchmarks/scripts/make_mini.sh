#!/usr/bin/env bash
# Create a subsampled "mini" copy of a benchmark dataset for fast dev/testing.
#
# The mini keeps ALL fragments for a random subset of called cells (a *real*
# subset, NOT pre-filtered to the bins), so it still exercises the full binning
# path — including cut-sites in bin gaps — in seconds rather than minutes.
#
# Usage: make_mini.sh <dataset_id> [n_cells]   (default 300)
set -euo pipefail

ID="${1:?usage: make_mini.sh <dataset_id> [n_cells]}"
N="${2:-300}"
ROOT="${SCATOOLS_BENCH_DATA:-$HOME/work/scatools_benchmark_data}"
SRC="$ROOT/$ID/atac"
DST="$ROOT/${ID}_mini/atac"
mkdir -p "$DST"

if [[ ! -s "$SRC/cells.txt" ]]; then
  echo "error: $SRC/cells.txt not found (need called-cell barcodes)" >&2
  exit 1
fi

shuf -n "$N" "$SRC/cells.txt" > "$DST/cells.txt"
echo "[${ID}_mini] subsampled $(wc -l < "$DST/cells.txt") cells"

zcat "$SRC/fragments.tsv.gz" \
  | awk -F'\t' 'FNR==NR{c[$1];next} /^#/{print;next} ($4 in c)' "$DST/cells.txt" - \
  | gzip > "$DST/fragments.tsv.gz"

echo "[${ID}_mini] fragments: $(zcat "$DST/fragments.tsv.gz" | grep -vc '^#') rows, $(du -h "$DST/fragments.tsv.gz" | cut -f1)"
