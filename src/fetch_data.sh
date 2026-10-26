#!/usr/bin/env bash
set -euo pipefail
DATA_DIR="data/raw"
mkdir -p "$DATA_DIR"

download_if_missing() {
    local url="$1"; local dest="$2"
    if [[ -f "$dest" ]]; then
        echo "[skip] $dest"
    else
        echo "[get ] $url"
        curl -L --fail --retry 3 -o "$dest" "$url"
    fi
}

download_if_missing \
  "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE5nnn/GSE5281/matrix/GSE5281_series_matrix.txt.gz" \
  "$DATA_DIR/GSE5281_series_matrix.txt.gz"

download_if_missing \
  "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE49nnn/GSE49036/matrix/GSE49036_series_matrix.txt.gz" \
  "$DATA_DIR/GSE49036_series_matrix.txt.gz"

echo "[ok ] data acquisition complete"
