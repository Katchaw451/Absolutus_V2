#!/usr/bin/env bash
set -euo pipefail

ANN_DIR="data/annotation"
mkdir -p "$ANN_DIR"

fetch_gpl() {
    local gpl="$1"
    local dest="$ANN_DIR/${gpl}.annot.gz"
    if [[ -f "$dest" ]]; then
        echo "[skip] $dest"
        return 0
    fi
    # GEO FTP pattern for GPL570 is GPLnnn (not GPL5nnn)
    local family="GPLnnn"
    local url="https://ftp.ncbi.nlm.nih.gov/geo/platforms/${family}/${gpl}/annot/${gpl}.annot.gz"
    echo "[get ] $url"
    curl -L --fail --retry 3 -o "$dest" "$url"
}

# AD platform: GSE5281 uses GPL570
fetch_gpl "GPL570"

# Detect PD platform from GSE49036 metadata
PD_PLATFORM=$(curl -s "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE49036&targ=self&form=text&view=brief" \
              | grep -oP '(?<=!Series_platform_id = )\S+' | head -1 || true)
echo "[info] PD platform detected: ${PD_PLATFORM:-unknown}"
if [[ -n "${PD_PLATFORM:-}" ]] && [[ "$PD_PLATFORM" != "GPL570" ]]; then
    fetch_gpl "$PD_PLATFORM" || echo "[warn] PD annotation unavailable; will fall back to AD mapping"
fi

echo "[ok ] annotation complete"
