#!/usr/bin/env bash
set -euo pipefail
SELECTION="${1:-results/quantum_selection.json}"
MIN_HUBS="${2:-2}"
[[ -f "$SELECTION" ]] || { echo "[err] not found: $SELECTION" >&2; exit 2; }
python src/validate.py "$SELECTION" --min-hubs "$MIN_HUBS"
