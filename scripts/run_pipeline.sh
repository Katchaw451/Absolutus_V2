#!/usr/bin/env bash
set -euo pipefail

K="${1:-4}"
SOLVER="${2:-exact}"

echo "[1/5] Fetching...";     bash src/fetch_data.sh
echo "[2/5] Preprocessing..."; python src/preprocess.py
echo "[3/5] Network...";       python src/build_network.py --top-k 16
echo "[4/5] Optimizing (k=$K, solver=$SOLVER)..."
python src/quantum_optimizer.py --k "$K" --solver "$SOLVER"
echo "[5/5] Validating...";    bash scripts/verify_against_literature.sh results/quantum_selection.json
