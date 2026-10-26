#!/usr/bin/env bash
set -euo pipefail

K="${1:-4}"
SOLVER="${2:-exact}"

echo "[1/6] Fetching data...";        bash src/fetch_data.sh
echo "[2/6] Fetching annotation...";  bash src/fetch_annotation.sh
echo "[3/6] Preprocessing...";        python src/preprocess.py
echo "[4/6] Network...";              python src/build_network.py --top-k 10
echo "[5/6] Optimizing (k=$K, solver=$SOLVER)..."
python src/quantum_optimizer.py --k "$K" --solver "$SOLVER"
echo "[6/6] Validating...";           bash scripts/verify_against_literature.sh results/quantum_selection.json
