#!/usr/bin/env bash
set -euo pipefail

K="${1:-4}"
SOLVER="${2:-exact}"
CAP="${3:-}"

echo "[1/6] Fetching data...";        bash src/fetch_data.sh
echo "[2/6] Fetching annotation...";  bash src/fetch_annotation.sh
echo "[3/6] Preprocessing...";        python src/preprocess.py
echo "[4/6] Network (with hubs)...";  python src/build_network.py \
  --top-k 12 --per-pathway 2 --include-hubs CALM3 CDC42 RAB3A
echo "[5/6] Optimizing (k=$K, solver=$SOLVER, cap=${CAP:-none})..."
if [[ -n "$CAP" ]]; then
  python src/quantum_optimizer.py --k "$K" --solver "$SOLVER" --max-per-pathway "$CAP"
else
  python src/quantum_optimizer.py --k "$K" --solver "$SOLVER"
fi
echo "[6/6] Validating...";           bash scripts/verify_against_literature.sh results/quantum_selection.json
