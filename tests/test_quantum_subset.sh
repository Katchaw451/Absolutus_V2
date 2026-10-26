#!/usr/bin/env bash
set -euo pipefail
SELECTION="results/quantum_selection.json"
[[ -f "$SELECTION" ]] || { echo "[skip] no selection"; exit 0; }
python - <<'PY' "$SELECTION"
import json, sys
p = json.load(open(sys.argv[1]))
assert p["k"] > 0
assert len(p["selected_features"]) == p["k"]
print("PASS: selection file is well-formed")
PY
