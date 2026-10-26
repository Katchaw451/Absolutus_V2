#!/usr/bin/env bash
set -euo pipefail
SCORE="${1:-shared_score}"
mkdir -p results/pareto
OUT="results/pareto/sweep_${SCORE}.csv"
printf "cap,objective,pathways,hubs,selected\n" > "$OUT"

for cap in 2 3 4 6 8; do
    if ! python src/quantum_optimizer.py --k 8 --solver exact \
            --max-per-pathway "$cap" --score-field "$SCORE" \
            --out "results/pareto/${SCORE}_cap_${cap}.json" >/dev/null 2>&1; then
        continue
    fi
    python - <<PY "results/pareto/${SCORE}_cap_${cap}.json" "$cap" >> "$OUT"
import json, sys
KNOWN = {"CALM3","CDC42","RAB3A"}
PW = {
  "mito":   {"NDUFA1","NDUFB2","SDHA","SDHB","COX5A","ATP5A1","ATP5B","UQCRC1","UQCRC2","CYCS"},
  "neuro":  {"IL6","TNF","IL1B","CXCL8","CCL2","IL10","IL4","TGFB1","NLRP3","TLR4"},
  "proteo": {"HSPA1A","HSPA8","PSMD1","UBB","SQSTM1","HSPB1","DNAJB1","BAG3","MAP1LC3B","ATG5"},
  "syn":    {"RAB3A","SNAP25","SYN1","CALM3","CDC42","STX1A","VAMP2","SNAP29","RAB5A","DNM1"},
}
d = json.load(open(sys.argv[1])); sel = set(d["selected_features"])
n_pw = sum(1 for m in PW.values() if m & sel)
print(f"{sys.argv[2]},{d['objective_value']:.4f},{n_pw},{len(KNOWN & sel)},{'|'.join(sorted(sel))}")
PY
done
echo "=== ${SCORE} ==="
column -t -s',' "$OUT"
