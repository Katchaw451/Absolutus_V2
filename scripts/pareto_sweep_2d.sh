#!/usr/bin/env bash
set -euo pipefail
mkdir -p results/pareto_2d
printf "k,cap,objective,pathways,hubs,selected\n" > results/pareto_2d/sweep.csv

for k in 4 6 8 10; do
  for cap in 2 3 4; do
    OUT="results/pareto_2d/k${k}_cap${cap}.json"
    if ! python src/quantum_optimizer.py --k "$k" --solver exact \
            --max-per-pathway "$cap" --out "$OUT" >/dev/null 2>&1; then
      continue
    fi
    python - <<PY "$OUT" "$k" "$cap" >> results/pareto_2d/sweep.csv
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
n_hubs = len(KNOWN & sel)
print(f"{sys.argv[2]},{sys.argv[3]},{d['objective_value']:.4f},{n_pw},{n_hubs},{'|'.join(sorted(sel))}")
PY
  done
done
column -t -s',' results/pareto_2d/sweep.csv
