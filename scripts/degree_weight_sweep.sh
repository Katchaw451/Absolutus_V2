#!/usr/bin/env bash
set -euo pipefail
mkdir -p results/degree_weight
OUT="results/degree_weight/sweep.csv"
printf "weight,objective,pathways,CALM3,CDC42,RAB3A,selected\n" > "$OUT"

for w in 0 0.1 0.25 0.5 1.0 2.0; do
    JSON="results/degree_weight/w${w}.json"
    python src/quantum_optimizer.py \
        --k 8 --solver exact --max-per-pathway 4 \
        --degree-weight "$w" --out "$JSON" >/dev/null
    python - <<PY "$JSON" "$w" >> "$OUT"
import json, sys
KNOWN = {"CALM3","CDC42","RAB3A"}
PW = {
  "mito":   {"NDUFA1","NDUFB2","SDHA","SDHB","COX5A","ATP5A1","ATP5B","UQCRC1","UQCRC2","CYCS"},
  "neuro":  {"IL6","TNF","IL1B","CXCL8","CCL2","IL10","IL4","TGFB1","NLRP3","TLR4"},
  "proteo": {"HSPA1A","HSPA8","PSMD1","UBB","SQSTM1","HSPB1","DNAJB1","BAG3","MAP1LC3B","ATG5"},
  "syn":    {"RAB3A","SNAP25","SYN1","CALM3","CDC42","STX1A","VAMP2","SNAP29","RAB5A","DNM1"},
}
d = json.load(open(sys.argv[1])); sel = set(d["selected_features"]); w = sys.argv[2]
n_pw = sum(1 for m in PW.values() if m & sel)
calm  = int("CALM3" in sel)
cdc42 = int("CDC42" in sel)
rab3a = int("RAB3A" in sel)
print(f"{w},{d['objective_value']:.4f},{n_pw},{calm},{cdc42},{rab3a},{'|'.join(sorted(sel))}")
PY
done
column -t -s',' "$OUT"
