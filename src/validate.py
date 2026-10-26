"""validate.py - check recovered hubs and pathway coverage."""
from __future__ import annotations
import argparse, json
from pathlib import Path

KNOWN_HUBS = {"CALM3", "CDC42", "RAB3A"}
PATHWAYS = {
    "mitochondrial_dysfunction": {"NDUFA1","NDUFB2","SDHA","SDHB","COX5A",
                                  "ATP5A1","ATP5B","UQCRC1","UQCRC2","CYCS"},
    "neuroinflammation": {"IL6","TNF","IL1B","CXCL8","CCL2","IL10","IL4",
                          "TGFB1","NLRP3","TLR4"},
    "proteostasis": {"HSPA1A","HSPA8","PSMD1","UBB","SQSTM1","HSPB1",
                     "DNAJB1","BAG3","MAP1LC3B","ATG5"},
    "synaptic_function": {"RAB3A","SNAP25","SYN1","CALM3","CDC42",
                          "STX1A","VAMP2","SNAP29","RAB5A","DNM1"},
}


def main():
    p = argparse.ArgumentParser()
    p.add_argument("selection")
    p.add_argument("--min-hubs", type=int, default=2)
    args = p.parse_args()

    sel = set(json.loads(Path(args.selection).read_text())["selected_features"])
    recovered = KNOWN_HUBS & sel

    print(f"Selected ({len(sel)}): {sorted(sel)}")
    print(f"Recovered hubs:        {sorted(recovered)}")
    print(f"Missing hubs:          {sorted(KNOWN_HUBS - sel)}")
    print()
    print("Pathway coverage:")
    covered = 0
    for pathway, members in PATHWAYS.items():
        hit = members & sel
        if hit:
            covered += 1
            print(f"  [x] {pathway}: {sorted(hit)}")
        else:
            print(f"  [ ] {pathway}")
    print()
    print(f"Pathways covered: {covered}/{len(PATHWAYS)}")

    if len(recovered) >= args.min_hubs:
        print(f"PASS: recovered >= {args.min_hubs} known hub(s).")
        raise SystemExit(0)
    if covered >= 3 and len(recovered) >= 1:
        print(f"PASS (weak): {len(recovered)} hub(s) + {covered} pathways. "
              f"Interpretable module but not hub-driven.")
        raise SystemExit(0)
    print("FAIL: insufficient coverage.")
    raise SystemExit(1)


if __name__ == "__main__":
    main()
