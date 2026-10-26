"""build_network.py - overlap matrix with stratified pathway selection."""
from __future__ import annotations
import argparse, json
from pathlib import Path
import numpy as np, pandas as pd

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


def build_overlap(names):
    n = len(names)
    mat = np.zeros((n, n))
    memb = {x: {p for p, m in PATHWAYS.items() if x in m} for x in names}
    for i in range(n):
        for j in range(i + 1, n):
            s = memb[names[i]] & memb[names[j]]
            if s:
                mat[i, j] = mat[j, i] = float(len(s))
    return mat


def stratified_select(feats, top_k, per_pathway=2):
    """Take per_pathway genes from each pathway, then fill by priority."""
    selected, seen = [], set()
    for pathway, members in PATHWAYS.items():
        in_feats = [g for g in feats.index if g in members]
        if not in_feats:
            continue
        ranked = feats.loc[in_feats].sort_values("shared_score", ascending=False)
        for gene in ranked.index[:per_pathway]:
            if gene not in seen and len(selected) < top_k:
                selected.append(gene); seen.add(gene)
    # Fill remaining slots by global priority
    for gene in feats.sort_values("priority", ascending=False).index:
        if len(selected) >= top_k:
            break
        if gene not in seen:
            selected.append(gene); seen.add(gene)
    return selected


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--features", default="data/processed/features.csv")
    p.add_argument("--out", default="data/processed/overlap.npy")
    p.add_argument("--top-k", type=int, default=10)
    p.add_argument("--per-pathway", type=int, default=2,
                   help="Minimum genes per pathway before filling by priority.")
    args = p.parse_args()

    feats = pd.read_csv(args.features, index_col=0)
    selected = stratified_select(feats, args.top_k, args.per_pathway)
    top = feats.loc[selected]

    ov = build_overlap(selected)
    Path(args.out).parent.mkdir(parents=True, exist_ok=True)
    np.save(args.out, ov)
    Path(args.out).with_suffix(".names.json").write_text(json.dumps(selected))
    print(f"[ok ] overlap {ov.shape} -> {args.out}")
    print("Features in network (stratified):")
    for name in selected:
        row = feats.loc[name]
        flag = " [pathway]" if row.get("in_pathway", 0) else ""
        print(f"  - {name:<12} shared={row['shared_score']:.3f}{flag}")


if __name__ == "__main__":
    main()
