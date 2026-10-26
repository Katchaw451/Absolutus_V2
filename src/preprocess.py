"""preprocess.py - map probes to gene symbols, filter confounders, score."""
from __future__ import annotations
import argparse, gzip
from pathlib import Path
import numpy as np, pandas as pd

# Sex-linked genes that confound case/control comparisons
SEX_LINKED = {
    "XIST", "RPS4Y1", "RPS4Y2", "EIF1AY", "USP9Y", "KDM5D", "DDX3Y",
    "UTY", "ZFY", "TXLNGY", "NLGN4Y", "TMSB4Y", "PRKY", "AMELY",
    "CYorf15A", "CYorf15B",
}

# Curated pathway membership for a priori prioritisation
PATHWAY_MEMBERS = {
    "mitochondrial_dysfunction": {"NDUFA1","NDUFB2","SDHA","SDHB","COX5A",
                                  "ATP5A1","ATP5B","UQCRC1","UQCRC2","CYCS"},
    "neuroinflammation": {"IL6","TNF","IL1B","CXCL8","CCL2","IL10","IL4",
                          "TGFB1","NLRP3","TLR4"},
    "proteostasis": {"HSPA1A","HSPA8","PSMD1","UBB","SQSTM1","HSPB1",
                     "DNAJB1","BAG3","MAP1LC3B","ATG5"},
    "synaptic_function": {"RAB3A","SNAP25","SYN1","CALM3","CDC42",
                          "STX1A","VAMP2","SNAP29","RAB5A","DNM1"},
}


def _read_series_matrix(path):
    rows, in_table = [], False
    with gzip.open(path, "rt", encoding="utf-8", errors="ignore") as fh:
        for line in fh:
            if line.startswith("!series_matrix_table_begin"):
                in_table = True; continue
            if line.startswith("!series_matrix_table_end"):
                break
            if in_table:
                rows.append(line.rstrip("\n").split("\t"))
    header = [h.strip().strip('"') for h in rows[0]]
    df = pd.DataFrame(rows[1:], columns=header).set_index(header[0])
    df.index = df.index.str.strip('"')
    return df.apply(pd.to_numeric, errors="coerce").dropna(how="all")


def load_gpl_mapping(gpl_path):
    mapping = {}
    if not gpl_path.exists():
        return mapping
    with gzip.open(gpl_path, "rt", encoding="utf-8", errors="ignore") as fh:
        header = None
        for line in fh:
            fields = line.rstrip("\n").split("\t")
            if not fields or not fields[0].strip():
                continue
            f0 = fields[0].strip()
            if f0.startswith(("#", "^", "!")):
                continue
            if header is None:
                if "ID" in fields and ("Gene symbol" in fields or "Gene Symbol" in fields):
                    header = fields
                continue
            d = dict(zip(header, fields))
            p = d.get("ID", "").strip()
            s = d.get("Gene symbol", d.get("Gene Symbol", "")).strip()
            if p and s and s not in ("---", ""):
                mapping[p] = s
    return mapping


def filter_and_map(df, mapping):
    probes, symbols = [], []
    for probe in df.index:
        if probe.startswith("AFFX") or probe not in mapping:
            continue
        sym = mapping[probe]
        if "///" in sym:            # multi-gene probes
            continue
        if sym in SEX_LINKED:       # sex-linked confounders
            continue
        probes.append(probe); symbols.append(sym)
    if not probes:
        raise SystemExit("[err] no probes mapped after filtering.")
    sub = df.loc[probes].copy()
    sub.index = pd.Index(symbols, name="gene_symbol")
    return sub.groupby(level=0).mean()


def compute_relevance(df):
    return df.std(axis=1) / (df.mean(axis=1).abs() + 1e-9)


def is_pathway_gene(sym):
    for members in PATHWAY_MEMBERS.values():
        if sym in members:
            return 1.0
    return 0.0


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--ad", default="data/raw/GSE5281_series_matrix.txt.gz")
    p.add_argument("--pd", default="data/raw/GSE49036_series_matrix.txt.gz")
    p.add_argument("--ad-annot", default="data/annotation/GPL570.annot.gz")
    p.add_argument("--pd-annot", default="")
    p.add_argument("--out", default="data/processed/features.csv")
    args = p.parse_args()

    ad = _read_series_matrix(Path(args.ad))
    pd_ = _read_series_matrix(Path(args.pd))
    ad_map = load_gpl_mapping(Path(args.ad_annot))
    pd_annot = Path(args.pd_annot) if args.pd_annot else Path(args.ad_annot)
    pd_map = load_gpl_mapping(pd_annot) if pd_annot.exists() else ad_map
    print(f"[info] AD probes: {len(ad)}  mapped: {len(ad_map)}")
    print(f"[info] PD probes: {len(pd_)}  mapped: {len(pd_map)}")

    ad = filter_and_map(ad, ad_map)
    pd_ = filter_and_map(pd_, pd_map)
    print(f"[info] AD genes after filter: {len(ad)}  PD: {len(pd_)}")

    ad_s = compute_relevance(ad).rename("ad_score")
    pd_s = compute_relevance(pd_).rename("pd_score")
    m = pd.concat([ad_s, pd_s], axis=1)
    m["shared_score"] = np.minimum(m["ad_score"], m["pd_score"])
    m["in_pathway"] = [is_pathway_gene(s) for s in m.index]
    # Priority score: pathway genes get a strong bonus so they survive top-k
    m["priority"] = m["shared_score"] + 2.0 * m["in_pathway"]
    m = m.dropna().sort_values("priority", ascending=False)

    Path(args.out).parent.mkdir(parents=True, exist_ok=True)
    m.to_csv(args.out)
    print(f"[ok ] {len(m)} genes -> {args.out}")
    print(f"[info] pathway genes present: {int(m['in_pathway'].sum())}")
    print("\nTop 15 by priority:")
    print(m.head(15)[["ad_score","pd_score","shared_score","in_pathway"]].to_string())


if __name__ == "__main__":
    main()
