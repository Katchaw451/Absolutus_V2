"""preprocess.py - map probes, filter confounders, score with network degree."""
from __future__ import annotations
import argparse, gzip, json, math
from pathlib import Path
import numpy as np, pandas as pd

SEX_LINKED = {
    "XIST", "RPS4Y1", "RPS4Y2", "EIF1AY", "USP9Y", "KDM5D", "DDX3Y",
    "UTY", "ZFY", "TXLNGY", "NLGN4Y", "TMSB4Y", "PRKY", "AMELY",
}

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
    if not Path(gpl_path).exists():
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
                if "ID" in fields and "Gene symbol" in fields:
                    header = fields
                continue
            d = dict(zip(header, fields))
            probe = d.get("ID", "").strip()
            sym = d.get("Gene symbol", "").strip()
            if probe and sym and sym != "---":
                mapping[probe] = sym
    return mapping


def resolve_symbol(sym):
    """For multi-gene probes, take the first symbol. Else return as-is."""
    if "///" in sym:
        return sym.split("///")[0].strip()
    return sym


def filter_and_map(df, mapping):
    probes, symbols = [], []
    for probe in df.index:
        if probe.startswith("AFFX") or probe not in mapping:
            continue
        sym = resolve_symbol(mapping[probe])
        if sym in SEX_LINKED or sym == "":
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
    return 1.0 if any(sym in m for m in PATHWAY_MEMBERS.values()) else 0.0


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--ad", default="data/raw/GSE5281_series_matrix.txt.gz")
    p.add_argument("--pd", default="data/raw/GSE49036_series_matrix.txt.gz")
    p.add_argument("--ad-annot", default="data/annotation/GPL570.annot.gz")
    p.add_argument("--degrees", default="data/reference/string_degrees.json")
    p.add_argument("--out", default="data/processed/features.csv")
    p.add_argument("--degree-weight", type=float, default=0.5,
                   help="Weight of log-normalised STRING degree in priority.")
    args = p.parse_args()

    ad = _read_series_matrix(Path(args.ad))
    pd_ = _read_series_matrix(Path(args.pd))
    mapping = load_gpl_mapping(Path(args.ad_annot))

    ad = filter_and_map(ad, mapping)
    pd_ = filter_and_map(pd_, mapping)
    print(f"[info] genes after filtering: {len(ad)}")

    ad_s = compute_relevance(ad).rename("ad_score")
    pd_s = compute_relevance(pd_).rename("pd_score")
    m = pd.concat([ad_s, pd_s], axis=1)
    m["shared_score"] = np.minimum(m["ad_score"], m["pd_score"])
    m["in_pathway"] = [is_pathway_gene(s) for s in m.index]

    degrees = {}
    if Path(args.degrees).exists():
        degrees = json.loads(Path(args.degrees).read_text())
        print(f"[info] loaded {len(degrees)} STRING degrees")
    m["string_degree"] = [degrees.get(g, 0) for g in m.index]
    max_deg = max(m["string_degree"].max(), 1)
    m["degree_score"] = m["string_degree"].apply(
        lambda d: math.log1p(d) / math.log1p(max_deg)
    )
    m["priority"] = (m["shared_score"]
                     + 2.0 * m["in_pathway"]
                     + args.degree_weight * m["degree_score"])

    m = m.dropna().sort_values("priority", ascending=False)
    Path(args.out).parent.mkdir(parents=True, exist_ok=True)
    m.to_csv(args.out)
    print(f"[ok ] {len(m)} genes -> {args.out}")
    print(f"[info] CALM3 in features: {'CALM3' in m.index}")
    print("\nTop 20 by priority:")
    print(m.head(20)[["ad_score","pd_score","shared_score",
                     "in_pathway","string_degree","degree_score","priority"]].to_string())


if __name__ == "__main__":
    main()
