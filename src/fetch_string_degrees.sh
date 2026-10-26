#!/usr/bin/env bash
set -euo pipefail

OUT="data/reference/string_degrees.json"
mkdir -p "$(dirname "$OUT")"

if [[ -f "$OUT" ]]; then
    echo "[skip] $OUT exists"
    exit 0
fi

GENES="NDUFA1 NDUFB2 SDHA SDHB COX5A ATP5A1 ATP5B UQCRC1 UQCRC2 CYCS \
IL6 TNF IL1B CXCL8 CCL2 IL10 IL4 TGFB1 NLRP3 TLR4 \
HSPA1A HSPA8 PSMD1 UBB SQSTM1 HSPB1 DNAJB1 BAG3 MAP1LC3B ATG5 \
RAB3A SNAP25 SYN1 CALM3 CDC42 STX1A VAMP2 SNAP29 RAB5A DNM1"

echo "{" > "$OUT"
first=1
for gene in $GENES; do
    n=$(curl -sS "https://string-db.org/api/tsv/interaction_partners?identifiers=${gene}&species=9606&required_score=700&limit=0" \
        | tail -n +2 | wc -l)
    printf '[fetch] %-10s %5d partners\n' "$gene" "$n" >&2
    if [[ $first -eq 1 ]]; then
        printf '  "%s": %d' "$gene" "$n" >> "$OUT"
        first=0
    else
        printf ',\n  "%s": %d' "$gene" "$n" >> "$OUT"
    fi
    sleep 0.3
done
printf '\n}\n' >> "$OUT"
echo "[ok ] -> $OUT"
