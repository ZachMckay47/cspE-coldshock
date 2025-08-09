#!/usr/bin/env bash
set -euo pipefail
CFG=config/cspE_config.json
FASTA=$(jq -r .fasta "$CFG")
SD_START=$(jq -r .sd_start "$CFG")
SD_END=$(jq -r .sd_end "$CFG")
B15=$(jq -r .temp_baseline "$CFG")
EXTRA=$(jq -r '.temps_extra[]' "$CFG" 2>/dev/null || true)
OUT=results/baseline_sd_accessibility.csv

echo "temp_c,p_unpaired_SD" > "$OUT"

for T in $B15 $EXTRA; do
  python3 scripts/sd_analysis.py \
    --fasta "$FASTA" \
    --sd-start "$SD_START" --sd-end "$SD_END" \
    --temp "$T" \
    --run-plfold >/dev/null

  P=$(awk -v end="$SD_END" '($1==end){print $7}' results/cspE_transcript_fold_lunp)
  echo "$T,$P" >> "$OUT"
done

echo "Wrote $OUT"
column -s, -t "$OUT" | sed 's/^/  /'
