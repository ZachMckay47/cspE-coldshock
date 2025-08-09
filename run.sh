#!/usr/bin/env bash
set -euo pipefail

python3 scripts/sd_analysis.py \
  --fasta data/cspE_transcript_fold.fasta \
  --sd-start 134 --sd-end 139 \
  --temp 15 \
  --run-plfold

echo
echo "p_unpaired_SD (u=6) at 15°C:"
awk '$1==139 {print $7}' results/cspE_transcript_fold_lunp
