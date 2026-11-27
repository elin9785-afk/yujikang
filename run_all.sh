#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")"
META="${1:-strains.csv}"
DATA="${2:-data}"
OUT="${3:-results}"
REPS="${4:-3}"
python test_all.py
python analyze.py --metadata "$META" --data "$DATA" --out "$OUT" --reps "$REPS" --seed 42
python compare.py "$OUT"
python plot.py "$OUT"
python check_output.py "$META" "$OUT" "$REPS"
echo "Pipeline completed: $OUT"
