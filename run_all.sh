#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")"
# 실행 예시: bash run_all.sh strains.csv data reproduced_results 3
META="${1:-strains.csv}"
DATA="${2:-data}"
OUT="${3:-results}"
REPS="${4:-3}"
# 기능 검사 → 유전체 분석 → 종 간 비교 → 그래프 저장 → 결과 검사
python test_all.py
python analyze.py --metadata "$META" --data "$DATA" --out "$OUT" --reps "$REPS" --seed 42
python compare.py "$OUT"
python plot.py "$OUT"
python check_output.py "$META" "$OUT" "$REPS"
echo "Pipeline completed: $OUT"
