#!/bin/bash -eu

SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null && pwd )"

if [ $# -ne 4 ]; then
    echo "Usage: $(basename "${BASH_SOURCE[0]}") num_cases max_ftt lf_dir out_dir"
    exit 1
fi
NUM_CASES="$1"
MAX_FTT="$2"
LF_DIR="$3"
OUT_DIR="$4"

echo -e "AvgFluidT\tAvgParticleT\tAvgCellOfParticleT"
for I in $(seq 0 "$((NUM_CASES-1))" ); do
    "$SCRIPT_DIR"/time_average.py "$MAX_FTT" "$LF_DIR"/case"$I".json "$OUT_DIR"/sample"$((I*2+1))"/probe0.csv
done
