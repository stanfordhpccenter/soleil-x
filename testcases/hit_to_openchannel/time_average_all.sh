#!/bin/bash -eu

SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null && pwd )"

if [ $# -ne 3 ]; then
    echo "Usage: $(basename "${BASH_SOURCE[0]}") num_cases lf_dir out_dir"
    exit 1
fi
NUM_CASES="$1"
LF_DIR="$2"
OUT_DIR="$3"

echo -e "AvgFluidT\tAvgParticleT\tAvgCellOfParticleT"
for I in $(seq 0 "$((NUM_CASES-1))" ); do
    "$SCRIPT_DIR"/time_average.py "$LF_DIR"/case"$I".json "$OUT_DIR"/sample"$((I*2+1))"/probe0.csv
done
