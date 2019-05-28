#!/bin/bash -eu

SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null && pwd )"

if [ $# -ne 4 ]; then
    echo "Usage: $(basename "${BASH_SOURCE[0]}") num_cases max_ftt old_dir new_dir"
    exit 1
fi
NUM_CASES="$1"
MAX_FTT="$2"
OLD_DIR="$3"
NEW_DIR="$4"

mkdir "$NEW_DIR"
"$SCRIPT_DIR"/downsample_levels.py "$MAX_FTT" "$OLD_DIR"/levels.dat > "$NEW_DIR"/levels.dat
for OLD_LF_DIR in "$OLD_DIR"/lf*; do
    NEW_LF_DIR="$NEW_DIR"/"$(basename "$OLD_LF_DIR")"
    mkdir "$NEW_LF_DIR"
    OLD_JOBOUT="$(ls "$OLD_LF_DIR"/*.out | tail -n 1)"
    OLD_OUT_DIR="$( head -n 1 "$OLD_JOBOUT" | awk '{print $4'} )"
    "$SCRIPT_DIR"/time_average_all.sh "$NUM_CASES" "$MAX_FTT" "$OLD_LF_DIR" "$OLD_OUT_DIR" > "$NEW_LF_DIR"/averages.csv
done
