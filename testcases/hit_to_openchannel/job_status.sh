#!/bin/bash -eu

SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null && pwd )"

COUNT_FAILURES="${COUNT_FAILURES:-0}"
AVERAGE="${AVERAGE:-0}"
NUM_CASES="${NUM_CASES:-32}"

for ARG in "$@"; do
    if [[ -d "$ARG" ]]; then
        if ! ls "$ARG"/*.out 1> /dev/null 2>&1; then
            echo "$ARG not started"
            continue
        fi
        JOBOUT="$(ls "$ARG"/*.out | tail -n 1)"
        LAUNCH_DIR="$ARG"
    else
        JOBOUT="$ARG"
        LAUNCH_DIR="$(dirname "$JOBOUT")"
    fi
    BASEOUT="$(basename "$JOBOUT")"
    JOBID="${BASEOUT%.out}"
    OUT_DIR="$( head -n 1 "$JOBOUT" | awk '{print $4'} )"

    NO_ERROR=0
    if [[ ! -e "$JOBOUT" ]]; then
        echo -n "$JOBOUT not started"
    elif grep -q 'CUDA_ERROR_OUT_OF_MEMORY' "$JOBOUT"; then
        echo -n "$JOBOUT cuda oom"
    elif grep -q 'TERM_ADMIN' "$JOBOUT"; then
        echo -n "$JOBOUT killed by admin"
    elif grep -q 'TERM_OWNER' "$JOBOUT"; then
        echo -n "$JOBOUT manually killed"
    elif grep -q 'TERM_RUNLIMIT' "$JOBOUT"; then
        echo -n "$JOBOUT timeout"
    elif grep -q 'Ran out of space while copying particles from other section' "$JOBOUT"; then
        echo -n "$JOBOUT channel section overflow"
    elif grep -q Successfully "$JOBOUT"; then
        echo -n "$JOBOUT done"
        NO_ERROR=1
        if [[ "$AVERAGE" == 1 ]]; then
            echo -n ", averaging"
            if [[ ! -e "$LAUNCH_DIR"/"$JOBID".csv ]]; then
                "$SCRIPT_DIR"/time_average_all.sh "$NUM_CASES" "$LAUNCH_DIR" "$OUT_DIR" > "$LAUNCH_DIR"/"$JOBID".csv
            fi
        fi
    elif grep -q summary "$JOBOUT"; then
        echo -n "$JOBOUT OTHER ERROR"
    else
        echo -n "$JOBOUT running"
        NO_ERROR=1
    fi
    if [[ "$COUNT_FAILURES" == 1 && "$NO_ERROR" == 1 ]]; then
        FAILURES=`tail -q -n 1 "$OUT_DIR"/sample*/console.txt | grep 'nan' | wc -l`
        echo " ($FAILURES samples diverged)"
    else
        echo
    fi
done
