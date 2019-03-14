#!/bin/bash -eu

SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null && pwd )"

COUNT_FAILURES="${COUNT_FAILURES:-0}"
AVERAGE="${AVERAGE:-0}"
NUM_CASES="${NUM_CASES:-32}"

for ARG in "$@"; do
    if [[ -d "$ARG" ]]; then
        if ! ls "$ARG"/*.out 1> /dev/null 2>&1; then
	    if (( "$#" > 1 )); then
		echo -n "$ARG: "
	    fi
            echo "not started"
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
    if (( "$#" > 1 )); then
	echo -n "$JOBOUT: "
    fi
    if [[ ! -e "$JOBOUT" ]]; then
        echo -n "not started"
    elif grep -q 'CUDA_ERROR_OUT_OF_MEMORY' "$JOBOUT"; then
        echo -n "cuda oom"
    elif grep -q 'TERM_ADMIN' "$JOBOUT"; then
        echo -n "killed by admin"
    elif grep -q 'TERM_OWNER' "$JOBOUT"; then
        echo -n "manually killed"
    elif grep -q 'TERM_RUNLIMIT' "$JOBOUT"; then
        echo -n "timeout"
    elif grep -q 'Ran out of space while copying particles from other section' "$JOBOUT"; then
        echo -n "channel section overflow"
    elif grep -q 'Cannot open your job file'  "$JOBOUT"; then
        echo -n "jobfile"
    elif grep -q Successfully "$JOBOUT"; then
        echo -n "done"
        NO_ERROR=1
        if [[ "$AVERAGE" == 1 ]]; then
            echo -n ", averaging"
            if [[ ! -e "$LAUNCH_DIR"/"$JOBID".csv ]]; then
                "$SCRIPT_DIR"/time_average_all.sh "$NUM_CASES" "$LAUNCH_DIR" "$OUT_DIR" > "$LAUNCH_DIR"/"$JOBID".csv
            fi
        fi
    elif grep -q summary "$JOBOUT"; then
        echo -n "OTHER ERROR"
    else
        echo -n "running"
        NO_ERROR=1
    fi
    if [[ "$COUNT_FAILURES" == 1 && "$NO_ERROR" == 1 ]]; then
        echo -n ", samples diverged: "
        FAILURES=`tail -q -n 1 "$OUT_DIR"/sample*/console.txt | grep 'nan' | wc -l`
        echo -n "$FAILURES"
    fi
    echo
done
