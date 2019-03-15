#!/bin/bash -eu

SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null && pwd )"

COUNT_FAILURES="${COUNT_FAILURES:-0}"
AVERAGE="${AVERAGE:-0}"
RESTART="${RESTART:-0}"
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

    function average() {
        if [[ "$AVERAGE" == 1 ]]; then
            echo -n ", averaging"
            if [[ ! -e "$LAUNCH_DIR"/"$JOBID".csv ]]; then
                "$SCRIPT_DIR"/time_average_all.sh "$NUM_CASES" "$LAUNCH_DIR" "$OUT_DIR" > "$LAUNCH_DIR"/"$JOBID".csv
            fi
        fi
    }
    function count_failures() {
        if [[ "$COUNT_FAILURES" == 1 ]]; then
            echo -n ", samples diverged: "
            FAILURES=`tail -q -n 1 "$OUT_DIR"/sample*/console.txt | grep 'nan' | wc -l`
            echo -n "$FAILURES"
        fi
    }
    function restart() {
        if [[ "$RESTART" == 1 ]]; then
            echo ", restarting"
            # Create a placeholder jobout, while the rerun job is in the queue
            touch "$LAUNCH_DIR"/"$(( JOBID + 1 ))".out
            cd "$LAUNCH_DIR"
            RANKS_PER_NODE=4 "$SOLEIL_DIR"/src/soleil.sh $(echo -m\ case{0..31}.json)
            cd ..
        else
	    echo
	fi
    }

    if (( "$#" > 1 )); then
        echo -n "$JOBOUT: "
    fi
    if [[ ! -s "$JOBOUT" ]]; then
        echo "not started"
    elif grep -q 'CUDA_ERROR_OUT_OF_MEMORY' "$JOBOUT"; then
        echo "cuda oom"
    elif grep -q 'TERM_ADMIN' "$JOBOUT"; then
        echo -n "killed by admin"
        restart
    elif grep -q 'TERM_OWNER' "$JOBOUT"; then
        echo "manually killed"
    elif grep -q 'TERM_RUNLIMIT' "$JOBOUT"; then
        echo "timeout"
    elif grep -q 'Ran out of space while copying particles from other section' "$JOBOUT"; then
        echo "channel section overflow"
    elif grep -q 'Cannot open your job file' "$JOBOUT"; then
        echo "cannot open jobfile"
    elif grep -q Successfully "$JOBOUT"; then
        echo -n "done"
        average
        count_failures
        echo
    elif grep -q summary "$JOBOUT"; then
        echo "OTHER ERROR"
    else
        echo -n "running"
        count_failures
        echo
    fi
done
