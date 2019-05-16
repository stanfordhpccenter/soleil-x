#!/bin/bash -eu

SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null && pwd )"

COUNT_FAILURES="${COUNT_FAILURES:-0}"
AVERAGE="${AVERAGE:-0}"
RESTART="${RESTART:-0}"
PATCH="${PATCH:-0}"
TIME_RUNS="${TIME_RUNS:-0}"

function clean_run() {
    COUNT_FAILURES=0 AVERAGE=0 RESTART=0 PATCH=0 "${BASH_SOURCE[0]}" "$1"
}

function pushd() {
    command pushd "$@" > /dev/null
}

function popd() {
    command popd "$@" > /dev/null
}

for ARG in "$@"; do
    # For directories, find the latest jobout
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

    # Possible actions
    function average() {
        if [[ "$AVERAGE" == 1 ]]; then
            echo -n ", averaging"
            if [[ ! -e "$LAUNCH_DIR"/"$JOBID".csv ]]; then
                "$SCRIPT_DIR"/time_average_all.sh 32 9999 "$LAUNCH_DIR" "$OUT_DIR" > "$LAUNCH_DIR"/"$JOBID".csv
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
            pushd "$LAUNCH_DIR"
            RANKS_PER_NODE=4 "$SOLEIL_DIR"/src/soleil.sh $(echo -m\ case{0..31}.json)
            popd
        else
            echo
        fi
    }
    function patch() {
        echo
        if [[ "$PATCH" == 1 ]]; then
            for I in {0..31}; do
                MAX_ITER=`grep maxIter "$LAUNCH_DIR"/case$I.json | head -1 | awk '{print $2}'`
                FINAL_ITER=`tail -n 1 "$OUT_DIR"/sample"$((I*2+1))"/console.txt | awk '{print $1}'`
                if (( "$MAX_ITER" != "$FINAL_ITER" )); then
                    echo -n "  case$I"
                    PATCH_LAUNCH_DIR="$LAUNCH_DIR/patch$I"
                    if [[ -d "$PATCH_LAUNCH_DIR" ]]; then
                        echo -n ", patch run: "
                        PATCH_RESULT="$(clean_run "$PATCH_LAUNCH_DIR")"
                        echo -n "$PATCH_RESULT"
                        if [[ "$PATCH_RESULT" == "done" ]]; then
                            echo -n ", postprocessing"
                            PATCH_JOBOUT="$(ls "$PATCH_LAUNCH_DIR"/*.out | tail -n 1)"
                            PATCH_OUT_DIR="$( head -n 1 "$PATCH_JOBOUT" | awk '{print $4'} )"
                            mv "$OUT_DIR"/sample"$((I*2))" "$OUT_DIR"/timeout_sample"$((I*2))"
                            mv "$OUT_DIR"/sample"$((I*2+1))" "$OUT_DIR"/timeout_sample"$((I*2+1))"
                            ln -s "$PATCH_OUT_DIR"/sample0 "$OUT_DIR"/sample"$((I*2))"
                            ln -s "$PATCH_OUT_DIR"/sample1 "$OUT_DIR"/sample"$((I*2+1))"
                        fi
                        echo
                    else
                        echo -n ", patching: "
                        mkdir "$PATCH_LAUNCH_DIR"
                        pushd "$PATCH_LAUNCH_DIR"
                        RANKS_PER_NODE=4 "$SOLEIL_DIR"/src/soleil.sh -m "../case$I.json"
                        popd
                    fi
                fi
            done
        fi
    }
    function time_runs() {
        if [[ "$TIME_RUNS" == 1 ]]; then
            echo -n ", timing"
            MAX_TIME=0.0
            for I in {0..31}; do
                CASE_TIME=`tail -n 1 "$OUT_DIR"/sample"$((I*2+1))"/console.txt | awk '{print $3}'`
                if (( $(echo "$CASE_TIME > $MAX_TIME" | bc -l) )); then
                    MAX_TIME="$CASE_TIME"
                fi
            done
            echo -n " $MAX_TIME"
        fi
    }

    # Identify job status and proceed accordingly
    if (( "$#" > 1 )); then
        echo -n "$JOBOUT: "
    fi
    if [[ ! -s "$JOBOUT" ]]; then
        echo "not started"
    elif grep -q 'TERM_RUNLIMIT' "$JOBOUT"; then
        echo -n "timeout"
        average
        count_failures
        time_runs
        patch
    elif grep -q Successfully "$JOBOUT"; then
        echo -n "done"
        average
        count_failures
        time_runs
        echo
    elif grep -q summary "$JOBOUT"; then
        echo -n "ERROR"
        restart
    else
        echo -n "running"
        count_failures
        echo
    fi
done
