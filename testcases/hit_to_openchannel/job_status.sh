#!/bin/bash -eu

VERBOSE="${VERBOSE:-0}"

for JOBOUT in "$@"; do
    COUNT_FAILURES=0
    JOBID="${JOBOUT%.out}"
    OUTDIR="$( head -n 1 "$JOBOUT" | awk '{print $4'} )"
    if [[ ! -e "$JOBOUT" ]]; then
        echo "$JOBID: not started"
    elif grep -q 'CUDA_ERROR_OUT_OF_MEMORY' "$JOBOUT"; then
        echo "$JOBID: cuda oom"
    elif grep -q 'TERM_ADMIN' "$JOBOUT"; then
        echo "$JOBID: killed by admin"
    elif grep -q 'TERM_OWNER' "$JOBOUT"; then
        echo "$JOBID: manually killed"
    elif grep -q 'TERM_RUNLIMIT' "$JOBOUT"; then
        echo "$JOBID: timeout"
    elif grep -q 'Cannot open your job file' "$JOBOUT"; then
        echo "$JOBID: cannot open job file"
    elif grep -q 'code 247' "$JOBOUT"; then
        echo "$JOBID: code 247"
    elif grep -q 'code 255' "$JOBOUT"; then
        echo "$JOBID: code 255"
    elif grep -q 'JSM server is not responding' "$JOBOUT"; then
        echo "$JOBID: jsm server not responding"
    elif grep -q 'id.is_event' "$JOBOUT"; then
        echo "$JOBID: realm assert: id.is_event"
    elif grep -q 'Ran out of space while copying particles from other section' "$JOBOUT"; then
        echo "$JOBID: channel section overflow"
    elif grep -q Successfully "$JOBOUT"; then
        echo -n "$JOBID: done"
	COUNT_FAILURES=1
    elif grep -q summary "$JOBOUT"; then
        echo "$JOBID: OTHER ERROR"
    else
        echo -n "$JOBID: running"
	COUNT_FAILURES=1
    fi
    if [[ "$VERBOSE" == 1 && "$COUNT_FAILURES" == 1 ]]; then
        FAILURES=`tail -q -n 1 "$OUTDIR"/sample*/console.txt | grep 'nan' | wc -l`
        echo " ($FAILURES cases diverged)"
    else
        echo
    fi
done
