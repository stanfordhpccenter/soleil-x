#!/bin/bash -eu

for JOBID in "$@"; do
    if [[ ! -e "$JOBID".out ]]; then
        echo "$JOBID: not started"
    elif grep -q 'CUDA_ERROR_OUT_OF_MEMORY' "$JOBID".out; then
        echo "$JOBID: cuda oom"
    elif grep -q 'TERM_ADMIN' "$JOBID".out; then
        echo "$JOBID: killed by admin"
    elif grep -q 'TERM_OWNER' "$JOBID".out; then
        echo "$JOBID: manually killed"
    elif grep -q 'TERM_RUNLIMIT' "$JOBID".out; then
        echo "$JOBID: timeout"
    elif grep -q 'Cannot open your job file' "$JOBID".out; then
        echo "$JOBID: cannot open job file"
    elif grep -q 'code 247' "$JOBID".out; then
        echo "$JOBID: code 247"
    elif grep -q 'code 255' "$JOBID".out; then
        echo "$JOBID: code 255"
    elif grep -q 'JSM server is not responding' "$JOBID".out; then
        echo "$JOBID: jsm server not responding"
    elif grep -q Successfully "$JOBID".out; then
        echo -n "$JOBID: done"
	FAILURES=`tail -q -n 1 "$SCRATCH"/"$JOBID"/sample*/console.txt | grep 'nan' | wc -l`
	echo " ($FAILURES/32 cases diverged)"
    elif grep -q summary "$JOBID".out; then
        echo "$JOBID: OTHER ERROR"
    else
        echo -n "$JOBID: running"
	FAILURES=`tail -q -n 1 "$SCRATCH"/"$JOBID"/sample*/console.txt | grep 'nan' | wc -l`
	echo " ($FAILURES/32 cases diverged)"
    fi
done
