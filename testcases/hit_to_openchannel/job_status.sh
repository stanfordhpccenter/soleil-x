#!/bin/bash -eu

for JOBOUT in "$@"; do
    JOBID="${JOBOUT%.out}"
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
    elif grep -q 'Ran out of space while copying particles from other section' "$JOBOUT"; then
        echo "$JOBID: channel section overflow"
    elif grep -q Successfully "$JOBOUT"; then
        echo -n "$JOBID: done"
	FAILURES=`tail -q -n 1 "$SCRATCH"/"$JOBID"/sample*/console.txt | grep 'nan' | wc -l`
	echo " ($FAILURES/32 cases diverged)"
    elif grep -q summary "$JOBOUT"; then
        echo "$JOBID: OTHER ERROR"
    else
        echo -n "$JOBID: running"
	FAILURES=`tail -q -n 1 "$SCRATCH"/"$JOBID"/sample*/console.txt | grep 'nan' | wc -l`
	echo " ($FAILURES/32 cases diverged)"
    fi
done
