#!/bin/bash -eu

TARGET_ITERS=5000
MAX_TIME=0.0
MIN_TIME=1000000.0
for CONSOLE in "$1"/sample*/console.txt; do
    LAST_LINE="$(tail -n 1 "$CONSOLE")"
    ITERS="$(echo "$LAST_LINE" | awk '{print $1}')"
    if [[ "$ITERS" -ne "$TARGET_ITERS" ]]; then
	echo "Unexpected iteration number"
	exit 1
    fi
    TIME="$(echo "$LAST_LINE" | awk '{print $3}')"
    if (( $(echo "$TIME > $MAX_TIME" | bc -l) )); then
        MAX_TIME="$TIME"
    fi
    if (( $(echo "$TIME < $MIN_TIME" | bc -l) )); then
        MIN_TIME="$TIME"
    fi
done
echo -e "$MIN_TIME\t$MAX_TIME"
