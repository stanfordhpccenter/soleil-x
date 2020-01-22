#!/bin/bash -eu

TARGET_ITER=5000
TARGET_LINE="$(( TARGET_ITER + 2 ))"

I=1
while true; do
    CONSOLE="$1"/sample"$I"/console.txt
    if [[ ! -e "$CONSOLE" ]]; then
	break
    fi
    sed "${TARGET_LINE}q;d" "$CONSOLE"
    I="$(( I + 2 ))"
done
