#!/bin/bash -eu

AFTER_FIRST=0
for FILE in "$@"; do
    if [[ "$AFTER_FIRST" == 1 ]]; then
	tail -n +3 "$FILE"
    else
	cat "$FILE"
	AFTER_FIRST=1
    fi
done
