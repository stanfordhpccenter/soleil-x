#!/bin/bash -eu

LEVELS_DAT="$1"
JOBID="$2"

LF="$(grep -o 'lf[0-9]*' "$JOBID".out | head -1)"
LFID="${LF#lf}"
LINE=$(( LFID + 2 ))
tail -n+"$LINE" "$LEVELS_DAT" | head -n1
