#!/bin/bash -eu

SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null && pwd )"

if [ $# -lt 2 ]; then
    echo "Usage: $(basename "${BASH_SOURCE[0]}") uncertainties.dat level.json [stagger=N]"
    exit 1
fi
UNCERTAINTIES_DAT="$1"
LEVEL_JSON="$2"
STAGGER="$3"

I=-1
while read -r LINE; do
    if (( "$I" >= 0 )); then
        "$SCRIPT_DIR"/make_case.py "$LEVEL_JSON" $LINE $STAGGER > case"$I".json
    fi
    I="$(( I + 1 ))"
done < "$UNCERTAINTIES_DAT"
