#!/bin/bash -eu

SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null && pwd )"

if [ $# -ne 2 ]; then
    echo "Usage: $(basename "${BASH_SOURCE[0]}") uncertainties.dat level.json"
    exit 1
fi
UNCERTAINTIES_DAT="$1"
LEVEL_JSON="$2"

I=-1
while read -r LINE; do
    if (( "$I" >= 0 )); then
        "$SCRIPT_DIR"/make_case.py "$LEVEL_JSON" $LINE > case"$I".json
    fi
    I="$(( I + 1 ))"
done < "$UNCERTAINTIES_DAT"
