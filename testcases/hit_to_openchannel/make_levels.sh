#!/bin/bash -eu

SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null && pwd )"

if [ $# -ne 2 ]; then
    echo "Usage: $(basename "${BASH_SOURCE[0]}") levels.dat hf.json"
    exit 1
fi
LEVELS_DAT="$1"
HF_JSON="$2"

I=-1
while read -r LINE; do
    if (( "$I" >= 0 )); then
        mkdir "lf$I"
        "$SCRIPT_DIR"/make_level.py "$HF_JSON" $LINE > "lf$I/base.json"
    fi
    I="$(( I + 1 ))"
done < "$LEVELS_DAT"

echo "LF cases created in lf0/ to lf$((I-1))/."
echo "You may wish to update each lf*/base.json with more appropriate settings for:"
echo "- configs[*].Mapping.tiles"
echo "- configs[*].Mapping.tilesPerRank"
echo "- configs[*].Mapping.wallTime"
echo "- configs[*].IO.restartEveryTimeSteps"
echo "- collocateSections"
echo "before continuing."
