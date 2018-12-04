#!/bin/bash -eu

SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null && pwd )"

if [ $# -ne 2 ]; then
    echo "Usage: $(basename "${BASH_SOURCE[0]}") <dat_file> <json_template>"
    exit 1
fi
DAT_FILE="$1"
JSON_TEMPLATE="$2"

I=0
while read -r LINE; do
    "$SCRIPT_DIR"/make_config.py $LINE --json_template "$JSON_TEMPLATE" > case"$I".json
    I="$(( I + 1 ))"
done < "$DAT_FILE"
