#!/bin/bash -eu

SOLEIL_SRC="$(cd "$(dirname "$(perl -MCwd -le 'print Cwd::abs_path(shift)' "${BASH_SOURCE[0]}")")" && pwd)"
cd "$SOLEIL_SRC"

if [[ $(uname -n) == *"titan"* ]]; then
    qsub soleil.pbs
else
    LD_LIBRARY_PATH="$LEGION_DIR"/bindings/terra/ ./soleil.exec \
        -i ../testcases/cavity/cavity_32x32.json
fi
