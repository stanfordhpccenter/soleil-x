#!/bin/bash -eu

SOLEIL_SRC="$(cd "$(dirname "$(perl -MCwd -le 'print Cwd::abs_path(shift)' "${BASH_SOURCE[0]}")")" && pwd)"

cd "$SOLEIL_SRC"

LD_LIBRARY_PATH=.:"$LEGION_PATH"/bindings/terra/ ./soleil.exec -i ../testcases/cavity/cavity_32x32.json
