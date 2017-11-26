#!/bin/bash -eu

SOLEIL_SRC="$(cd "$(dirname "$(readlink -f "${BASH_SOURCE[0]}")")" && pwd)"

cd "$SOLEIL_SRC"

LD_LIBRARY_PATH=.:"$LEGION_PATH"/bindings/terra/ ./soleil.exec ../testcases/cavity/cavity_32x32.json

cd -
