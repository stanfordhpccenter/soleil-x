#!/bin/bash -eu

SOLEIL_SRC="$(cd "$(dirname "$(perl -MCwd -le 'print Cwd::abs_path(shift)' "${BASH_SOURCE[0]}")")" && pwd)"
cd "$SOLEIL_SRC"

export LD_LIBRARY_PATH="${LD_LIBRARY_PATH:-}:$LEGION_DIR/bindings/regent/"
if [ ! -z "${HDF_ROOT:-}" ]; then
    export LD_LIBRARY_PATH="${LD_LIBRARY_PATH:-}:$HDF_ROOT/lib"
fi

./dom_host.exec "$@" -ll:csize 9000
