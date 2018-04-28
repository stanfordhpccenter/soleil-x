#!/bin/bash -eu

if [ "$(uname)" == "Darwin" ]; then
    export DYLD_LIBRARY_PATH="${DYLD_LIBRARY_PATH:-}:$LEGION_DIR/bindings/regent/"
    if [[ ! -z "${HDF_ROOT:-}" ]]; then
        export DYLD_LIBRARY_PATH="${DYLD_LIBRARY_PATH:-}:$HDF_ROOT/lib"
    fi
else
    export LD_LIBRARY_PATH="${LD_LIBRARY_PATH:-}:$LEGION_DIR/bindings/regent/"
    if [[ ! -z "${HDF_ROOT:-}" ]]; then
        export LD_LIBRARY_PATH="${LD_LIBRARY_PATH:-}:$HDF_ROOT/lib"
    fi
fi

"$SOLEIL_DIR"/src/dom_host.exec "$@" -ll:csize 9000
