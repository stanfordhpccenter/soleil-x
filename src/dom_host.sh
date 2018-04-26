#!/bin/bash -eu

export LD_LIBRARY_PATH="${LD_LIBRARY_PATH:-}:$LEGION_DIR/bindings/regent/"
export DYLD_LIBRARY_PATH="${LD_LIBRARY_PATH:-}:$LEGION_DIR/bindings/regent/"
if [ ! -z "${HDF_ROOT:-}" ]; then
    export LD_LIBRARY_PATH="${LD_LIBRARY_PATH:-}:$HDF_ROOT/lib"
fi

"$SOLEIL_DIR"/src/dom_host.exec "$@" -ll:csize 9000
