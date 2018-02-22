#!/bin/bash -eu

SOLEIL_SRC="$(cd "$(dirname "$(perl -MCwd -le 'print Cwd::abs_path(shift)' "${BASH_SOURCE[0]}")")" && pwd)"
cd "$SOLEIL_SRC"

export LD_LIBRARY_PATH="$LEGION_DIR"/bindings/regent/
if [ ! -z "${HDF_ROOT:-}" ]; then
    export LD_LIBRARY_PATH="$LD_LIBRARY_PATH:$HDF_ROOT/lib"
fi

if [[ $(uname -n) == *"titan"* ]]; then
    qsub -v LD_LIBRARY_PATH titan.pbs
elif [[ $(uname -n) == *"certainty"* ]]; then
    qsub -v LD_LIBRARY_PATH certainty.pbs
elif [[ $(uname -n) == *"sapling"* ]]; then
    mpiexec -H n0000 --bind-to none -x LD_LIBRARY_PATH \
        ./soleil.exec \
        -i ../testcases/tgv_64x64x64.json \
        -ll:cpu 0 -ll:ocpu 1 -ll:onuma 0 -ll:okindhack -ll:othr 8 -ll:gpu 1
else
    ./soleil.exec \
        -i ../testcases/tgv_64x64x64.json \
        -ll:cpu 0 -ll:ocpu 1 -ll:onuma 0 -ll:okindhack -ll:othr 3 \
        -ll:csize 9000
fi
