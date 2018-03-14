#!/bin/bash -eu

SOLEIL_SRC="$(cd "$(dirname "$(perl -MCwd -le 'print Cwd::abs_path(shift)' "${BASH_SOURCE[0]}")")" && pwd)"
cd "$SOLEIL_SRC"

###############################################################################

function quit {
    echo "$1" >&2
    exit 1
}

function get_walltime {
    python -c "
import json
config = json.load(open('$1'))
print config['Mapping']['wallTime']"
}

function get_num_nodes {
    python -c "
import json
config = json.load(open('$1'))
xTiles = int(config['Mapping']['xTiles'])
yTiles = int(config['Mapping']['yTiles'])
zTiles = int(config['Mapping']['zTiles'])
print xTiles * yTiles * zTiles"
}

###############################################################################

for (( i = 1; i <= $#; i++ )); do
    if [[ "${!i}" != "${!i%[[:space:]]*}" ]]; then
        quit "Cannot handle spaces in command line arguments"
    fi
done
export ARGS=$@

WALLTIME=0
export NUM_NODES=0
for (( i = 1; i <= $#; i++ )); do
    if [[ "${!i}" == "-i" ]] && (( $i < $# )); then
        j=$((i+1))
        _WALLTIME=$(get_walltime "${!j}")
        WALLTIME=$(( WALLTIME > _WALLTIME ? WALLTIME : _WALLTIME ))
        _NUM_NODES=$(get_num_nodes "${!j}")
        export NUM_NODES=$(( NUM_NODES + _NUM_NODES ))
    fi
done
if (( NUM_NODES < 1 )); then
    quit "Usage: $0 -i <config1.json> [-i <config2.json> ...]"
fi
WALLTIME=$((WALLTIME/60)):$((WALLTIME%60)):00

export LD_LIBRARY_PATH="${LD_LIBRARY_PATH:-}:$LEGION_DIR/bindings/regent/"
if [ ! -z "${HDF_ROOT:-}" ]; then
    export LD_LIBRARY_PATH="${LD_LIBRARY_PATH:-}:$HDF_ROOT/lib"
fi

###############################################################################

if [[ $(uname -n) == *"titan"* ]]; then
    qsub -v LD_LIBRARY_PATH,ARGS,NUM_NODES \
        -l nodes=$NUM_NODES -l walltime=$WALLTIME -q debug \
        titan.pbs
elif [[ $(uname -n) == *"certainty"* ]]; then
    qsub -v LD_LIBRARY_PATH,ARGS,NUM_NODES \
        -l nodes=$NUM_NODES:ppn=24 -l walltime=$WALLTIME -q gpu \
        certainty.pbs
elif [[ $(uname -n) == *"sapling"* ]]; then
    if (( NUM_NODES > 4 )); then
        quit "Too many nodes requested"
    fi
    NODES=n0000
    for (( i = 1; i < NUM_NODES; i++ )); do
        NODES=$NODES,n000$i
    done
    mpiexec -H $NODES --bind-to none -x LD_LIBRARY_PATH \
        ./soleil.exec $ARGS \
        -ll:cpu 0 -ll:ocpu 1 -ll:onuma 0 -ll:okindhack -ll:othr 8 -ll:gpu 1 \
        -ll:csize 20000 -ll:fsize 2048
else
    ./soleil.exec $ARGS \
        -ll:cpu 0 -ll:ocpu 1 -ll:onuma 0 -ll:okindhack -ll:othr 3 \
        -ll:csize 9000
fi
