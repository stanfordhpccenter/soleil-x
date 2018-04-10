#!/bin/bash -eu

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

# Command-line arguments are passed directly to the job script.
for (( i = 1; i <= $#; i++ )); do
    if [[ "${!i}" != "${!i%[[:space:]]*}" ]]; then
        quit "Cannot handle spaces in command line arguments"
    fi
done
export ARGS=$@

# Total wall-clock time is the maximum across all samples.
# Total number of nodes is the sum of all sample node requirements.
MINUTES=0
export NUM_NODES=0
for (( i = 1; i <= $#; i++ )); do
    if [[ "${!i}" == "-i" ]] && (( $i < $# )); then
        j=$((i+1))
        _MINUTES="$(get_walltime "${!j}")"
        MINUTES=$(( MINUTES > _MINUTES ? MINUTES : _MINUTES ))
        _NUM_NODES="$(get_num_nodes "${!j}")"
        export NUM_NODES=$(( NUM_NODES + _NUM_NODES ))
    fi
done
if (( NUM_NODES < 1 )); then
    quit "Usage: $0 -i <config1.json> [-i <config2.json> ...]"
fi
WALLTIME="$(printf "%02d:%02d:00" $((MINUTES/60)) $((MINUTES%60)))"

export LD_LIBRARY_PATH="${LD_LIBRARY_PATH:-}:$LEGION_DIR/bindings/regent/"
if [[ ! -z "${HDF_ROOT:-}" ]]; then
    export LD_LIBRARY_PATH="${LD_LIBRARY_PATH:-}:$HDF_ROOT/lib"
fi

###############################################################################

function run_titan {
    export QUEUE="${QUEUE:-debug}"
    if [[ "$QUEUE" == "batch" ]]; then
	true
    elif [[ "$QUEUE" == "killable" ]]; then
	true
    elif [[ "$QUEUE" == "debug" ]]; then
        true
    else
	quit "Unrecognized queue $QUEUE"
    fi
    export CURR_DIR="$(pwd)"
    qsub -v LD_LIBRARY_PATH,ARGS,NUM_NODES,CURR_DIR \
        -l nodes="$NUM_NODES" -l walltime="$WALLTIME" -q "$QUEUE" \
        "$SOLEIL_DIR"/src/titan.pbs
}

function run_certainty {
    export QUEUE="${QUEUE:-gpu}"
    # TODO: Only supporting blacklist on gpu queue
    if [[ "$QUEUE" == "default" ]]; then
	if (( MINUTES > 2 * 24 * 60 )); then quit "Walltime too long"; fi
	NODES="$NUM_NODES:ppn=32"
    elif [[ "$QUEUE" == "debug" ]]; then
	if (( MINUTES > 4 * 60 )); then quit "Walltime too long"; fi
	NODES="$NUM_NODES:ppn=24"
    elif [[ "$QUEUE" == "gpu" ]]; then
        # HACK: Torque doesn't support node exclusion, so we just list all free
        # GPU nodes, pick $NUM_NODES that aren't on the blacklist, and request
        # those specifically.
        NUM_AVAIL="$(pbsnodes -l free | grep gpu | awk '{print $1}' | sort |
                     comm -23 - "$SOLEIL_DIR"/src/blacklist/certainty.txt |
                     wc -l)"
        if (( NUM_AVAIL < NUM_NODES )); then
	    quit "Too many nodes requested"
	fi
        NODES="$(pbsnodes -l free | grep gpu | awk '{print $1}' | sort |
                 comm -23 - "$SOLEIL_DIR"/src/blacklist/certainty.txt |
                 head -n $NUM_NODES | paste -sd '+' -)"
    elif [[ "$QUEUE" == "largemem" ]]; then
	if (( MINUTES > 7 * 24 * 60 )); then quit "Walltime too long"; fi
	NODES="$NUM_NODES:ppn=24"
    elif [[ "$QUEUE" == "long" ]]; then
	NODES="$NUM_NODES:ppn=24"
    elif [[ "$QUEUE" == "quick" ]]; then
	if (( MINUTES > 24 * 60 )); then quit "Walltime too long"; fi
	NODES="$NUM_NODES:ppn=16"
    else
	quit "Unrecognized queue $QUEUE"
    fi
    qsub -v LD_LIBRARY_PATH,ARGS,NUM_NODES,QUEUE \
        -l nodes="$NODES" -l walltime="$WALLTIME" -q "$QUEUE" \
        "$SOLEIL_DIR"/src/certainty.pbs
}

function run_sapling {
    # Allocate up to 4 nodes, from n0000 up to n0003
    if (( NUM_NODES > 4 )); then quit "Too many nodes requested"; fi
    NODES=n0000
    for (( i = 1; i < NUM_NODES; i++ )); do
        NODES="$NODES,n000$i"
    done
    mpiexec -H "$NODES" --bind-to none -x LD_LIBRARY_PATH \
        "$SOLEIL_DIR"/src/soleil.exec $ARGS \
        -ll:cpu 0 -ll:ocpu 1 -ll:onuma 0 -ll:okindhack -ll:othr 8 -ll:gpu 1 \
        -ll:csize 20000 -ll:fsize 2048
}

function run_local {
    "$SOLEIL_DIR"/src/soleil.exec $ARGS \
        -ll:cpu 0 -ll:ocpu 1 -ll:onuma 0 -ll:okindhack -ll:othr 3 \
        -ll:csize 9000
}

###############################################################################

if [[ "$(uname -n)" == *"titan"* ]]; then
    run_titan
elif [[ "$(uname -n)" == *"certainty"* ]]; then
    run_certainty
elif [[ "$(uname -n)" == *"sapling"* ]]; then
    run_sapling
else
    run_local
fi
