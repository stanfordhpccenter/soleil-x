#!/bin/bash -eu

# Inputs
export QUEUE="${QUEUE:-}"
export USE_CUDA="${USE_CUDA:-1}"

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

function get_num_ranks {
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
# Total number of ranks is the sum of all sample rank requirements.
MINUTES=0
NUM_RANKS=0
function parse_config {
    _MINUTES="$(get_walltime "$1")"
    MINUTES=$(( MINUTES > _MINUTES ? MINUTES : _MINUTES ))
    _NUM_RANKS="$(get_num_ranks "$1")"
    NUM_RANKS=$(( NUM_RANKS + _NUM_RANKS ))
}
for (( i = 1; i <= $#; i++ )); do
    if [[ "${!i}" == "-i" ]] && (( $i < $# )); then
        j=$((i+1))
        parse_config "${!j}"
    elif [[ "${!i}" == "-I" ]] && (( $i < $# )); then
        j=$((i+1))
        while IFS='' read -r LINE || [[ -n "$LINE" ]]; do
            parse_config "$LINE"
        done < "${!j}"
    fi
done
if (( NUM_RANKS < 1 )); then
    quit "No configuration files provided"
fi
WALLTIME="$(printf "%02d:%02d:00" $((MINUTES/60)) $((MINUTES%60)))"
export NUM_RANKS

export LD_LIBRARY_PATH="${LD_LIBRARY_PATH:-}:$LEGION_DIR/bindings/regent/"
if [[ ! -z "${HDF_ROOT:-}" ]]; then
    export LD_LIBRARY_PATH="${LD_LIBRARY_PATH:-}:$HDF_ROOT/lib"
fi

export CURR_DIR="$(pwd)"

###############################################################################

function run_titan {
    export QUEUE="${QUEUE:-debug}"
    qsub -V \
        -l nodes="$NUM_RANKS" -l walltime="$WALLTIME" -q "$QUEUE" \
        "$SOLEIL_DIR"/src/titan.pbs
}

function run_summit {
    export QUEUE="${QUEUE:-batch}"
    EXCLUDED="$(sed 's/^/hname!=/' "$SOLEIL_DIR"/src/blacklist/summit.txt |
                paste -sd '&' | sed 's/&/ && /g')"
    NUM_NODES="$(( NUM_RANKS/2 + NUM_RANKS%2 ))"
    # TODO: Support blacklist -- this requires Expert Mode (-csm y).
    bsub -J soleil -P CSC275IACCARINO -alloc_flags smt4 \
        -nnodes "$NUM_NODES" -W "$MINUTES" -q "$QUEUE" \
        "$SOLEIL_DIR"/src/summit.lsf
}

function run_certainty {
    export QUEUE="${QUEUE:-gpu}"
    RESOURCES=
    if [[ "$QUEUE" == "gpu" ]]; then
	RESOURCES="gpu:4"
    fi
    EXCLUDED="$(paste -sd ',' "$SOLEIL_DIR"/src/blacklist/certainty.txt)"
    sbatch --export=ALL \
        -N "$NUM_RANKS" -t "$WALLTIME" -p "$QUEUE" --gres="$RESOURCES" \
	--exclude="$EXCLUDED" \
        "$SOLEIL_DIR"/src/certainty.slurm
}

function run_sapling {
    # Allocate up to 4 nodes, from n0000 up to n0003
    if (( NUM_RANKS > 4 )); then quit "Too many nodes requested"; fi
    NODES=n0000
    for (( i = 1; i < NUM_RANKS; i++ )); do
        NODES="$NODES,n000$i"
    done
    GPU_OPTS=
    if [[ "$USE_CUDA" == 1 ]]; then
        GPU_OPTS="-ll:gpu 1 -ll:fsize 2048"
    fi
    mpiexec -H "$NODES" --bind-to none -x LD_LIBRARY_PATH \
        "$SOLEIL_DIR"/src/soleil.exec $ARGS $GPU_OPTS \
        -ll:cpu 0 -ll:ocpu 1 -ll:onuma 0 -ll:okindhack -ll:othr 8 \
        -ll:csize 20000
    # Resources:
    # 40230MB RAM per node
    # 2 NUMA domains per node
    # 4 cores per NUMA domain
    # 2 Tesla C2070 GPUs per node
    # 6GB FB per GPU
}

function run_local {
    "$SOLEIL_DIR"/src/soleil.exec $ARGS \
        -ll:cpu 0 -ll:ocpu 1 -ll:onuma 0 -ll:okindhack -ll:othr 3 \
        -ll:csize 9000
}

###############################################################################

if [[ "$(uname -n)" == *"titan"* ]]; then
    run_titan
elif [[ "$(dnsdomainname)" == *"summit"* ]]; then
    run_summit
elif [[ "$(uname -n)" == *"certainty"* ]]; then
    run_certainty
elif [[ "$(uname -n)" == *"sapling"* ]]; then
    run_sapling
else
    run_local
fi
