#!/bin/bash -eu

###############################################################################
# Inputs
###############################################################################

# Which executable to invoke
export EXECUTABLE="${EXECUTABLE:-}"

# Which group to submit jobs under (if a scheduler is available)
export GROUP="${GROUP:-}"

# Which queue/partition to use (if a scheduler is available)
export QUEUE="${QUEUE:-}"

# Which job to wait for before starting (if a scheduler is available)
export AFTER="${AFTER:-}"

# Whether to use GPUs (if available)
export USE_CUDA="${USE_CUDA:-1}"

# Whether to emit Legion profiler logs
export PROFILE="${PROFILE:-0}"

# Whether to freeze Legion execution on crash
export DEBUG="${DEBUG:-0}"

# How many ranks to instantiate per node
export RANKS_PER_NODE="${RANKS_PER_NODE:-1}"

# How many cores per rank to reserve for the runtime
export RESERVED_CORES="${RESERVED_CORES:-8}"

# Whether to dump additional HDF files, for debugging cross-section copying
export DEBUG_COPYING="${DEBUG_COPYING:-0}"

###############################################################################
# Helper functions
###############################################################################

function quit {
    echo "$1" >&2
    exit 1
}

function read_json {
    python2 -c "
import json
def wallTime(sample):
  return int(sample['Mapping']['wallTime'])
def numRanks(sample):
  tiles = sample['Mapping']['tiles']
  tilesPerRank = sample['Mapping']['tilesPerRank']
  xRanks = int(tiles[0]) / int(tilesPerRank[0])
  yRanks = int(tiles[1]) / int(tilesPerRank[1])
  zRanks = int(tiles[2]) / int(tilesPerRank[2])
  return xRanks * yRanks * zRanks
f = json.load(open('$1'))
if '$2' == 'single':
  print wallTime(f), numRanks(f)
elif '$2' == 'dual':
  print max(wallTime(f['configs'][0]), wallTime(f['configs'][1])), \
        max(numRanks(f['configs'][0]), numRanks(f['configs'][1]))  \
        if f['collocateSections'] else                             \
        numRanks(f['configs'][0]) + numRanks(f['configs'][1])
else:
  assert(false)"
}

###############################################################################
# Derived options
###############################################################################

# We build the command line in a string before executing it, and it's very hard
# to get this to work if the executable name contains spaces, so we punt.
if [[ "$EXECUTABLE" != "${EXECUTABLE%[[:space:]]*}" ]]; then
    quit "Cannot handle spaces in executable name"
fi
# Command-line arguments are passed directly to the job script. We need to
# accept multiple arguments separated by whitespace, and pass them through the
# environment. It is very hard to properly handle spaces in arguments in this
# mode, so we punt.
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
    read -r _MINUTES _NUM_RANKS <<<"$(read_json "$@")"
    MINUTES=$(( MINUTES > _MINUTES ? MINUTES : _MINUTES ))
    NUM_RANKS=$(( NUM_RANKS + _NUM_RANKS ))
}
for (( i = 1; i <= $#; i++ )); do
    j=$((i+1))
    if [[ "${!i}" == "-i" ]] && (( $i < $# )); then
        parse_config "${!j}" "single"
    elif [[ "${!i}" == "-m" ]] && (( $i < $# )); then
        parse_config "${!j}" "dual"
    fi
done
if (( NUM_RANKS < 1 )); then
    quit "No configuration files provided"
fi
export MINUTES
export NUM_RANKS
export WALLTIME="$(printf "%02d:%02d:00" $((MINUTES/60)) $((MINUTES%60)))"

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

# Make sure the number of requested ranks is divisible by the number of nodes.
export NUM_NODES=$(( NUM_RANKS / RANKS_PER_NODE ))
if (( NUM_RANKS % RANKS_PER_NODE > 0 )); then
    export NUM_NODES=$(( NUM_NODES + 1 ))
fi
export NUM_RANKS=$(( NUM_NODES * RANKS_PER_NODE ))

export LOCAL_RUN=0

export REALM_BACKTRACE=1

export LEGION_FREEZE_ON_ERROR="$DEBUG"

###############################################################################
# Machine-specific handling
###############################################################################

function run_titan {
    GROUP="${GROUP:-CSC188}"
    export QUEUE="${QUEUE:-batch}"
    DEPS=
    if [[ ! -z "$AFTER" ]]; then
        DEPS="-W depend=afterok:$AFTER"
    fi
    qsub -V -A "$GROUP" \
        -l nodes="$NUM_NODES" -l walltime="$WALLTIME" -q "$QUEUE" $DEPS \
        "$SOLEIL_DIR"/src/titan.pbs
}

function run_summit {
    GROUP="${GROUP:-CSC335}"
    export QUEUE="${QUEUE:-batch}"
    EXCLUDED="$(cat "$SOLEIL_DIR"/src/blacklist/summit.txt |
                sed 's/^/ \&\& (hname != /'  | sed 's/$/)/' |
                paste -sd '')"
    NUM_CORES=$(( NUM_NODES * 42 ))
    RESOURCES="1*{select[LN]span[hosts=1]} +
               $NUM_CORES*{select[CN$EXCLUDED]span[ptile=42]}"
    DEPS=
    if [[ ! -z "$AFTER" ]]; then
        DEPS="-w done($AFTER)"
    fi
    bsub -csm y -J soleil -P "$GROUP" -alloc_flags smt1 \
        -R "$RESOURCES" -W "$MINUTES" -q "$QUEUE" $DEPS \
        "$SOLEIL_DIR"/src/summit.lsf
}

function run_lassen {
    GROUP="${GROUP:-stanford}"
    export QUEUE="${QUEUE:-pbatch}"
    DEPS=
    if [[ ! -z "$AFTER" ]]; then
        DEPS="-w done($AFTER)"
    fi
    bsub -G "$GROUP" \
        -nnodes "$NUM_NODES" -W "$MINUTES" -q "$QUEUE" $DEPS \
        "$SOLEIL_DIR"/src/lassen.lsf
}

function run_pizdaint {
    export QUEUE="${QUEUE:-normal}"
    DEPS=
    if [[ ! -z "$AFTER" ]]; then
        DEPS="-d afterok:$AFTER"
    fi
    sbatch --export=ALL \
        -N "$NUM_NODES" -t "$WALLTIME" -p "$QUEUE" $DEPS \
        "$SOLEIL_DIR"/src/pizdaint.slurm
}

function run_certainty {
    export QUEUE="${QUEUE:-all}"
    EXCLUDED="$(paste -sd ',' "$SOLEIL_DIR"/src/blacklist/certainty.txt)"
    DEPS=
    if [[ ! -z "$AFTER" ]]; then
        DEPS="-d afterok:$AFTER"
    fi
    sbatch --export=ALL \
        -N "$NUM_NODES" -t "$WALLTIME" -p "$QUEUE" $DEPS \
	--exclude="$EXCLUDED" \
        "$SOLEIL_DIR"/src/certainty.slurm
}

function run_sapling {
    export QUEUE="${QUEUE:-cpu}"
    if [[ "$USE_CUDA" == 1 ]]; then
        export QUEUE=gpu
    fi
    DEPS=
    if [[ ! -z "$AFTER" ]]; then
        DEPS="-d afterok:$AFTER"
    fi
    sbatch --export=ALL \
        -N "$NUM_NODES" -t "$WALLTIME" -p "$QUEUE" $DEPS \
        "$SOLEIL_DIR"/src/sapling.slurm
}

function run_local {
    if (( NUM_NODES > 1 )); then
        quit "Too many nodes requested"
    fi
    # Overrides for local, non-GPU run
    LOCAL_RUN=1
    USE_CUDA=0
    RESERVED_CORES=2
    # Synthesize final command
    CORES_PER_NODE="$(grep ^cpu\\scores /proc/cpuinfo | uniq |  awk '{print $4}')"
    RAM_PER_NODE="$(free -m | head -2 | tail -1 | awk '{print $2}')"
    RAM_PER_NODE=$(( RAM_PER_NODE / 2 ))
    source "$SOLEIL_DIR"/src/jobscript_shared.sh
    # Emit final command
    $COMMAND
}

###############################################################################
# Switch on machine
###############################################################################

if [[ "$(uname -n)" == *"titan"* ]]; then
    run_titan
elif [[ "$(hostname -d)" == *"summit"* ]]; then
    run_summit
elif [[ "$(uname -n)" == *"lassen"* || "$(uname -n)" == *"sierra"* ]]; then
    run_lassen
elif [[ "$(uname -n)" == *"daint"* ]]; then
    run_pizdaint
elif [[ "$(uname -n)" == *"certainty"* ]]; then
    run_certainty
elif [[ "$(uname -n)" == *"sapling"* ]]; then
    run_sapling
else
    echo 'Hostname not recognized; assuming local machine run w/o GPUs'
    run_local
fi
