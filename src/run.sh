#!/bin/bash -eu

###############################################################################
# Helper functions
###############################################################################

function quit {
    echo "$1" >&2
    exit 1
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
    GROUP="${GROUP:-CSC275IACCARINO}"
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
    bsub -csm y -J soleil -P "$GROUP" -alloc_flags smt4 \
        -R "$RESOURCES" -W "$MINUTES" -q "$QUEUE" $DEPS \
        "$SOLEIL_DIR"/src/summit.lsf
}

function run_lassen {
    GROUP="${GROUP:-guests}"
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
    export QUEUE="${QUEUE:-gpu}"
    RESOURCES=
    if [[ "$QUEUE" == "gpu" ]]; then
	RESOURCES="gpu:4"
    fi
    EXCLUDED="$(paste -sd ',' "$SOLEIL_DIR"/src/blacklist/certainty.txt)"
    DEPS=
    if [[ ! -z "$AFTER" ]]; then
        DEPS="-d afterok:$AFTER"
    fi
    sbatch --export=ALL \
        -N "$NUM_NODES" -t "$WALLTIME" -p "$QUEUE" --gres="$RESOURCES" $DEPS \
	--exclude="$EXCLUDED" \
        "$SOLEIL_DIR"/src/certainty.slurm
}

function run_sherlock {
    RESOURCES=
    if [[ "$USE_CUDA" == 1 ]]; then
        RESOURCES="gpu:4"
    fi
    DEPS=
    if [[ ! -z "$AFTER" ]]; then
        DEPS="-d afterok:$AFTER"
    fi
    sbatch --export=ALL \
        -N "$NUM_NODES" -t "$WALLTIME" --gres="$RESOURCES" $DEPS \
        "$SOLEIL_DIR"/src/sherlock.slurm
}

function run_sapling {
    # Allocate up to 2 nodes, from n0002 up to n0003
    if (( NUM_NODES == 1 )); then
        NODES="n0002"
    elif (( NUM_NODES == 2 )); then
        NODES="n0002,n0003"
    else
        quit "Too many nodes requested"
    fi
    # Synthesize final command
    CORES_PER_NODE=12
    RAM_PER_NODE=30000
    GPUS_PER_NODE=2
    FB_PER_GPU=5000
    source "$SOLEIL_DIR"/src/jobscript_shared.sh
    # Emit final command
    mpiexec -H "$NODES" --bind-to none \
        -x LD_LIBRARY_PATH -x SOLEIL_DIR -x GASNET_BACKTRACE -x SCRATCH \
        $COMMAND
    # Resources:
    # 40230MB RAM per node
    # 2 NUMA domains per node
    # 6 cores per NUMA domain
    # 2-way SMT per core
    # 2 Tesla C2070 GPUs per node
    # 6GB FB per GPU
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
    CORES_PER_NODE=4
    RAM_PER_NODE="$(free -m | head -2 | tail -1 | awk '{print $2}')"
    RAM_PER_NODE=$(( RAM_PER_NODE / 2 ))
    source "$SOLEIL_DIR"/src/jobscript_shared.sh
    # Emit final command
    LEGION_FREEZE_ON_ERROR=1 $COMMAND
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
elif [[ "$(uname -n)" == *"sh-ln"* ]]; then
    run_sherlock
elif [[ "$(uname -n)" == *"sapling"* ]]; then
    run_sapling
else
    echo 'Hostname not recognized; assuming local machine run w/o GPUs'
    run_local
fi
