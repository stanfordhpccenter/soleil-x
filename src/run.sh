#!/bin/bash -eu

###############################################################################
# Derived options
###############################################################################

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

###############################################################################
# Machine-specific handling
###############################################################################

function run_titan {
    export QUEUE="${QUEUE:-debug}"
    qsub -V \
        -l nodes="$NUM_RANKS" -l walltime="$WALLTIME" -q "$QUEUE" \
        "$SOLEIL_DIR"/src/titan.pbs
}

function run_summit {
    export QUEUE="${QUEUE:-batch}"
    EXCLUDED="$(cat "$SOLEIL_DIR"/src/blacklist/summit.txt |
                sed 's/^/ \&\& (hname != /'  | sed 's/$/)/' |
                paste -sd '')"
    NUM_NODES="$(( NUM_RANKS/2 + NUM_RANKS%2 ))"
    NUM_CORES="$(( NUM_NODES * 42 ))"
    RESOURCES="1*{select[LN]span[hosts=1]} +
               $NUM_CORES*{select[CN$EXCLUDED]span[ptile=42]}"
    bsub -csm y -J soleil -P CSC275IACCARINO -alloc_flags smt4 \
        -R "$RESOURCES" -W "$MINUTES" -q "$QUEUE" \
        "$SOLEIL_DIR"/src/summit.lsf
}

function run_pizdaint {
    export QUEUE="${QUEUE:-debug}"
    sbatch --export=ALL \
        -N "$NUM_RANKS" -t "$WALLTIME" -p "$QUEUE" \
        "$SOLEIL_DIR"/src/pizdaint.slurm
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

function run_sherlock {
    RESOURCES=
    if [[ "$USE_CUDA" == 1 ]]; then
        RESOURCES="gpu:4"
    fi
    sbatch --export=ALL \
        -N "$NUM_RANKS" -t "$WALLTIME" --gres="$RESOURCES" \
        "$SOLEIL_DIR"/src/sherlock.slurm
}

function run_sapling {
    # Allocate up to 3 nodes, from n0001 up to n0003
    if (( NUM_RANKS > 3 )); then quit "Too many nodes requested"; fi
    NODES=n0001
    for (( i = 2; i <= NUM_RANKS; i++ )); do
        NODES="$NODES,n000$i"
    done
    GPU_OPTS=
    if [[ "$USE_CUDA" == 1 ]]; then
        GPU_OPTS="-ll:gpu 1 -ll:fsize 5500"
    fi
    mpiexec -H "$NODES" --bind-to none \
        -x LD_LIBRARY_PATH -x SOLEIL_DIR -x REALM_BACKTRACE \
        "$EXECUTABLE" $ARGS \
        -ll:cpu 0 -ll:ocpu 1 -ll:onuma 0 -ll:okindhack -ll:othr 8 \
        $GPU_OPTS -ll:util 2 -ll:dma 2 -ll:ahandlers 2 \
        -ll:csize 36000 -ll:ostack 8 -lg:sched -1
    # Resources:
    # 40230MB RAM per node
    # 2 NUMA domains per node
    # 4 cores per NUMA domain
    # 2 Tesla C2070 GPUs per node
    # 6GB FB per GPU
}

function run_local {
    "$EXECUTABLE" $ARGS \
        -ll:cpu 0 -ll:ocpu 1 -ll:onuma 0 -ll:okindhack -ll:othr 3 \
        -ll:csize 9000 -ll:ostack 8 -lg:sched -1
}

###############################################################################
# Switch on machine
###############################################################################

if [[ "$(uname -n)" == *"titan"* ]]; then
    run_titan
elif [[ "$(hostname -d)" == *"summit"* ]]; then
    run_summit
elif [[ "$(uname -n)" == *"daint"* ]]; then
    run_pizdaint
elif [[ "$(uname -n)" == *"certainty"* ]]; then
    run_certainty
elif [[ "$(uname -n)" == *"sh-ln"* ]]; then
    run_sherlock
elif [[ "$(uname -n)" == *"sapling"* ]]; then
    run_sapling
else
    run_local
fi
