#!/bin/bash -eu
#SBATCH -J soleil

USE_CUDA=0
if [[ "$QUEUE" == "all" ]]; then
    CORES_PER_NODE=24
    RAM_PER_NODE=25000
elif [[ "$QUEUE" == "32core" ]]; then
    CORES_PER_NODE=32
    RAM_PER_NODE=50000
    # 64GB RAM per node
    # 2 NUMA domains per node
    # 16 cores per NUMA domain
    # 2-way SMT per core
elif [[ "$QUEUE" == "24core" ]]; then
    CORES_PER_NODE=24
    RAM_PER_NODE=25000
    # 36GB RAM per node
    # 2 NUMA domains per node
    # 12 cores per NUMA domain
    # 2-way SMT per core
elif [[ "$QUEUE" == "largemem" ]]; then
    CORES_PER_NODE=24
    RAM_PER_NODE=120000
    # 144GB RAM per node
    # 2 NUMA domains per node
    # 12 cores per NUMA domain
    # 2-way SMT per core
else
    echo "Unrecognized queue $QUEUE" >&2
    exit 1
fi
source "$SOLEIL_DIR"/src/jobscript_shared.sh

mpiexec -np "$NUM_RANKS" --map-by ppr:"$RANKS_PER_NODE":node --bind-to none \
    -x LD_LIBRARY_PATH -x SOLEIL_DIR -x REALM_BACKTRACE -x LEGION_FREEZE_ON_ERROR -x DEBUG_COPYING \
    $COMMAND
