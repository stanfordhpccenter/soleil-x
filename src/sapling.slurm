#!/bin/bash -eu
#SBATCH --exclusive
#SBATCH --job-name=soleil

source "$SOLEIL_DIR"/src/sapling_reload_modules.sh

CORES_PER_NODE=20
RAM_PER_NODE=200000
GPUS_PER_NODE=4
FB_PER_GPU=14000
source "$SOLEIL_DIR"/src/jobscript_shared.sh

mpirun -n "$NUM_RANKS" -npernode "$RANKS_PER_NODE" --bind-to none \
    -x LD_LIBRARY_PATH -x SOLEIL_DIR -x REALM_BACKTRACE -x LEGION_FREEZE_ON_ERROR -x DEBUG_COPYING \
    $COMMAND

# Resources:
# 256GB RAM per node
# 2 NUMA domains per node
# 10 cores per NUMA domain
# 4 Tesla P100 GPUs per node (4/8 nodes)
# 16GB FB per GPU
