#!/bin/bash -eu
#SBATCH --exclusive
#SBATCH --job-name=soleil
#SBATCH --partition=aaiken

CORES_PER_NODE=20
RAM_PER_NODE=200000
GPUS_PER_NODE=4
FB_PER_GPU=15000
source "$SOLEIL_DIR"/src/jobscript_shared.sh

srun -n "$NUM_RANKS" --ntasks-per-node="$RANKS_PER_NODE" --cpus-per-task="$CORES_PER_RANK" \
    --export=ALL \
    $COMMAND

# Resources:
# 256GB RAM per node
# 2 NUMA domains per node
# 10 cores per NUMA domain
# 4 Tesla P100 GPUs per node (4/8 nodes)
# 16GB FB per GPU
