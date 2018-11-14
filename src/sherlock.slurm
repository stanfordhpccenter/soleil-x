#!/bin/bash -eu
#SBATCH --exclusive
#SBATCH --job-name=soleil
#SBATCH --partition=aaiken

source "$SOLEIL_DIR"/src/jobscript_shared.sh

# Process CUDA configuration
GPU_OPTS=
if [[ "$USE_CUDA" == 1 ]]; then
    GPU_OPTS="-ll:gpu 4 -ll:fsize 15000 -ll:zsize 1024 -ll:ib_zsize 0"
fi

srun -N "$NUM_RANKS" --ntasks-per-node=1 --cpus-per-task=20 \
    --export=ALL \
    "$EXECUTABLE" $ARGS \
    -ll:cpu 0 -ll:ocpu 1 -ll:onuma 0 -ll:okindhack -ll:othr 19 \
    $GPU_OPTS -ll:dma 2 -ll:ahandlers 2 \
    -ll:csize 240000 -ll:rsize 1024 -ll:ib_rsize 1024 -ll:gsize 0 \
    -ll:stacksize 8 -ll:ostack 8 -lg:sched -1

# Resources:
# 256GB RAM per node
# 2 NUMA domains per node
# 10 cores per NUMA domain
# 4 Tesla P100 GPUs per node (4/8 nodes)
# 16GB FB per GPU
