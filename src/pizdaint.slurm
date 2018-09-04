#!/bin/bash -eu
#SBATCH --job-name=soleil
#SBATCH --constraint=gpu

source "$SOLEIL_DIR"/src/jobscript_shared.sh

GPU_OPTS=
if [[ "$USE_CUDA" == 1 ]]; then
    GPU_OPTS="-ll:gpu 1 -ll:fsize 15000"
fi

srun -N "$NUM_RANKS" -n "$NUM_RANKS" \
    --export=ALL \
    "$EXECUTABLE" $ARGS \
    -ll:cpu 0 -ll:ocpu 1 -ll:onuma 0 -ll:okindhack -ll:othr 11 \
    $GPU_OPTS -ll:dma 2 -ll:ahandlers 2 \
    -ll:csize 58000 -ll:ostack 8 -lg:sched -1

# Resources:
# 12 cores per node
# 2-way SMT per core
# 64GB RAM per node
# 1 Tesla P100 GPU per node
# 16GB FB per GPU
