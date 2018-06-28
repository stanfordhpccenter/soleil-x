#!/bin/bash -eu
#SBATCH --exclusive
#SBATCH --job-name=soleil
#SBATCH --partition=aaiken

GPU_OPTS=
if [[ "$USE_CUDA" == 1 ]]; then
    GPU_OPTS="-ll:gpu 4 -ll:fsize 15000"
fi

srun -N "$NUM_RANKS" --ntasks-per-node=1 --cpus-per-task=20 \
    --export=ALL \
    "$SOLEIL_DIR"/src/soleil.exec $ARGS $GPU_OPTS \
    -ll:cpu 0 -ll:ocpu 1 -ll:onuma 0 -ll:okindhack -ll:othr 19 \
    -ll:csize 240000 -ll:ostack 8

# Resources:
# 256GB RAM per node
# 2 NUMA domains per node
# 10 cores per NUMA domain
# 4 Tesla P100 GPUs per node (4/8 nodes)
# 16GB FB per GPU
