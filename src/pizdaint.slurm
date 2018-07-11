#!/bin/bash -eu
#SBATCH --job-name=soleil
#SBATCH --constraint=gpu

# Create a directory on the scratch filesystem, for output (if the user hasn't
# specified an output directory explicitly).
OUT_DIR=
OUT_DIR_FOLLOWS=false
for ARG in $ARGS; do
    if [[ "$OUT_DIR_FOLLOWS" == true ]]; then
        OUT_DIR="$ARG"
        break
    elif [[ "$ARG" == "-o" ]]; then
        OUT_DIR_FOLLOWS=true
    fi
done
if [[ -z "$OUT_DIR" ]]; then
    OUT_DIR="$SCRATCH/$SLURM_JOBID"
    mkdir "$OUT_DIR"
    ARGS="$ARGS -o $OUT_DIR"
    echo "Redirecting output to $OUT_DIR"
fi

GPU_OPTS=
if [[ "$USE_CUDA" == 1 ]]; then
    GPU_OPTS="-ll:gpu 1 -ll:fsize 15000"
fi

srun -N "$NUM_RANKS" -n "$NUM_RANKS" \
    --export=ALL \
    "$SOLEIL_DIR"/src/soleil.exec $ARGS $GPU_OPTS \
    -ll:cpu 0 -ll:ocpu 1 -ll:onuma 0 -ll:okindhack -ll:othr 11 \
    -ll:csize 58000 -ll:ostack 8

# Resources:
# 12 cores per node
# 2-way SMT per core
# 64GB RAM per node
# 1 Tesla P100 GPU per node
# 16GB FB per GPU
