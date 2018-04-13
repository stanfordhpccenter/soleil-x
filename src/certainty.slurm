#!/bin/bash -eu
#SBATCH -J soleil

GPUS=
if [[ "$QUEUE" == "32core" ]]; then
    NUM_THREADS=15
    # 64GB RAM per node
    # 2 NUMA domains per node
    # 16 cores per NUMA domain
    # 2-way SMT per core
elif [[ "$QUEUE" == "24core" ]]; then
    NUM_THREADS=11
    # 36GB RAM per node
    # 2 NUMA domains per node
    # 12 cores per NUMA domain
    # 2-way SMT per core
elif [[ "$QUEUE" == "largemem" ]]; then
    NUM_THREADS=11
    # 144GB RAM per node
    # 2 NUMA domains per node
    # 12 cores per NUMA domain
    # 2-way SMT per core
elif [[ "$QUEUE" == "gpu" ]]; then
    NUM_THREADS=11
    GPUS="-ll:gpu 4"
    # 36GB RAM per node
    # 2 NUMA domains per node
    # 8 cores per NUMA domain
    # 2-way SMT per core
    # 4 Kepler M2050 GPUs per node
    # 3GB FB per GPU
else
    echo "Unrecognized queue $QUEUE" >&2
    exit 1
fi

mpiexec -np "$SLURM_JOB_NUM_NODES" --bind-to none -x LD_LIBRARY_PATH \
    "$SOLEIL_DIR"/src/soleil.exec $ARGS \
    -ll:cpu 0 -ll:ocpu 1 -ll:onuma 0 -ll:okindhack -ll:othr "$NUM_THREADS" $GPUS \
    -ll:csize 28000 -ll:fsize 2500
