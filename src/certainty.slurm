#!/bin/bash -eu
#SBATCH -J soleil

source "$SOLEIL_DIR"/src/jobscript_shared.sh

GPU_OPTS=
if [[ "$QUEUE" == "all" ]]; then
    NUM_THREADS=11
elif [[ "$QUEUE" == "32core" ]]; then
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
    if [[ "$USE_CUDA" == 1 ]]; then
        GPU_OPTS="-ll:gpu 4 -ll:fsize 2500 -ll:zsize 1024 -ll:ib_zsize 0"
    fi
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

mpiexec -np "$NUM_RANKS" --map-by ppr:1:node --bind-to none \
    -x LD_LIBRARY_PATH -x SOLEIL_DIR -x REALM_BACKTRACE \
    "$EXECUTABLE" $ARGS \
    -ll:cpu 0 -ll:ocpu 1 -ll:onuma 0 -ll:okindhack -ll:othr "$NUM_THREADS" \
    $GPU_OPTS -ll:dma 2 -ll:ahandlers 2 \
    -ll:csize 28000 -ll:rsize 1024 -ll:ib_rsize 1024 -ll:gsize 0 \
    -ll:stacksize 8 -ll:ostack 8 -lg:sched -1
