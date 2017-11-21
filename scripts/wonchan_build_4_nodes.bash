#!/bin/bash
export TERRA_PATH=${LISZT_PATH}/include/?.t
export OUTPUT_PATH=.
echo terra path ${TERRA_PATH}
export USE_HDF=0
export SAVEOBJ=1
export OBJNAME=${OUTPUT_PATH}/taylor-512x512x256-shardsize1
COMMAND="${LEGION_PATH}/language/regent.py ${SOLEIL_PATH}/src/soleil-x.t -i ${SOLEIL_PATH}/testcases/taylor_green_vortex/taylor_green_vortex_512_512_256.lua    -fcuda 0 -fopenmp 1 -fflow-spmd 1 -fflow-spmd-shardsize 1 -fparallelize-dop 2,2,1 "
echo $COMMAND
$COMMAND
