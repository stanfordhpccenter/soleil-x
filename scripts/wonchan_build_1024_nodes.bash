#!/bin/bash
time USE_HDF=0 TERRA_PATH=liszt-legion/include/?.t SAVEOBJ=1 OBJNAME="$OUTPUT_PATH"/taylor-4096x2048x2048-shardsize1 ./regent.py soleil-master/src/soleil-x.t -i soleil-master/testcases/taylor_green_vortex/taylor_green_vortex_4096_2048_2048.lua -fcuda 0 -fopenmp 1 -fflow-spmd 1 -fflow-spmd-shardsize 1 -fparallelize-dop 16,8,8
