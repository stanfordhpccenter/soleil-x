#!/bin/bash
LD_FLAGS=-lpmi2 USE_CUDA=1 USE_OPENMP=1 USE_GASNET=1 USE_HDF=1 srun -N 1 -c 10 -p aaiken --gres=gpu:4 scripts/setup_env.py --llvm-version 38

