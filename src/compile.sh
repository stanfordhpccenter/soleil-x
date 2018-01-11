#!/bin/bash -eu

SOLEIL_SRC="$(cd "$(dirname "$(perl -MCwd -le 'print Cwd::abs_path(shift)' "${BASH_SOURCE[0]}")")" && pwd)"
cd "$SOLEIL_SRC"

# Translation options
export HDF_LIBNAME=hdf5_serial
export HDF_HEADER=hdf5/serial/hdf5.h
export USE_HDF=1
export DEBUG=1
export OBJNAME=soleil.exec

# Regent options
export TERRA_PATH=liszt/?.t
export LIBRARY_PATH="${LIBRARY_PATH:-}:."
export INCLUDE_PATH="."

gcc -g -O2 -c -o json.o json.c
ar rcs libjsonparser.a json.o

"$LEGION_PATH"/language/regent.py soleil-x.t \
    -i ../testcases/cavity/cavity_32x32.json \
    -fflow 0 -fflow-spmd 0 -fcuda 0 -fopenmp 0 \
    -fdebug 1 1> soleil.out 2> soleil.err

./make_parsable.py soleil.out > soleil.rg
