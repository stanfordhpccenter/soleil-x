#!/bin/bash -eu

SOLEIL_SRC="$(cd "$(dirname "$(perl -MCwd -le 'print Cwd::abs_path(shift)' "${BASH_SOURCE[0]}")")" && pwd)"
cd "$SOLEIL_SRC"

# Translator options
if [[ $(uname -n) == *"titan"* ]]; then
    export HDF_LIBNAME=hdf5
    export HDF_HEADER=hdf5.h
else
    export HDF_LIBNAME=hdf5_serial
    export HDF_HEADER=hdf5/serial/hdf5.h
fi
export USE_HDF=1
export DUMP_REGENT=1
export OBJNAME=soleil.exec

# Regent options
if [[ $(uname -n) == *"titan"* ]]; then
    export LIBRARY_PATH=".:/opt/cray/hdf5/1.10.0.1/GNU/4.9/lib"
    export INCLUDE_PATH=".;/opt/cray/hdf5/1.10.0.1/GNU/4.9/include"
else
    export LIBRARY_PATH="."
    export INCLUDE_PATH="."
fi
export TERRA_PATH=liszt/?.t

# Build libraries
gcc -g -O2 -c -o json.o json.c
ar rcs libjsonparser.a json.o

# Compile Liszt
"$LEGION_DIR"/language/regent.py soleil-x.t \
    -i ../testcases/cavity/cavity_32x32.json \
    -fflow 0 -fflow-spmd 0 -fcuda 0 -fopenmp 0 \
    1> soleil.out 2> soleil.err

# Post-process dumped Regent
if [[ $DUMP_REGENT == 1 ]]; then
    ./make_parsable.py soleil.out > soleil.rg
fi
