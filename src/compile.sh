#!/bin/bash -eu

SOLEIL_SRC="$(cd "$(dirname "$(perl -MCwd -le 'print Cwd::abs_path(shift)' "${BASH_SOURCE[0]}")")" && pwd)"
cd "$SOLEIL_SRC"

# Translator options
export HDF_HEADER="${HDF_HEADER:-hdf5.h}"
export HDF_LIBNAME="${HDF_LIBNAME:-hdf5}"
export USE_HDF=0
export DUMP_REGENT=1
export OBJNAME=soleil.exec

# Regent options
export INCLUDE_PATH="."
export LIBRARY_PATH="."
if [ ! -z "${HDF_ROOT:-}" ]; then
    export INCLUDE_PATH="$INCLUDE_PATH;$HDF_ROOT/include"
    export LIBRARY_PATH="$LIBRARY_PATH:$HDF_ROOT/lib"
fi
export TERRA_PATH=liszt/?.t

# Build libraries
gcc -g -O2 -c -o json.o json.c
ar rcs libjsonparser.a json.o

# Compile Liszt
"$LEGION_DIR"/language/regent.py soleil-x.t \
    -i ../testcases/tgv_64x64x64.json \
    -fflow 0 -fflow-spmd 0 -fcuda 0 -fopenmp 1 \
    1> soleil.out 2> soleil.err

# Post-process dumped Regent
if [[ $DUMP_REGENT == 1 ]]; then
    ./make_parsable.py soleil.out > soleil.rg
fi
