#!/bin/bash -eu

SOLEIL_SRC="$(cd "$(dirname "$(perl -MCwd -le 'print Cwd::abs_path(shift)' "${BASH_SOURCE[0]}")")" && pwd)"
cd "$SOLEIL_SRC"

# Soleil options
export HDF_HEADER="${HDF_HEADER:-hdf5.h}"
export HDF_LIBNAME="${HDF_LIBNAME:-hdf5}"
export USE_HDF=0

# Regent options
export INCLUDE_PATH="."
export LIBRARY_PATH="."
if [ ! -z "${HDF_ROOT:-}" ]; then
    export INCLUDE_PATH="$INCLUDE_PATH;$HDF_ROOT/include"
    export LIBRARY_PATH="$LIBRARY_PATH:$HDF_ROOT/lib"
fi
export REGENT_FLAGS="-fflow 0 -fopenmp 1 -fcuda 1 -fcuda-offline 1"

# Build libraries
gcc -g -O2 -c -o json.o json.c
"$LEGION_DIR"/language/regent.py config_helper.rg config_schema.lua

# Compile Soleil
"$LEGION_DIR"/language/regent.py soleil.rg $REGENT_FLAGS
