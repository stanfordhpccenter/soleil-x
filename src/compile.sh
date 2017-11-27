#!/bin/bash -eu

SOLEIL_SRC="$(cd "$(dirname "$(perl -MCwd -le 'print Cwd::abs_path(shift)' "${BASH_SOURCE[0]}")")" && pwd)"

export HDF_LIBNAME=hdf5_serial
export HDF_HEADER=hdf5/serial/hdf5.h

cd "$SOLEIL_SRC"

USE_HDF=1 DEBUG=1 OBJNAME=soleil.exec \
    "$LISZT_PATH"/liszt-legion.sh soleil-x.t \
    -i ../testcases/cavity/cavity_32x32.json \
    1> soleil.out 2> soleil.err

"$LISZT_PATH"/make_parsable.py soleil.out > soleil.rg
