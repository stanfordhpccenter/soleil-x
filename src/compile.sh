#!/bin/bash -eu

SOLEIL_SRC="$(cd "$(dirname "$(readlink -f "${BASH_SOURCE[0]}")")" && pwd)"

cd "$SOLEIL_SRC"

gcc -g -O2 -c json.c
ar rcs libjsonparser.a json.o

USE_HDF=1 HDF_LIBNAME=hdf5_serial HDF_HEADER=hdf5/serial/hdf5.h \
       LIBRARY_PATH="." \
       INCLUDE_PATH="." \
       DEBUG=1 OBJNAME=soleil.exec \
       $LISZT_PATH/liszt-legion.sh soleil-x.t \
       -i ../testcases/cavity/cavity_32x32.json \
       -fdebug 1 1> soleil.out 2> soleil.err

cd -
