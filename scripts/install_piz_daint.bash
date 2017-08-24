#!/bin/bash
module unload PrgEnv-cray
module load PrgEnv-gnu
module load daint-gpu
unset LG_RT_DIR

cd ${LEGION_PATH}/language
rm -rf terra.build terra terra.build.master
USE_OPENMP=1 CC=cc CXX=CC HOST_CC=gcc HOST_CXX=g++ scripts/setup_env.py 
rm -rf terra.build.master
mv terra.build terra.build.master
git clone -b luajit2.1 https://github.com/magnatelee/terra.git terra.build
rm -rf terra
ln -s terra.build terra
cd terra.build
CC=gcc CXX=g++ make LLVM_CONFIG=`readlink -f ../llvm/install/bin/llvm-config` CLANG=`readlink -f ../llvm/install/bin/clang` -j
cd ..


