#!/bin/bash
mkdir LLVM
cd LLVM
wget http://releases.llvm.org/4.0.0/llvm-4.0.0.src.tar.xz
tar xJf llvm-4.0.0.src.tar.xz
mkdir build
cd build
module load cmake
cmake -DLLVM_BUILD_LLVM_DYLIB=1 -DLLVM_TARGETS_TO_BUILD=X86 ../llvm-4.0.0.src
make
cmake -DCMAKE_INSTALL_PREFIX=$LLVM_DIR/install -P cmake_install.cmake
