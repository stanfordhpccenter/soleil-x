Setup (generic)
===============

### Add to shell startup

Normally you'd need to edit file `~/.bashrc`.

```
# Module loads (if necessary)
...
# Legion build config (if necessary)
...
# Path setup (mandatory)
export LEGION_DIR=???
export LG_RT_DIR="$LEGION_DIR"/runtime
export SOLEIL_DIR=???
# HDF config (if using HDF, and the machine has an unusual HDF installation)
export HDF_HEADER=???
export HDF_LIBNAME=???
export HDF_ROOT=???
# CUDA config (if using CUDA code generation)
export CUDA_HOME=???
export CUDA="$CUDA_HOME"
export GPU_ARCH=???
```

### Download software

```
$ git clone -b master https://github.com/StanfordLegion/legion.git "$LEGION_DIR"
$ git clone https://github.com/stanfordhpccenter/soleil-x.git "$SOLEIL_DIR"
```

### Install Legion

```
$ cd "$LEGION_DIR"/language
$ unset LG_RT_DIR
$ USE_CUDA=? USE_OPENMP=? USE_GASNET=? USE_HDF=? scripts/setup_env.py --terra-url 'https://github.com/StanfordLegion/terra.git' --terra-branch 'luajit2.1'
```

or (if LLVM is installed system-wide, and you're not using GASNET):

```
$ cd "$LEGION_DIR"/language
$ git clone -b luajit2.1 https://github.com/StanfordLegion/terra.git terra.build
$ ln -s terra.build terra
$ cd "$LEGION_DIR"/language/terra
$ make all -j
$ cd "$LEGION_DIR"/language
$ USE_CUDA=? USE_OPENMP=? USE_GASNET=? USE_HDF=? ./install.py
```

See [Elliott's instructions](https://docs.google.com/document/d/1Qkl6r-1ZIb8WyH1f_WZbKgjp3due_Q8UiWKLh_nG1ec/edit) for more help.

Setup (local Ubuntu machine w/o GPU)
====================================

### Add to shell startup

```
# Path setup
export LEGION_DIR=???
export LG_RT_DIR="$LEGION_DIR"/runtime
export SOLEIL_DIR=???
# HDF config
export HDF_HEADER=hdf5/serial/hdf5.h
export HDF_LIBNAME=hdf5_serial
```

### Download software

```
$ sudo apt-get install llvm-3.8-dev clang-3.8
$ sudo apt-get install libhdf5-serial-dev
$ git clone -b master https://github.com/StanfordLegion/legion.git "$LEGION_DIR"
$ git clone https://github.com/stanfordhpccenter/soleil-x.git "$SOLEIL_DIR"
```

### Install Legion

```
$ cd "$LEGION_DIR"/language
$ git clone -b luajit2.1 https://github.com/StanfordLegion/terra.git terra.build
$ ln -s terra.build terra
$ cd "$LEGION_DIR"/language/terra
$ make all -j
$ cd "$LEGION_DIR"/language
$ USE_CUDA=0 USE_OPENMP=1 USE_GASNET=0 USE_HDF=1 ./install.py
```

Setup (Sapling @ Stanford)
==========================

### Add to shell startup

```
# Module loads
module load mpi/openmpi/1.8.2
module load gasnet/1.22.4-openmpi
module load cuda/8.0
# Legion build config
export CONDUIT=ibv
# Path setup
export LEGION_DIR=???
export LG_RT_DIR="$LEGION_DIR"/runtime
export SOLEIL_DIR=???
# CUDA config
export CUDA_HOME=/usr/local/cuda-8.0
export CUDA="$CUDA_HOME"
export GPU_ARCH=fermi
```

### Download software

```
$ git clone -b master https://github.com/StanfordLegion/legion.git "$LEGION_DIR"
$ git clone https://github.com/stanfordhpccenter/soleil-x.git "$SOLEIL_DIR"
```

### Install Legion

```
$ cd "$LEGION_DIR"/language
$ git clone -b luajit2.1 https://github.com/StanfordLegion/terra.git terra.build
$ ln -s terra.build terra
$ cd "$LEGION_DIR"/language/terra
$ make all -j
$ cd "$LEGION_DIR"/language
$ USE_CUDA=1 USE_OPENMP=1 USE_GASNET=1 USE_HDF=1 ./install.py
```

Setup (Certainty @ Stanford)
==========================

There's an issue currently on Certainty preventing the use of git over HTTPS. Until that is fixed, you will need to use git over SSH. This requires you to [set up an SSH key pair on github](https://help.github.com/articles/connecting-to-github-with-ssh/).

### Add to shell startup

```
# Module loads
module load gnu/5.2
module load openmpi/1.10.2-gnu-5.2
module load cuda/7.5.18
module load hdf5/1.10.1-gnu-5.2
module load python/2.7.8
# Legion build config
export CC=gcc
export CXX=g++
export CONDUIT=ibv
# Path setup
export LEGION_DIR=???
export LG_RT_DIR="$LEGION_DIR"/runtime
export SOLEIL_DIR=???
# CUDA config
export CUDA_HOME=/usr/local/cuda-7.5/
export CUDA="$CUDA_HOME"
export GPU_ARCH=fermi
```

### Download software

```
$ git clone -b master git@github.com:StanfordLegion/legion.git "$LEGION_DIR"
$ git clone git@github.com:stanfordhpccenter/soleil-x.git "$SOLEIL_DIR"
```

### Install Legion

```
$ cd "$LEGION_DIR"/language
$ unset LG_RT_DIR
$ git clone -b luajit2.1 git@github.com:StanfordLegion/terra.git terra.build
$ sed -i 's|https://github.com/|git@github.com:|g' terra.build/Makefile
$ sed -i 's|https://github.com/|git@github.com:|g' scripts/setup_env.py
$ sed -i 's|~eslaught|~manolis|g' scripts/setup_env.py
$ sed -i 's|85faf7cbd518dabeafc4d3f7e909338fc1dab3c4|c29dc69803e07c43c3654ace58e98368857e1669|g' scripts/setup_env.py
$ USE_CUDA=1 USE_OPENMP=1 USE_GASNET=1 USE_HDF=1 scripts/setup_env.py --llvm-version 35 --terra-url 'git@github.com:StanfordLegion/terra.git' --terra-branch 'luajit2.1'
```

Setup (Titan @ Oak Ridge)
=========================

### Add to shell startup

```
# Module loads
module load python/2.7.9
module load cray-hdf5/1.10.0.1
module load cudatoolkit/7.5.18-1.0502.10743.2.1
module swap PrgEnv-pgi PrgEnv-gnu
# Legion build config
export MARCH=barcelona
export CC=cc
export CXX=CC
export HOST_CC=gcc
export HOST_CXX=g++
# Path setup
export LEGION_DIR=???
export LG_RT_DIR="$LEGION_DIR"/runtime
export SOLEIL_DIR=???
# CUDA config
export CUDA_HOME=/opt/nvidia/cudatoolkit7.5/7.5.18-1.0502.10743.2.1
export CUDA="$CUDA_HOME"
export GPU_ARCH=kepler
```

### Download software

```
$ git clone -b master https://github.com/StanfordLegion/legion.git "$LEGION_DIR"
$ git clone https://github.com/stanfordhpccenter/soleil-x.git "$SOLEIL_DIR"
```

### Install Legion

```
$ cd "$LEGION_DIR"/language
$ unset LG_RT_DIR
$ USE_CUDA=1 USE_OPENMP=1 USE_GASNET=1 USE_HDF=1 scripts/setup_env.py --terra-url 'https://github.com/StanfordLegion/terra.git' --terra-branch 'luajit2.1'
```

Running
=======

```
$ cd "$SOLEIL_DIR"/src
$ ./compile.sh
$ ./run.sh
```

Edit the `compile.sh` and `run.sh` scripts to choose a different testcase, edit runtime options etc.
