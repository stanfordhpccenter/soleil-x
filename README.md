Setup (generic)
===============

### Add to shell startup

Normally you'd need to edit file `~/.bashrc`.

```
# Any required module loads:
...
# Mandatory:
export LEGION_DIR=???
export LG_RT_DIR="$LEGION_DIR"/runtime
export SOLEIL_DIR=???
# Optional, if using HDF and the machine has an unusual HDF installation:
export HDF_HEADER=???
export HDF_LIBNAME=???
export HDF_ROOT=???
# Optional, if using CUDA code generation:
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
$ [CONDUIT=???] USE_CUDA=? USE_OPENMP=? USE_GASNET=? USE_HDF=? scripts/setup_env.py --terra-url 'https://github.com/magnatelee/terra.git' --terra-branch 'luajit2.1'
```

or (if not using GASNET, or it's already installed system-wide):

```
$ cd "$LEGION_DIR"/language
$ git clone -b luajit2.1 https://github.com/magnatelee/terra.git terra.build
$ ln -s terra.build terra
$ cd "$LEGION_DIR"/language/terra
$ make all -j
$ cd "$LEGION_DIR"/language
$ [CONDUIT=???] USE_CUDA=? USE_OPENMP=? USE_GASNET=? USE_HDF=? ./install.py
```

See [Elliott's instructions](https://docs.google.com/document/d/1Qkl6r-1ZIb8WyH1f_WZbKgjp3due_Q8UiWKLh_nG1ec/edit) for more help.

Setup (local Ubuntu machine w/o GPU)
====================================

### Add to shell startup

```
export LEGION_DIR=???
export LG_RT_DIR="$LEGION_DIR"/runtime
export SOLEIL_DIR=???
export HDF_HEADER=hdf5/serial/hdf5.h
export HDF_LIBNAME=hdf5_serial
```

### Download software

```
$ sudo apt-get install libhdf5-serial-dev
$ git clone -b master https://github.com/StanfordLegion/legion.git "$LEGION_DIR"
$ git clone https://github.com/stanfordhpccenter/soleil-x.git "$SOLEIL_DIR"
```

### Install Legion

```
$ cd "$LEGION_DIR"/language
$ git clone -b luajit2.1 https://github.com/magnatelee/terra.git terra.build
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
module load mpi/openmpi/1.8.2
module load gasnet/1.22.4-openmpi
module load cuda/8.0
export LEGION_DIR=???
export LG_RT_DIR="$LEGION_DIR"/runtime
export SOLEIL_DIR=???
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
$ git clone -b luajit2.1 https://github.com/magnatelee/terra.git terra.build
$ ln -s terra.build terra
$ cd "$LEGION_DIR"/language/terra
$ make all -j
$ cd "$LEGION_DIR"/language
$ CONDUIT=ibv USE_CUDA=1 USE_OPENMP=1 USE_GASNET=1 USE_HDF=1 ./install.py
```

Setup (Titan @ Oak Ridge)
=========================



### Add to shell startup

```
module load python/2.7.9
module load cray-hdf5/1.10.0.1
module load cudatoolkit/7.5.18-1.0502.10743.2.1
module swap PrgEnv-pgi PrgEnv-gnu
export LEGION_DIR=???
export LG_RT_DIR="$LEGION_DIR"/runtime
export SOLEIL_DIR=???
export MARCH=barcelona
export CC=cc
export CXX=CC
export HOST_CC=gcc
export HOST_CXX=g++
export HDF_ROOT=/opt/cray/hdf5/1.10.0.1/GNU/4.9
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
$ USE_CUDA=1 USE_OPENMP=1 USE_GASNET=1 USE_HDF=1 scripts/setup_env.py --terra-url 'https://github.com/magnatelee/terra.git' --terra-branch 'luajit2.1'
```

Running
=======

```
$ cd "$SOLEIL_DIR"/src
$ ./compile.sh
$ ./run.sh
```

Edit the `compile.sh` and `run.sh` scripts to choose a different testcase, edit runtime options etc.
