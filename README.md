Setup (generic)
===============

See below for instructions targeting specific systems.

### Prerequisites

* Legion (latest version -- Soleil-X follows the developments in Legion)
* GCC 4.9+ (we need a working `std::regex` library)
* CUDA 7.5-9.0 (Regent's CUDA codegen won't work with later versions)
* Python 2.X

The following are automatically installed during Legion installation:

* LLVM 3.8 (for CUDA 8.0-9.0) or 3.5 (for CUDA 7.5)
* GASNET (custom version)
* Terra (custom version -- we need to use LuaJIT2.1 instead of LuaJIT2.0, because the latter exhibits a spurious out-of-memory error when compiling large Regent programs)
* HDF5 (any recent version)

### Add to shell startup

Normally you'd need to edit file `~/.bashrc`. Replace the `???` depending on your system.

```
# Module loads (if necessary)
...
# Build config (if necessary, for Legion or Soleil)
...
# Path setup (mandatory)
export LEGION_DIR=???
export HDF_ROOT="$LEGION_DIR"/language/hdf/install
export SOLEIL_DIR=???
# CUDA config (if using CUDA code generation)
export CUDA_HOME=???
export CUDA="$CUDA_HOME"
export GPU_ARCH=???
```

### Download software

```
git clone https://gitlab.com/StanfordLegion/legion.git "$LEGION_DIR"
git clone https://github.com/stanfordhpccenter/soleil-x.git "$SOLEIL_DIR"
```

### Install Legion

Replace the `?` depending on your system's capabilities.

```
cd "$LEGION_DIR"/language
USE_CUDA=? USE_OPENMP=? USE_GASNET=? USE_HDF=? scripts/setup_env.py --llvm-version 38 --terra-url 'https://github.com/StanfordLegion/terra.git' --terra-branch 'luajit2.1'
```

See [Elliott's instructions](https://docs.google.com/document/d/1Qkl6r-1ZIb8WyH1f_WZbKgjp3due_Q8UiWKLh_nG1ec/edit) for more help.

Setup (local Ubuntu machine w/o GPU)
====================================

### Add to shell startup

```
# Build config
export CC=gcc
export CXX=g++
# Path setup
export LEGION_DIR=???
export HDF_ROOT="$LEGION_DIR"/language/hdf/install
export SOLEIL_DIR=???
```

### Download software

```
git clone https://gitlab.com/StanfordLegion/legion.git "$LEGION_DIR"
git clone https://github.com/stanfordhpccenter/soleil-x.git "$SOLEIL_DIR"
```

### Install Legion

```
cd "$LEGION_DIR"/language
USE_CUDA=0 USE_OPENMP=1 USE_GASNET=0 USE_HDF=1 scripts/setup_env.py --llvm-version 38 --terra-url 'https://github.com/StanfordLegion/terra.git' --terra-branch 'luajit2.1'
```

Setup (Sapling @ Stanford)
==========================

### Add to shell startup

```
# Module loads
module load mpi/openmpi/1.8.2
module load cuda/8.0
# Build config
export CONDUIT=ibv
export CC=gcc-4.9
export CXX=g++-4.9
# Path setup
export LEGION_DIR=???
export HDF_ROOT="$LEGION_DIR"/language/hdf/install
export SOLEIL_DIR=???
# CUDA config
export CUDA_HOME=/usr/local/cuda-8.0
export CUDA="$CUDA_HOME"
export GPU_ARCH=fermi
```

### Download software

```
git clone https://gitlab.com/StanfordLegion/legion.git "$LEGION_DIR"
git clone https://github.com/stanfordhpccenter/soleil-x.git "$SOLEIL_DIR"
```

### Install Legion

```
cd "$LEGION_DIR"/language
USE_CUDA=1 USE_OPENMP=1 USE_GASNET=1 USE_HDF=1 scripts/setup_env.py --llvm-version 38 --terra-url 'https://github.com/StanfordLegion/terra.git' --terra-branch 'luajit2.1'
```

Setup (Certainty @ Stanford)
============================

### Add to shell startup

```
# Module loads
module load gnu7/7.2.0
module load cuda/8.0
module load openmpi3/3.0.0
# Build config
export CONDUIT=ibv
export CC=gcc
export CXX=g++
# Path setup
export LEGION_DIR=???
export HDF_ROOT="$LEGION_DIR"/language/hdf/install
export SOLEIL_DIR=???
# CUDA config
export CUDA_HOME=/usr/local/cuda-8.0
export CUDA="$CUDA_HOME"
export GPU_ARCH=fermi
```

### Download software

```
git clone https://gitlab.com/StanfordLegion/legion.git "$LEGION_DIR"
git clone https://github.com/stanfordhpccenter/soleil-x.git "$SOLEIL_DIR"
```

### Install Legion

Legion's CUDA codegen is currently not working on Certainty, so we have to build without CUDA support.

```
# Disable PMI in GASnet, because the PMI library is missing on Certainty.
git clone https://github.com/StanfordLegion/gasnet.git $LEGION_DIR/language/gasnet
cd "$LEGION_DIR"/language/gasnet
sed -i 's|../$(GASNET_VERSION)/configure|../$(GASNET_VERSION)/configure --disable-pmi|g' Makefile
make
# Rest of compilation as normal
cd "$LEGION_DIR"/language
USE_CUDA=0 USE_OPENMP=1 USE_GASNET=1 USE_HDF=1 scripts/setup_env.py --llvm-version 38 --terra-url 'https://github.com/StanfordLegion/terra.git' --terra-branch 'luajit2.1'
```

Setup (Titan @ ORNL)
====================

### Add to shell startup

Install Legion and Soleil-X on the `/ccs/proj` filesystem, not your home directory (home directories are not mounted on the compute nodes).

We set `MARCH`, the processor architecture that Legion will be built for, to one that is compatible with both the login node and the compute nodes. This is done so that `libregent.so` can be used by both the final executable, when it's executing on a compute node, and the Regent compiler, which needs to link against it while compiling Soleil on the login node.

```
# Module loads
module load python/2.7.9
module load cudatoolkit/7.5.18-1.0502.10743.2.1
module swap PrgEnv-pgi PrgEnv-gnu
# Build config
export MARCH=barcelona
export CC=cc
export CXX=CC
export HOST_CC=gcc
export HOST_CXX=g++
# Path setup
export LEGION_DIR=???
export HDF_ROOT="$LEGION_DIR"/language/hdf/install
export SOLEIL_DIR=???
# CUDA config
export CUDA_HOME=/opt/nvidia/cudatoolkit7.5/7.5.18-1.0502.10743.2.1
export CUDA="$CUDA_HOME"
export GPU_ARCH=kepler
```

### Download software

```
git clone https://gitlab.com/StanfordLegion/legion.git "$LEGION_DIR"
git clone https://github.com/stanfordhpccenter/soleil-x.git "$SOLEIL_DIR"
```

### Install Legion

```
cd "$LEGION_DIR"/language
USE_CUDA=1 USE_OPENMP=1 USE_GASNET=1 USE_HDF=1 scripts/setup_env.py --llvm-version 35 --terra-url 'https://github.com/StanfordLegion/terra.git' --terra-branch 'luajit2.1'
```

Running
=======

```
cd "$SOLEIL_DIR"/src
[USE_CUDA=0] [USE_HDF=0] make
[USE_CUDA=0] [QUEUE=???] ./soleil.sh ...
```

The `soleil.sh` script forwards all arguments to the `soleil.exec` executable. This includes any options that Soleil itself expects, and any additional options to the Legion runtime.

Currently, Soleil reads the following options:

* `-i <sample>.json`: Provide a configuration file, to be run as an additional sample. See [src/config_schema.lua](src/config_schema.lua) for documentation on the available configuration options.
* `-I <samples>.csv`: Provide a file listing multiple configuration files to run, one per line.
