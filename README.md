Setup (generic)
===============

See below for instructions targeting specific systems.

### Prerequisites

* Legion (latest version -- Soleil-X follows the developments in Legion)
* GCC 4.9+ (we need a working `std::regex` library)
* CUDA 7.5+
* Python 2.X

The following are automatically installed during Legion installation:

* LLVM 3.8 (for CUDA 8.0+) or 3.5 (for CUDA 7.5, and better debug info)
* GASNET (custom Legion-specific version)
* Terra (custom Legion-specific version)
* HDF5 (any recent version)

### Add to shell startup

Normally you'd need to edit file `~/.bashrc`. Replace the `???` depending on your system.

```
# Module loads (if necessary)
...
# Build config (if necessary, for Legion or Soleil-X)
...
# Path setup (mandatory)
export LEGION_DIR=???
export HDF_ROOT="$LEGION_DIR"/language/hdf/install
export SOLEIL_DIR=???
[export SCRATCH=???]
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
USE_CUDA=? USE_OPENMP=? USE_GASNET=? USE_HDF=? scripts/setup_env.py --llvm-version 38
```

See [Elliott's instructions](https://docs.google.com/document/d/1Qkl6r-1ZIb8WyH1f_WZbKgjp3due_Q8UiWKLh_nG1ec/edit) for more help.

### Compile Soleil-X

```
cd "$SOLEIL_DIR"/src
[USE_CUDA=0] [USE_HDF=0] make
```

Running
=======

```
cd "$SOLEIL_DIR"/src
./soleil.sh ...
```

The [src/soleil.sh](src/soleil.sh) script accepts some options through the environment (see the top of that file for details), and forwards all command-line arguments to the Soleil-X executable and the Legion runtime (each will ignore options it doesn't recognize).

Currently, Soleil-X reads the following options:

* `-i <config>.json`: Provide a case configuration file, to be run as an additional sample. See [src/config_schema.lua](src/config_schema.lua) for documentation on the available options (`Config` struct).
* `-m <multi-config>.json`: Provide a two-case configuration file, to be run as two connected samples. See [src/config_schema.lua](src/config_schema.lua) for documentation on the available options (`MultiConfig` struct).
* `-o <out_dir>`: Specify an output directory for the executable (if not defined, we use a new directory under `$SCRATCH` if that is defined, otherwise we use the current directory).

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
USE_CUDA=0 USE_OPENMP=1 USE_GASNET=0 USE_HDF=1 scripts/setup_env.py --llvm-version 38
```

### Compile Soleil-X

```
cd "$SOLEIL_DIR"/src
make
```

Setup (Sapling @ Stanford)
==========================

See [Elliott's guide](https://github.com/StanfordLegion/sapling-guide) for general guidance on using sapling.

### Add to shell startup

```
# Module loads
module load slurm/20.11.4
# Build config
export CONDUIT=ibv
export CC=gcc
export CXX=g++
export CC_FLAGS=-std=c++14
export NVCC_FLAGS=-std=c++14
# Path setup
export LEGION_DIR=???
export HDF_ROOT="$LEGION_DIR"/language/hdf/install
export SOLEIL_DIR=???
# CUDA config
export CUDA_HOME=/usr/local/cuda-11.1
export CUDA="$CUDA_HOME"
export GPU_ARCH=pascal
```

### Download software

```
git clone https://gitlab.com/StanfordLegion/legion.git "$LEGION_DIR"
git clone https://github.com/stanfordhpccenter/soleil-x.git "$SOLEIL_DIR"
```

### Install Legion

We can't build Legion on the head node, because it has a different architecture than the compute nodes. Currently we have to manually reload modules when moving to a compute node, because the version installed on the head node is incompatible with that on the compute nodes.

```
cd "$LEGION_DIR"/language
USE_CUDA=1 USE_OPENMP=1 USE_GASNET=1 USE_HDF=1 srun --exclusive -N 1 -p gpu bash -c '. "$SOLEIL_DIR"/src/sapling_reload_modules.sh; scripts/setup_env.py --llvm-version 60'
```

### Compile Soleil-X

Similarly, we build Soleil on a compute node.

```
cd "$SOLEIL_DIR"/src
srun --exclusive -N 1 -p gpu bash -c '. "$SOLEIL_DIR"/src/sapling_reload_modules.sh; make'
```

Setup (Certainty @ Stanford)
============================

Certainty no longer has a GPU queue.

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

```
# Disable PMI in GASnet, because the PMI library is missing on Certainty.
git clone https://github.com/StanfordLegion/gasnet.git $LEGION_DIR/language/gasnet
cd "$LEGION_DIR"/language/gasnet
sed -i 's|$(GASNET_VERSION)/configure --prefix=|$(GASNET_VERSION)/configure --disable-pmi --prefix=|g' Makefile
make
# Rest of compilation as normal
cd "$LEGION_DIR"/language
USE_CUDA=0 USE_OPENMP=1 USE_GASNET=1 USE_HDF=1 scripts/setup_env.py --llvm-version 38
```

### Compile Soleil-X

```
cd "$SOLEIL_DIR"/src
USE_CUDA=0 make
```

Setup (Titan @ ORNL)
====================

### Add to shell startup

Install Legion and Soleil-X on the `/ccs/proj` filesystem, not your home directory (home directories are not mounted on the compute nodes).

We set `MARCH`, the processor architecture that Legion will be built for, to one that is compatible with both the login node and the compute nodes. This is done so that `libregent.so` can be used by both the final executable, when it's executing on a compute node, and the Regent compiler, which needs to link against it while compiling Soleil-X on the login node.

```
# Module loads
module load python/2.7.9
module load cudatoolkit/9.1.85_3.10-1.0502.df1cc54.3.1
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
export SCRATCH="$PROJWORK"/csc188/stanford
# CUDA config
export CUDA_HOME=/opt/nvidia/cudatoolkit9.1/9.1.85_3.10-1.0502.df1cc54.3.1
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
USE_CUDA=1 USE_OPENMP=1 USE_GASNET=1 USE_HDF=1 scripts/setup_env.py --llvm-version 38
```

### Compile Soleil-X

```
cd "$SOLEIL_DIR"/src
make
```

Setup (Summit @ ORNL)
=====================

### Add to shell startup

Install Legion and Soleil-X on the `/ccs/proj` filesystem, not your home directory (home directories are not mounted on the compute nodes).

```
# Module loads
module load gcc/6.4.0
module load cuda/9.0.184
# Build config
export CC=gcc
export CXX=g++
export CONDUIT=ibv
# Path setup
export LEGION_DIR=???
export HDF_ROOT="$LEGION_DIR"/language/hdf/install
export SOLEIL_DIR=???
export SCRATCH=/gpfs/alpine/proj-shared/csc335/stanford
# CUDA config
export CUDA_HOME=/sw/summit/cuda/9.0.184
export CUDA="$CUDA_HOME"
export GPU_ARCH=volta
```

### Download software

```
git clone https://gitlab.com/StanfordLegion/legion.git "$LEGION_DIR"
git clone https://github.com/stanfordhpccenter/soleil-x.git "$SOLEIL_DIR"
```

### Install Legion

The system restricts the amount of memory we can use on a login node, so we need to build on a compute node.

```
cd "$LEGION_DIR"/language
CC_FLAGS='-DMAX_NUM_NODES=4096' TERRA_USE_PUC_LUA=1 USE_CUDA=1 USE_OPENMP=1 USE_GASNET=1 USE_HDF=1 bsub -W 60 -nnodes 1 -P CSC335 -I scripts/setup_env.py --llvm-version 38 --terra-branch 'puc_lua_master'
```

### Compile Soleil-X

```
cd "$SOLEIL_DIR"/src
make
```

Setup (PizDaint @ ETH)
======================

### Add to shell startup

```
# Module loads
module swap PrgEnv-cray PrgEnv-gnu
module load daint-gpu
module load cudatoolkit
# Build config
export CC=cc
export CXX=CC
export HOST_CC=gcc
export HOST_CXX=g++
# Path setup
export LEGION_DIR=???
export HDF_ROOT="$LEGION_DIR"/language/hdf/install
export SOLEIL_DIR=???
# CUDA config
export CUDA_HOME="$CUDATOOLKIT_HOME"
export CUDA="$CUDA_HOME"
export GPU_ARCH=pascal
```

### Download software

```
git clone https://gitlab.com/StanfordLegion/legion.git "$LEGION_DIR"
git clone https://github.com/stanfordhpccenter/soleil-x.git "$SOLEIL_DIR"
```

### Install Legion

```
cd "$LEGION_DIR"/language
USE_CUDA=1 USE_OPENMP=1 USE_GASNET=1 USE_HDF=1 scripts/setup_env.py --llvm-version 60
```

### Compile Soleil-X

```
cd "$SOLEIL_DIR"/src
make
```

Setup (Lassen/Sierra @ LLNL)
============================

### Add to shell startup

Install Legion and Soleil-X on the `/usr/workspace` filesystem; your home directory has a low quota.

```
# Module loads
module load gcc/7.3.1
module load cuda/9.2.148
# Build config
export CC=gcc
export CXX=g++
export CONDUIT=ibv
# Path setup
export LEGION_DIR=???
export HDF_ROOT="$LEGION_DIR"/language/hdf/install
export SOLEIL_DIR=???
export SCRATCH=/p/gpfs1/`whoami`
# CUDA config
export CUDA_HOME=/usr/tce/packages/cuda/cuda-9.2.148
export CUDA="$CUDA_HOME"
export GPU_ARCH=volta
```

### Download software

```
git clone https://gitlab.com/StanfordLegion/legion.git "$LEGION_DIR"
git clone https://github.com/stanfordhpccenter/soleil-x.git "$SOLEIL_DIR"
```

### Install Legion

We need to go through the `lalloc` utility script to build on a compute node.

```
cd "$LEGION_DIR"/language
CC_FLAGS='-DMAX_NUM_NODES=4096' TERRA_USE_PUC_LUA=1 USE_CUDA=1 USE_OPENMP=1 USE_GASNET=1 USE_HDF=1 lalloc 1 -W 720 scripts/setup_env.py --llvm-version 38 --terra-branch 'puc_lua_master'
```

### Compile Soleil-X

Soleil-X must similarly be built on a compute node.

```
cd "$SOLEIL_DIR"/src
lalloc 1 -W 720 make
```
