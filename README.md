Setup (generic)
===============

See below for instructions targeting specific systems.

### Prerequisites

* Legion (latest version -- Soleil-X follows the developments in Legion)
* GCC 4.9+ (we need a working `std::regex` library)
* CUDA 7.0+
* Python 2.X

The following are automatically installed during Legion installation:

* LLVM 6.0 (for CUDA 8.0+) or 3.5 (for CUDA 7.X, and better debug info)
* GASNET (custom version)
* Terra (custom version -- we need to use LuaJIT2.1 instead of the default LuaJIT2.0, because the latter exhibits a spurious out-of-memory error when compiling large Regent programs)
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
USE_CUDA=? USE_OPENMP=? USE_GASNET=? USE_HDF=? scripts/setup_env.py --llvm-version 38 --terra-url 'https://github.com/StanfordLegion/terra.git' --terra-branch 'luajit2.1'
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
* `-o <out_dir>`: Specify an output directory for the executable (default is current directory).

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
USE_CUDA=0 USE_OPENMP=1 USE_GASNET=0 USE_HDF=1 scripts/setup_env.py --llvm-version 35 --terra-url 'https://github.com/StanfordLegion/terra.git' --terra-branch 'luajit2.1'
```

### Compile Soleil-X

```
cd "$SOLEIL_DIR"/src
make
```

Setup (Sapling @ Stanford)
==========================

### Add to shell startup

```
# Module loads
module load mpi/openmpi/1.8.2
module load cuda/7.0
# Build config
export CONDUIT=ibv
export CC=gcc-4.9
export CXX=g++-4.9
# Path setup
export LEGION_DIR=???
export HDF_ROOT="$LEGION_DIR"/language/hdf/install
export SOLEIL_DIR=???
export SCRATCH=/scratch/oldhome/`whoami`
# CUDA config
export CUDA_HOME=/usr/local/cuda-7.0
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
USE_CUDA=1 USE_OPENMP=1 USE_GASNET=1 USE_HDF=1 scripts/setup_env.py --llvm-version 35 --terra-url 'https://github.com/StanfordLegion/terra.git' --terra-branch 'luajit2.1'
```

### Compile Soleil-X

```
cd "$SOLEIL_DIR"/src
make
```

Setup (Sherlock @ Stanford)
===========================

### Add to shell startup

```
# Module loads
module load gcc/6.3.0
module load cuda/9.2.148
module load openmpi/2.0.2
# Build config
export CONDUIT=ibv
export CC=gcc
export CXX=g++
export LD_FLAGS=-lpmi2
# Path setup
export LEGION_DIR=???
export HDF_ROOT="$LEGION_DIR"/language/hdf/install
export SOLEIL_DIR=???
# CUDA config
export CUDA_HOME=/share/software/user/open/cuda/9.2.148
export CUDA="$CUDA_HOME"
export GPU_ARCH=pascal
```

### Download software

```
git clone https://gitlab.com/StanfordLegion/legion.git "$LEGION_DIR"
git clone https://github.com/stanfordhpccenter/soleil-x.git "$SOLEIL_DIR"
```

### Install Legion

We build Legion in a SLURM job, because processes on the login node are restricted to 1 core, and we also require a node with a proper CUDA installation.

```
cd "$LEGION_DIR"/language
USE_CUDA=1 USE_OPENMP=1 USE_GASNET=1 USE_HDF=1 srun -N 1 -c 10 -p aaiken --gres=gpu:4 scripts/setup_env.py --llvm-version 60 --terra-url 'https://github.com/StanfordLegion/terra.git' --terra-branch 'luajit2.1'
```

### Compile Soleil-X

Soleil-X must similarly be built in a SLURM job:

```
cd "$SOLEIL_DIR"/src
srun -N 1 -c 10 -p aaiken --gres=gpu:4 make
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

```
# Disable PMI in GASnet, because the PMI library is missing on Certainty.
git clone https://github.com/StanfordLegion/gasnet.git $LEGION_DIR/language/gasnet
cd "$LEGION_DIR"/language/gasnet
sed -i 's|$(GASNET_VERSION)/configure --prefix=|$(GASNET_VERSION)/configure --disable-pmi --prefix=|g' Makefile
make
# Rest of compilation as normal
cd "$LEGION_DIR"/language
USE_CUDA=1 USE_OPENMP=1 USE_GASNET=1 USE_HDF=1 scripts/setup_env.py --llvm-version 38 --terra-url 'https://github.com/StanfordLegion/terra.git' --terra-branch 'luajit2.1'
```

### Compile Soleil-X

```
cd "$SOLEIL_DIR"/src
make
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
USE_CUDA=1 USE_OPENMP=1 USE_GASNET=1 USE_HDF=1 scripts/setup_env.py --llvm-version 60 --terra-url 'https://github.com/StanfordLegion/terra.git' --terra-branch 'luajit2.1'
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
export SCRATCH=/gpfs/alpinetds/proj-shared/csc275/stanford
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

```
cd "$LEGION_DIR"/language
USE_CUDA=1 USE_OPENMP=1 USE_GASNET=1 USE_HDF=1 scripts/setup_env.py --llvm-version 38 --terra-url 'https://github.com/StanfordLegion/terra.git' --terra-branch 'luajit2.1-ppc64'
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
module load cudatoolkit/9.2.148_3.19-6.0.7.1_2.1__g3d9acc8
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
export CUDA_HOME=/opt/nvidia/cudatoolkit9.2/9.2.148_3.19-6.0.7.1_2.1__g3d9acc8
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
USE_CUDA=1 USE_OPENMP=1 USE_GASNET=1 USE_HDF=1 scripts/setup_env.py --llvm-version 60 --terra-url 'https://github.com/StanfordLegion/terra.git' --terra-branch 'luajit2.1'
```

### Compile Soleil-X

```
cd "$SOLEIL_DIR"/src
make
```
