#!/bin/bash -l
#SBATCH --job-name=make_soleil
#SBATCH --mail-user=aheirich@slac.stanford.edu
#SBATCH --time=8:00:00
#SBATCH --nodes=1
#SBATCH --partition=normal
#SBATCH --constraint=gpu

source ~/setup.bash
cd $SOLEIL_DIR
cd src
make
