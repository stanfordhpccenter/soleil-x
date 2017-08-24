#!/bin/bash -l
#SBATCH --job-name=256build
#SBATCH --mail-user=aheirich@stanford.edu
#SBATCH --time=01:00:00
#SBATCH --nodes=1
#SBATCH --partition=normal

ROOT=/users/aheirich
cd $ROOT
source setup.bash
cd PSAAP
source soleil-m/scripts/do.bash 7


OUTDIR=$SOLEIL_PATH/src/piz_daint_jobs
mkdir -p ${OUTDIR}
cd src
rm -rf piz_daint_jobs/Job_256
$SOLEIL_PATH/scripts/build_one_piz_daint_job.bash 256 8,8,4 taylor_green_vortex_2048_2048_1024.lua 1 8 00:10:00 7 ${OUTDIR}

