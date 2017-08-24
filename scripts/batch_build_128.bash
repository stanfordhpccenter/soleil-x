#!/bin/bash -l
#SBATCH --job-name=128build
#SBATCH --mail-user=aheirich@stanford.edu
#SBATCH --time=01:00:00
#SBATCH --nodes=1
#SBATCH --partition=normal

ROOT=/users/aheirich
cd $ROOT
source setup.bash
cd PSAAP
source soleil-m/scripts/do.bash 6


OUTDIR=$SOLEIL_PATH/src/piz_daint_jobs
mkdir -p ${OUTDIR}
cd src
rm -rf piz_daint_jobs/Job_128
$SOLEIL_PATH/scripts/build_one_piz_daint_job.bash 128 8,4,4 taylor_green_vortex_2048_1024_1024.lua 1 7 00:10:00 6 ${OUTDIR}

