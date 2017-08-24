#!/bin/bash -l
#SBATCH --job-name=512build
#SBATCH --mail-user=aheirich@stanford.edu
#SBATCH --time=01:00:00
#SBATCH --nodes=1
#SBATCH --partition=normal

ROOT=/users/aheirich
cd $ROOT
source setup.bash
cd PSAAP
source soleil-m/scripts/do.bash 8


OUTDIR=$SOLEIL_PATH/src/piz_daint_jobs
mkdir -p ${OUTDIR}
cd src
rm -rf piz_daint_jobs/Job_512
$SOLEIL_PATH/scripts/build_one_piz_daint_job.bash 512 8,8,8 taylor_green_vortex_2048_2048_2048.lua 1 9 00:10:00 8 ${OUTDIR}

