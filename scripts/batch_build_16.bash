#!/bin/bash -l
#SBATCH --job-name=16build
#SBATCH --mail-user=aheirich@stanford.edu
#SBATCH --time=01:00:00
#SBATCH --nodes=1
#SBATCH --partition=normal

ROOT=/users/aheirich
cd $ROOT
source setup.bash
cd PSAAP
source soleil-m/scripts/do.bash 3


OUTDIR=$SOLEIL_PATH/src/piz_daint_jobs
mkdir -p ${OUTDIR}
cd src
rm -rf piz_daint_jobs/Job_16
$SOLEIL_PATH/scripts/build_one_piz_daint_job.bash 16 4,2,2 taylor_green_vortex_1024_512_512.lua 1 4 00:10:00 3 ${OUTDIR}

