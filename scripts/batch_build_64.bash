#!/bin/bash -l
#SBATCH --job-name=64build
#SBATCH --mail-user=aheirich@stanford.edu
#SBATCH --time=01:00:00
#SBATCH --nodes=1
#SBATCH --partition=normal

ROOT=/users/aheirich
cd $ROOT
source setup.bash
cd PSAAP
source soleil-m/scripts/do.bash 5


OUTDIR=$SOLEIL_PATH/src/piz_daint_jobs
mkdir -p ${OUTDIR}
cd src
rm -rf piz_daint_jobs/Job_64
$SOLEIL_PATH/scripts/build_one_piz_daint_job.bash 64 4,4,4 taylor_green_vortex_1024_1024_1024.lua 1 6 00:10:00 5 ${OUTDIR}

