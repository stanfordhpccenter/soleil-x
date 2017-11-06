#!/bin/bash -l
#SBATCH --job-name=32build
#SBATCH --mail-user=aheirich@stanford.edu
#SBATCH --time=01:00:00
#SBATCH --nodes=1
#SBATCH --partition=normal

ROOT=/users/aheirich
cd $ROOT
source setup.bash
cd PSAAP
source soleil-master/scripts/do.bash 4


OUTDIR=$SOLEIL_PATH/src/piz_daint_jobs
mkdir -p ${OUTDIR}
cd src
rm -rf piz_daint_jobs/Job_32
$SOLEIL_PATH/scripts/build_one_piz_daint_job.bash 32 4,4,2 taylor_green_vortex_1024_1024_512.lua 1 5 00:10:00 4 ${OUTDIR}

