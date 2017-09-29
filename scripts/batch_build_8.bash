#!/bin/bash -l
#SBATCH --job-name=8build
#SBATCH --mail-user=aheirich@stanford.edu
#SBATCH --time=01:00:00
#SBATCH --nodes=1
#SBATCH --partition=normal

ROOT=/users/aheirich
cd $ROOT
source setup.bash
cd PSAAP
source soleil-m/scripts/do.bash 2


OUTDIR=$SOLEIL_PATH/src/piz_daint_jobs
mkdir -p ${OUTDIR}
cd src
rm -rf piz_daint_jobs/Job_8
$SOLEIL_PATH/scripts/build_one_piz_daint_job.bash 8 2,2,2 taylor_green_vortex_512_512_512.lua 1 3 00:10:00 2 ${OUTDIR}

