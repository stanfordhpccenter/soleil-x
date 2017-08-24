#!/bin/bash -l
#SBATCH --job-name=JOB_ID
#SBATCH --mail-user=aheirich@stanford.edu
#SBATCH --time=TIME_LIMIT
#SBATCH --nodes=NODES
#SBATCH --ntasks-per-core=1
#SBATCH --ntasks-per-node=1
#SBATCH --partition=normal
#SBATCH --constraint=gpu
#SBATCH --account=d51

srun -n 1 -N 1 lscpu

ROOT=/users/aheirich
cd $ROOT/PSAAP
source soleil-m/scripts/do.bash SOLEIL_DIR

export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
export GASNET_BACKTRACE=1
export LD_LIBRARY_PATH="$SOLEIL_PATH/src/piz_daint_jobs/JOB_ID/:$LEGION_PATH/bindings/terra/:$OSMESA_PATH/lib"
export REALM_BACKTRACE=1
export LEGION_BACKTRACE=1
export LEGION_FREEZE_ON_ERROR=1
RUNDIR=/scratch/snx3000/aheirich/runX/JOB_ID

rm -rf ${RUNDIR}
mkdir -p ${RUNDIR}
pushd ${RUNDIR}
mkdir out

echo === copy this script to run dir ===
cp $0 ${RUNDIR}/

COMMAND=" \
srun -n NODES \
        -N NODES \
        --ntasks-per-node 1 \
        --cpu_bind none \
        /lib64/ld-linux-x86-64.so.2 \
	$SOLEIL_PATH/src/piz_daint_jobs/JOB_ID/EXEC \
        -ll:cpu 2 \
        -ll:ocpu 1 \
        -ll:othr 8 \
        -ll:util 2 \
        -ll:dma 2 \
	-ll:io 1 \
        -ll:csize 40000 \
	-ll:rsize 8192 \
        -hl:sched -1 \
        -level legion_prof=2,5 \
	-logfile soleil_prof_% \
        -hl:prof 4 \
	-hl:serializer ascii \
        | tee JOB_ID.log "


echo === run the simulation ===
echo ${COMMAND}
echo ${COMMAND} | /bin/bash

echo === elapsed time ===
sacct -j ${SLURM_JOBID} -o elapsed

echo === compress images ===
pushd ${RUNDIR}/out ; gzip *.ppm ; popd

echo === execution complete, saving slurm log ===
cp $SOLEIL_PATH/src/piz_daint_jobs/JOB_ID/slurm-${SLURM_JOBID}.out ${RUNDIR}/

popd

