#!/bin/bash

NODES=$1
MESH=$2
TESTCASE=$3
NUM_FRAGMENTS=$4
NUM_TREE_LEVELS=$5
TIME_LIMIT=$6
SOLEIL_DIR=$7
OUTDIR=$8

if [[ "$TESTCASE" == "" ]]
then
  TESTCASE=taylor_green_vortex_256_256_256.lua
fi

cd $SOLEIL_PATH/src

$SOLEIL_PATH/scripts/build_multinode_openmp.bash $NUM_FRAGMENTS $NUM_TREE_LEVELS $SOLEIL_PATH/testcases/taylor_with_smaller_particles/$TESTCASE $MESH

JOB_ID=Job_$NODES
mkdir -p ${OUTDIR}/${JOB_ID}
EXEC=soleil_${NUM_FRAGMENTS}_${NUM_TREE_LEVELS}.exec
mv ./${EXEC} ${OUTDIR}/${JOB_ID}
mv ./*.so ${OUTDIR}/${JOB_ID}
SCRIPT="${OUTDIR}/${JOB_ID}/${JOB_ID}_piz_daint.bash"
COMMAND="cat $SOLEIL_PATH/scripts/piz_daint_multinode_job.bash | \
  sed -e \"s/TIME_LIMIT/${TIME_LIMIT}/g\" | \
  sed -e \"s/TESTCASE/${TESTCASE}/g\" | \
  sed -e \"s/NODES/${NODES}/g\" | \
  sed -e \"s/JOB_ID/${JOB_ID}/g\" | \
  sed -e \"s/EXEC/${EXEC}/g\" | \
  sed -e \"s:SCRIPT:${SCRIPT}:g\" | \
  sed -e \"s/SOLEIL_DIR/${SOLEIL_DIR}/g\" \
  > ${SCRIPT} "
echo ${COMMAND}
echo ${COMMAND} | /bin/bash


