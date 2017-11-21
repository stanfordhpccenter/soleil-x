#!/bin/bash

if [[ "`pwd`" != "$SOLEIL_PATH/src" ]]
then
  echo === Run this from $SOLEIL_PATH/src ===
  exit -1
fi 
NUM_FRAGMENTS=$1
if [[ "$NUM_FRAGMENTS" == "" ]]
then
  NUM_FRAGMENTS=2
fi
NUM_TREE_LEVELS=$2
if [[ "$NUM_TREE_LEVELS" == "" ]]
then
  NUM_TREE_LEVELS=2
fi
TESTCASE=$3
if [[ "$TESTCASE" == "" ]]
then
  CASEDIR=taylor_with_smaller_particles
  # CASENAME=taylor_green_vortex_512_512_256.lua
  CASENAME=taylor_green_vortex_256_256_256.lua
  TESTCASE=${SOLEIL_PATH}/testcases/${CASEDIR}/${CASENAME}
fi
MESH=$4
if [[ "$MESH" == "" ]]
then
  MESH=2,2,1
fi
MAX_ITER=$5
if [[ "$MAX_ITER" == "" ]]
then
  MAX_ITER=9
fi


echo === Build ${NUM_FRAGMENTS} fragments, ${NUM_TREE_LEVELS} tree levels, test case ${TESTCASE} MESH ${MESH} MAX_ITER ${MAX_ITER} ===

START=`date`
echo start ${START}
cd $SOLEIL_PATH/src 
EXEC=soleil_${NUM_FRAGMENTS}_${NUM_TREE_LEVELS}.exec
rm -f testcase.lua
cat ${TESTCASE} | sed -e "s:max_iter =.*:max_iter = ${MAX_ITER},:" > testcase.lua

cat tmp_src/viz.rg | sed -e "s/NUM_FRAGMENTS/${NUM_FRAGMENTS}/" > viz.rg

COMPILE_COMMAND=\
"CC=gcc CXX=g++ USE_HDF=0 DEBUG=0 SAVEOBJ=1 OBJNAME=${EXEC} \
	$LISZT_PATH/liszt-legion.sh $SOLEIL_PATH/src/soleil-master.t \
	-i testcase.lua \
	-fparallelize 1 \
	-fparallelize-dop ${MESH} \
	-fopenmp 1 \
        -fopenmp-strict 1 \
	-fflow-spmd 1 \
        -fflow-spmd-shardsize 1 "

echo $COMPILE_COMMAND
echo $COMPILE_COMMAND | /bin/bash


echo start ${START}
echo end `date`

