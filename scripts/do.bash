#!/bin/bash
if [[ "$1" == "" ]]
then
  echo give me a number 2..9
  exit
fi

source soleil-m/scripts/setup.bash
export SOLEIL_PATH=${SOLEIL_PATH}-$1
cd $SOLEIL_PATH
git stash
git pull

