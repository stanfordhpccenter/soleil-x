#!/bin/bash -eu

export QUEUE=
export USE_CUDA=1
export REALM_BACKTRACE=1
export EXECUTABLE="$SOLEIL_DIR"/src/dom_host.exec
export MINUTES=10
export NUM_RANKS=1

source "$SOLEIL_DIR"/src/run.sh
