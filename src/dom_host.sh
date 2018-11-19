#!/bin/bash -eu

export QUEUE="${QUEUE:-}"
export AFTER="${AFTER:-}"
export USE_CUDA="${USE_CUDA:-1}"
export PROFILED_RANKS="${PROFILED_RANKS:-0}"
export PROFILE="${PROFILE:-0}"
export GASNET_BACKTRACE="${GASNET_BACKTRACE:-1}"

export EXECUTABLE="$SOLEIL_DIR"/src/dom_host.exec
export MINUTES=10
export NUM_RANKS=1

source "$SOLEIL_DIR"/src/run.sh
