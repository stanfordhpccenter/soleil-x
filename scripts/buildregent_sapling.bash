#!/bin/bash
pushd $LEGION_DIR/language
CC_FLAGS=-std=c++11 DEBUG=1 scripts/setup_env.py --llvm-version 38 --terra-url 'https://github.com/StanfordLegion/terra.git' --terra-branch 'luajit2.1'
popd
