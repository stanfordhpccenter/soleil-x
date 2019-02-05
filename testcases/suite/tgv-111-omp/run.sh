#!/bin/bash -eu

USE_CUDA=0 "$SOLEIL_DIR"/src/soleil.sh -i tgv.json &> test.out
