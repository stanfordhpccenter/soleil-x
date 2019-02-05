#!/bin/bash -eu

for DIR in */; do
    echo "$DIR"
    cd "$DIR"
    ./run.sh
    cd ..
done
