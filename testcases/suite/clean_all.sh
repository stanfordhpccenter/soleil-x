#!/bin/bash -eu

for DIR in */; do
    echo "$DIR"
    cd "$DIR"
    ./clean.sh
    cd ..
done
