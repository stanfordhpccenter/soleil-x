#!/bin/bash -u

for DIR in */; do
    echo "$DIR"
    cd "$DIR"
    ./check.sh
    cd ..
done
