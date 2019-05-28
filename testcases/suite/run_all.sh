#!/bin/bash -u

for DIR in */; do
    echo "$DIR"
    cd "$DIR"
    if ! ./run.sh; then
	echo "FAIL"
    fi
    cd ..
done
