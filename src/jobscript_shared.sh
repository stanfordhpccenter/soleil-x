#!/bin/bash -eu

if [[ ! -z "${PBS_JOBID:-}" ]]; then
    JOBID="$PBS_JOBID"
elif [[ ! -z "${SLURM_JOBID:-}" ]]; then
    JOBID="$SLURM_JOBID"
elif [[ ! -z "${LSB_JOBID:-}" ]]; then
    JOBID="$LSB_JOBID"
else
    JOBID="$(date +%s)"
fi

# Create a directory on the scratch filesystem, for output (if SCRATCH is
# defined, and the user hasn't specified an output directory explicitly).
OUT_DIR=
OUT_DIR_FOLLOWS=false
for ARG in $ARGS; do
    if [[ "$OUT_DIR_FOLLOWS" == true ]]; then
        OUT_DIR="$ARG"
        break
    elif [[ "$ARG" == "-o" ]]; then
        OUT_DIR_FOLLOWS=true
    fi
done
if [[ ! -z "${SCRATCH:-}" && -z "$OUT_DIR" ]]; then
    OUT_DIR="$SCRATCH"/"$JOBID"
    mkdir "$OUT_DIR"
    ARGS="$ARGS -o $OUT_DIR"
    echo "Redirecting output to $OUT_DIR"
fi
