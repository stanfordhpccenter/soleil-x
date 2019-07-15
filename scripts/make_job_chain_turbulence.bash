#!/bin/bash
mkdir -p $SCRATCH/job_chain
rm -rf $SCRATCH/job_chain/*
cd $SOLEIL_DIR/testcases/hit_to_openchannel && ./chain_jobs.py $SCRATCH/job_chain viz_turbulence.json 1540 | sed -e "s/NODE=4/NODE=1/" > __job_commands.txt
