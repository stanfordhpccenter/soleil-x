#!/bin/bash
mkdir -p $SCRATCH/job_chain
rm -rf $SCRATCH/job_chain/*
cd $SOLEIL_DIR/testcases/hit_to_openchannel && ./chain_jobs.py $SCRATCH/job_chain use_this_for_turbulence_visualization.json 1500 > __job_commands.txt
