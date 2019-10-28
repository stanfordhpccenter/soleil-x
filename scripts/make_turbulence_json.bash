#!/bin/bash
module load cray-python/3.6.5.3
cd $SOLEIL_DIR/testcases/hit_to_openchannel && ./make_cases.sh uncertainties_32.dat hf-pizdaint-24hrs.json --stagger=140

