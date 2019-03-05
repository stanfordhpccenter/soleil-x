#!/bin/bash

# save current working directory
CWD=$(pwd)

# Make the dom code if needed
cd $SOLEIL_DIR/src
echo "########################    compile if needed     ############################"
make dom_host.exec
echo ""

cd $CWD

# remove any old data
echo "########################     remove old data      ############################"
rm volume_solution.txt
rm volume_solution.npz
echo ""

# Run the DOM code
echo "########################     run the dom code     ############################"
$SOLEIL_DIR/src/my_dom_host.sh $SOLEIL_DIR/testcases/verification/radiation/1D_radiation/dom_test.json
#$SOLEIL_DIR/src/dom_host.sh $SOLEIL_DIR/testcases/verification/radiation/1D_radiation/dom_test.json
echo ""

# Convert the output to numpy format
echo "########################   convert output to numpy   #########################"
python txt_to_numpy.py .
echo ""

# view the output 
echo "########################           plot              #########################"
#python viz_dom_output_3d.py ./volume_solution.npz

python compare_1D_postprocess.py volume_solution.npz
