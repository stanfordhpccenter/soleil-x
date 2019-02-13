# save current working directory
CWD=$(pwd)

# Make the dom code if needed
cd $SOLEIL_DIR/src
echo "########################    compile if needed     ############################"
make dom_host.exec
echo ""

# Run the DOM code
cd $CWD
echo "########################     run the dom code     ############################"
#$SOLEIL_DIR/src/my_dom_host.sh $SOLEIL_DIR/testcases/verification/1D_radiation/dom_test.json
$SOLEIL_DIR/src/dom_host.sh $SOLEIL_DIR/testcases/verification/1D_radiation/dom_test.json
echo ""

# Convert the output to numpy format
echo "########################   convert output to numpy   #########################"
python txt_to_numpy.py .
echo ""

# view the output 
echo "########################           plot              #########################"
python viz_dom_output_3d.py ./volume_solution.npz

