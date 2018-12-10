#!/bin/bash
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$LEGION_DIR/bindings/regent:$LEGION_DIR/language/hdf/install/lib
mpiexec -N 1 bug.exec -i ../testcases/tgv.json
