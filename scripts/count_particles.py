#!/usr/bin/env python2

import argparse
import h5py
import numpy as np

parser = argparse.ArgumentParser()
parser.add_argument('hdf_file')
args = parser.parse_args()

with h5py.File(args.hdf_file, 'r') as f:
    print np.count_nonzero(f['__valid'])
