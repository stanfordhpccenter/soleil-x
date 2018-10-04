#!/usr/bin/env python2

import argparse
import h5py

parser = argparse.ArgumentParser()
parser.add_argument('hdf_file')
parser.add_argument('density', type=float)
args = parser.parse_args()

f = h5py.File(args.hdf_file)
f['density'] = len(f['temperature']) * [args.density]
f.close()
