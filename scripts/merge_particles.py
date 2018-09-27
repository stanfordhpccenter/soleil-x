#!/usr/bin/env python

import argparse
import h5py
import os
import re

parser = argparse.ArgumentParser()
parser.add_argument('hdf_file', nargs='+')
args = parser.parse_args()

# Read input metadata
num_files = len(args.hdf_file)
bounds = [] # array(num_files,(int,int))
for i in range(num_files):
    base = os.path.basename(args.hdf_file[i])
    pat = r'([0-9]+)-([0-9]+).hdf'
    m = re.match(pat, base)
    assert(m is not None)
    bounds.append((int(m.group(1)), int(m.group(2))))

# Sanity checks
sorted_bounds = sorted(bounds)
for (prev,next) in zip(sorted_bounds[:-1],sorted_bounds[1:]):
    assert prev[1] == next[0] - 1

# Combine actual data into output HDF file
shape = (sorted_bounds[-1][1] - sorted_bounds[0][0] + 1,)
name = '%s-%s.hdf' % (sorted_bounds[0][0], sorted_bounds[-1][1])
with h5py.File(name, 'w') as fout:
    with h5py.File(args.hdf_file[0], 'r') as fin:
        for fld in fin:
            fout.create_dataset(fld, shape, dtype=fin[fld].dtype)
    for i in range(num_files):
        with h5py.File(args.hdf_file[i], 'r') as fin:
            for fld in fout:
                fout[fld][bounds[i][0]:bounds[i][1]+1] = fin[fld][:]
