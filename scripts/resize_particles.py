#!/usr/bin/env python2

import argparse
import h5py
import numpy as np

parser = argparse.ArgumentParser()
parser.add_argument('hdf_in')
parser.add_argument('hdf_out')
parser.add_argument('new_size', type=int)
args = parser.parse_args()

with h5py.File(args.hdf_in, 'r') as fin:
    for fld in fin:
        old_size = fin[fld].shape[0]
        break
    for fld in fin:
        assert old_size == fin[fld].shape[0]
    if old_size > args.new_size:
        assert np.count_nonzero(fin['__valid'][args.new_size:]) == 0
    with h5py.File(args.hdf_out, 'w') as fout:
        for fld in fin:
            fout.create_dataset(fld, shape=(args.new_size,), dtype=fin[fld].dtype)
            if old_size > args.new_size:
                fout[fld][:] = fin[fld][:args.new_size]
            else:
                fout[fld][:old_size] = fin[fld][:]
                if fld == '__valid': # otherwise the fill value is unspecified
                    fout[fld][old_size:] = 0
