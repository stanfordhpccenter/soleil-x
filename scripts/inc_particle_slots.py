#!/usr/bin/env python2

import argparse
import h5py

parser = argparse.ArgumentParser()
parser.add_argument('hdf_in')
parser.add_argument('hdf_out')
parser.add_argument('new_size', type=int)
args = parser.parse_args()

with h5py.File(args.hdf_in, 'r') as fin:
    with h5py.File(args.hdf_out, 'w') as fout:
        for fld in fin:
            old_size = fin[fld].shape[0]
            assert old_size <= args.new_size
            fout.create_dataset(fld, shape=(args.new_size,), dtype=fin[fld].dtype)
            if fld == '__valid':
                fout[fld][:] = 0
            fout[fld][:old_size] = fin[fld][:]
