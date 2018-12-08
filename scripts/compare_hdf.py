#!/usr/bin/env python2

import argparse
import h5py
import itertools
import numpy

parser = argparse.ArgumentParser()
parser.add_argument('f1')
parser.add_argument('f2')
args = parser.parse_args()

f1 = h5py.File(args.f1, 'r')
f2 = h5py.File(args.f2, 'r')

for fld in f1:
    assert fld in f2 and f1[fld].shape == f2[fld].shape
for fld in f2:
    assert fld in f1 and f1[fld].shape == f2[fld].shape

for fld in f1:
    max_diff = 0
    def check(x1, x2):
        global max_diff
        if type(x1) == numpy.ndarray or type(x1) == h5py.Dataset:
            assert type(x1) == type(x2) and x1.shape == x2.shape
            for (sub1,sub2) in itertools.izip(x1,x2):
                check(sub1, sub2)
        else:
            max_diff = max(max_diff, abs(x1 - x2))
    check(f1[fld], f2[fld])
    print '%s: max diff = %s' % (fld, max_diff)

f2.close()
f1.close()
