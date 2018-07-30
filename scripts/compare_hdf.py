#!/usr/bin/env python

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

def check():
    for fld in f1:
        if fld not in f2:
            return False
        if f1[fld].shape != f2[fld].shape:
            return False
        for (x1,x2) in itertools.izip(f1[fld],f2[fld]):
            if type(x1) == numpy.ndarray:
                if x1.shape != x2.shape:
                     return False
                for (e1,e2) in itertools.izip(x1,x2):
                     if e1 != e2:
                         return False
            else:
                if x1 != x2:
                    return False
    for fld in f2:
        if fld not in f1:
            return False
    return True

if check():
    print 'same'
else:
    print 'different'

f2.close()
f1.close()
