#!/usr/bin/env python2

import argparse
import json

EQUALITY_TYPES = [float, int, bool, str, unicode]

def compare(v1, v2, path):
    for typ in EQUALITY_TYPES:
        if type(v1) == typ and type(v2) == typ and v1 == v2:
            return
    if type(v1) == list and type(v2) == list and len(v1) == len(v2):
        for i in range(0, len(v1)):
            compare(v1[i], v2[i], path + '[' + str(i) + ']')
    elif type(v1) == dict and type(v2) == dict:
        for k in v1:
            if k not in v2:
                print 'Only in 1st file:', path + '.' + k
            else:
                compare(v1[k], v2[k], path + '.' + k)
        for k in v2:
            if k not in v1:
                print 'Only in 2nd file:', path + '.' + k
    else:
        print 'Difference in', path, ':', v1, 'vs', v2

parser = argparse.ArgumentParser()
parser.add_argument('f1', type=argparse.FileType('r'))
parser.add_argument('f2', type=argparse.FileType('r'))
args = parser.parse_args()

compare(json.load(args.f1), json.load(args.f2), '<root>')
