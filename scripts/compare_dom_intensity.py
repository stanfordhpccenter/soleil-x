#!/usr/bin/env python2

import argparse
import itertools
import sys

EPSILON = 10e-10

parser = argparse.ArgumentParser()
parser.add_argument('-v', '--verbose', action='store_true')
parser.add_argument('f1', type=argparse.FileType('r'))
parser.add_argument('f2', type=argparse.FileType('r'))
args = parser.parse_args()

lineno = 0
maxdiff = 0.0
for (l1,l2) in itertools.izip_longest(args.f1, args.f2):
    lineno += 1
    if l2 is None:
        print '%s is longer than %s' % (args.f1.name, args.f2.name)
        sys.exit(1)
    if l1 is None:
        print '%s is shorter than %s' % (args.f1.name, args.f2.name)
        sys.exit(1)
    if len(l1) > 1:
        x1 = float(l1[:-1])
        x2 = float(l2[:-1])
        reldiff = abs((x1-x2)/x1) if x1 != 0.0 else x2
        maxdiff = max(maxdiff, reldiff)
        if reldiff > EPSILON:
            print 'Line %s: %s vs %s' % (lineno, x1, x2)
            sys.exit(1)
if args.verbose:
    print 'Max difference: %s' % maxdiff
