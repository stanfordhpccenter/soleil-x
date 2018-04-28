#!/usr/bin/env python

import argparse
import itertools

parser = argparse.ArgumentParser()
parser.add_argument('f1')
parser.add_argument('f2')
args = parser.parse_args()

total = 0
maxDiff = 0.0
avgDiff = 0.0
with open(args.f1) as f1:
  with open(args.f2) as f2:
    for l1,l2 in itertools.izip(f1, f2):
      if len(l1) > 1:
        x1 = float(l1[:-1])
        x2 = float(l2[:-1])
        diff = abs((x1 - x2) / x1)
        maxDiff = max(maxDiff, diff)
        avgDiff += diff
        total += 1
avgDiff /= total
print "Max diff: " + str(maxDiff)
print "Avg diff: " + str(avgDiff)
