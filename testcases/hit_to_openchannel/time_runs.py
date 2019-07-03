#!/usr/bin/env python2

import argparse
import glob
import numpy as np
import os
import subprocess

parser = argparse.ArgumentParser()
parser.add_argument('-t', '--target-iters', type=int, default=5000)
parser.add_argument('jobids', nargs='+')
args = parser.parse_args()

tmin_all = []
tmax_all = []
for j in args.jobids:
    tmin = float('inf')
    tmax = 0.0
    for console in glob.glob('%s/%s/sample*/console.txt' % (os.environ['SCRATCH'], j)):
        last_line = subprocess.check_output(['tail', '-n', '1', console])
        toks = last_line.split()
        iters = int(toks[0])
        assert iters == args.target_iters
        t = float(toks[2])
        tmin = min(t, tmin)
        tmax = max(t, tmax)
    tmin_all.append(tmin)
    tmax_all.append(tmax)
    print '%s: %.3f %.3f' % (j, tmin, tmax)
print 'Average: %.3f, StdDev: %.3f' % (np.mean(tmax_all), np.std(tmax_all, ddof=1))
