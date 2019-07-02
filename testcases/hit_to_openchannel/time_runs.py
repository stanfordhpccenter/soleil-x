#!/usr/bin/env python2

import argparse
import glob
import os
import subprocess

parser = argparse.ArgumentParser()
parser.add_argument('-t', '--target-iters', type=int, default=5000)
parser.add_argument('jobids', nargs='+')
args = parser.parse_args()

avg_tmin = 0.0
avg_tmax = 0.0
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
    avg_tmin += tmin / len(args.jobids)
    avg_tmax += tmax / len(args.jobids)
    print '%s: %.3f %.3f' % (j, tmin, tmax)
print 'Average: %.3f %.3f' % (avg_tmin, avg_tmax)
