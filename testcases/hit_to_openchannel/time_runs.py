#!/usr/bin/env python2

import argparse
import glob
import numpy as np
import os
import subprocess
import sys

parser = argparse.ArgumentParser()
parser.add_argument('-t', '--target-iters', type=int, default=5000)
parser.add_argument('-v', '--verbose', action='store_true')
parser.add_argument('jobids', nargs='+')
args = parser.parse_args()

tmin_all = []
tmax_all = []
for j in args.jobids:
    if args.verbose:
        sys.stdout.write('%s: ' % j)
        sys.stdout.flush()
    status = subprocess.check_output([os.path.join(sys.path[0], 'job_status.sh'), j + '.out'])
    if status != 'done\n':
        if args.verbose:
            sys.stdout.write(status)
            sys.stdout.flush()
        continue
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
    if args.verbose:
        print 'Min: %.3f, Max: %.3f' % (tmin, tmax)
print 'Runs: %d, MinAvg: %.3f ~ %.2f%%, MaxAvg: %.3f ~ %.2f%%' % (len(tmax_all), np.mean(tmin_all), 100.0*np.std(tmin_all,ddof=1)/np.mean(tmin_all), np.mean(tmax_all), 100.0*np.std(tmax_all,ddof=1)/np.mean(tmax_all))
