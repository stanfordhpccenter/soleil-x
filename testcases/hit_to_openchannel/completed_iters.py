#!/usr/bin/env python2

import argparse
import glob
import itertools
import os
import re
import subprocess
import sys

parser = argparse.ArgumentParser()
parser.add_argument('target_time', type=float)
parser.add_argument('jobids', nargs='+')
args = parser.parse_args()

for j in args.jobids:
    sys.stdout.write('%s: ' % j)
    sys.stdout.flush()
    status = subprocess.check_output([os.path.join(sys.path[0], 'job_status.sh'), j + '.out'])
    if status != 'done\n':
        sys.stdout.write(status)
        sys.stdout.flush()
        continue
    iters = {}
    for console in glob.glob('%s/%s/sample*/console.txt' % (os.environ['SCRATCH'], j)):
        m = re.search(r'sample([0-9]+)', console)
        sample = int(m.group(1))
        with open(console) as fin:
            for line in itertools.islice(fin, 1, None):
                toks = line.split()
                t = float(toks[2])
                if t > args.target_time:
                    iters[sample] = int(toks[0])
                    break
    max_sample = max(iters.keys())
    print ' '.join([str(iters[i]) for i in range(max_sample+1)])
