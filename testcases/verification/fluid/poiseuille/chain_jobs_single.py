#!/usr/bin/env python2

import argparse
import json
import math
import os

parser = argparse.ArgumentParser()
parser.add_argument('out_dir')
parser.add_argument('json_file', nargs='+', type=argparse.FileType('r'))
parser.add_argument('iters_per_job', type=int)
args = parser.parse_args()

cases = []
max_iter = []
for fin in args.json_file:
    case = json.load(fin)
    cases.append(case)
    max_iter.append(int(case['Integrator']['maxIter']))

total_max_iter = max(max_iter)
num_jobs = int(math.ceil(float(total_max_iter) / args.iters_per_job))
for i in range(0, num_jobs):
    os.mkdir('job%d' % i)
    os.mkdir('%s/job%d' % (args.out_dir, i))
    for j in range(0, len(cases)):
        iter_i = min(i * args.iters_per_job, max_iter[j])
        iter_i_1 = min((i+1) * args.iters_per_job, max_iter[j])
        config = cases[j]
        dt = float(config['Integrator']['fixedDeltaTime'])
        config['Integrator']['startIter'] = iter_i
        config['Integrator']['startTime'] = iter_i * dt
        config['Integrator']['maxIter'] = iter_i_1
        if i > 0:
            config['Flow']['initCase'] = 'Restart'
            config['Flow']['restartDir'] = '%s/job%d/sample%s/fluid_iter%010d' % (args.out_dir, i-1, j, iter_i)
            config['Particles']['initCase'] = 'Restart'
            config['Particles']['restartDir'] = '%s/job%d/sample%s/particles_iter%010d' % (args.out_dir, i-1, j, iter_i)
        with open('job%d/case%d.json' % (i, j), 'w') as fout:
            json.dump(cases[j], fout, indent=4)

print "Testcases created; when you're ready, run:"
for i in range(0, num_jobs):
    after_str = 'AFTER=$JOBID ' if i > 0 else ''
    print 'JOBOUT=`%s$SOLEIL_DIR/src/soleil.sh $(echo -i\ job%d/case{0..%d}.json) -o %s/job%d`; JOBID="${JOBOUT//[!0-9]/}"; echo $JOBID' % (after_str, i, len(cases) - 1, args.out_dir, i)
