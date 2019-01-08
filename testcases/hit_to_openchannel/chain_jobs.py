#!/usr/bin/env python2

import argparse
import json
import math
import os

parser = argparse.ArgumentParser()
parser.add_argument('out_dir')
parser.add_argument('json_file', nargs='+', type=argparse.FileType('r'))
args = parser.parse_args()

iters_per_job = None
target_iter = 0
cases = []
for fin in args.json_file:
    case = json.load(fin)
    cases.append(case)
    assert('configs' in case and len(case['configs']) == 2)
    if iters_per_job is None:
        iters_per_job = int(case['configs'][0]['Integrator']['maxIter'])
    dt = float(case['configs'][0]['Integrator']['fixedDeltaTime'])
    ftts = int(case['flowThroughTimes'])
    u_0 = float(case['configs'][1]['BC']['xBCLeftInflowProfile']['addedVelocity'])
    length = float(case['configs'][1]['Grid']['xWidth'])
    target_iter = max(target_iter, int(math.ceil(ftts * length / u_0 / dt)))
    for config in case['configs']:
        assert(int(config['Integrator']['startIter']) == 0)
        assert(float(config['Integrator']['startTime']) == 0.0)
        assert(int(config['Integrator']['maxIter']) == iters_per_job)
        assert(float(config['Integrator']['cfl']) < 0.0)
        assert(float(config['Integrator']['fixedDeltaTime']) == dt)
        assert(config['Flow']['initCase'] != 'Restart')
        assert(config['Particles']['initCase'] != 'Restart')
        assert(bool(config['IO']['wrtRestart']))

num_jobs = int(math.ceil(float(target_iter) / iters_per_job))
for i in range(0, num_jobs):
    os.mkdir('job%04d' % i)
    os.mkdir('%s/job%04d' % (args.out_dir, i))
    for j in range(0, len(cases)):
        if i > 0:
            for k in range(0, 2):
                config = cases[j]['configs'][k]
                dt = float(config['Integrator']['fixedDeltaTime'])
                config['Integrator']['startIter'] = i * iters_per_job
                config['Integrator']['startTime'] = i * iters_per_job * dt
                config['Integrator']['maxIter'] = (i+1) * iters_per_job
                config['Flow']['initCase'] = 'Restart'
                config['Flow']['restartDir'] = '%s/job%04d/sample%s/fluid_iter%010d' % (args.out_dir, i-1, 2*j+k, i*iters_per_job)
                config['Particles']['initCase'] = 'Restart'
                config['Particles']['restartDir'] = '%s/job%04d/sample%s/particles_iter%010d' % (args.out_dir, i-1, 2*j+k, i*iters_per_job)
        with open('job%04d/case%04d.json' % (i, j), 'w') as fout:
            json.dump(cases[j], fout, indent=4)

print "Testcases created; when you're ready, run:"
for i in range(0, num_jobs):
    after_str = 'AFTER=$JOBID ' if i > 0 else ''
    cases_str = ' '.join(['-m job%04d/case%04d.json' % (i,j) for j in range(0,len(cases))])
    print 'JOBID=`%sQUEUE=batch $SOLEIL_DIR/src/soleil.sh %s -o %s/job%04d`; echo $JOBID' % (after_str, cases_str, args.out_dir, i)
