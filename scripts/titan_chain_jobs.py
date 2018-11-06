#!/usr/bin/env python2

import argparse
import json
import os

def mkdirp(directory):
    if not os.path.exists(directory):
        os.makedirs(directory)

parser = argparse.ArgumentParser()
parser.add_argument('num_runs', type=int)
parser.add_argument('out_dir')
parser.add_argument('json_file', nargs='+', type=argparse.FileType('r'))
args = parser.parse_args()

cases = []
for fin in args.json_file:
    case = json.load(fin)
    cases.append(case)
    assert('configs' in case and len(case['configs']) == 2)
    iters = int(case['configs'][0]['Integrator']['maxIter'])
    dt = float(case['configs'][0]['Integrator']['fixedDeltaTime'])
    for config in case['configs']:
        assert(int(config['Integrator']['startIter']) == 0)
        assert(float(config['Integrator']['startTime']) == 0.0)
        assert(int(config['Integrator']['maxIter']) == iters)
        assert(float(config['Integrator']['cfl']) < 0.0)
        assert(float(config['Integrator']['fixedDeltaTime']) == dt)
        assert(config['Flow']['initCase'] != 'Restart')
        assert(config['Particles']['initCase'] != 'Restart')
        assert(config['IO']['wrtRestart'])

for i in range(0, args.num_runs):
    mkdirp('run%s' % i)
    mkdirp('%s/run%s' % (args.out_dir, i))
    for j in range(0, len(cases)):
        if i > 0:
            for k in range(0, 2):
                config = cases[j]['configs'][k]
                iters = int(config['Integrator']['maxIter']) - int(config['Integrator']['startIter'])
                dt = float(config['Integrator']['fixedDeltaTime'])
                config['Integrator']['startIter'] = i * iters
                config['Integrator']['startTime'] = i * iters * dt
                config['Integrator']['maxIter'] = (i+1) * iters
                config['Flow']['initCase'] = 'Restart'
                config['Flow']['restartDir'] = '%s/run%s/sample%s/fluid_iter%010d' % (args.out_dir, i-1, 2*j+k, i*iters)
                config['Particles']['initCase'] = 'Restart'
                config['Particles']['restartDir'] = '%s/run%s/sample%s/particles_iter%010d' % (args.out_dir, i-1, 2*j+k, i*iters)
        with open('run%s/case%s.json' % (i, j), 'w') as fout:
            json.dump(cases[j], fout, indent=4)

print "Testcases created; when you're ready, run:"
for i in range(0, args.num_runs):
    after_str = 'AFTER=$JOBID ' if i > 0 else ''
    cases_str = ' '.join(['-m run%s/case%s.json' % (i,j) for j in range(0,len(cases))])
    print 'JOBID=`%sQUEUE=batch $SOLEIL_DIR/src/soleil.sh %s -o %s/run%s`; echo $JOBID' % (after_str, cases_str, args.out_dir, i)
