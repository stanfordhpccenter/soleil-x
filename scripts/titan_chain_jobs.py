#!/usr/bin/env python2

import argparse
import json
import os

def mkdirp(directory):
    if not os.path.exists(directory):
        os.makedirs(directory)

parser = argparse.ArgumentParser()
parser.add_argument('json_file')
parser.add_argument('out_dir')
parser.add_argument('num_runs', type=int)
args = parser.parse_args()

config = json.load(open(args.json_file))
samples = config['configs'] if 'configs' in config else [config]
iters_per_run = int(samples[0]['Integrator']['maxIter'])
dt = float(samples[0]['Integrator']['fixedDeltaTime'])
for j in range(0, len(samples)):
    assert(samples[j]['IO']['wrtRestart'])
    assert(float(samples[j]['Integrator']['cfl']) < 0.0)
    assert(iters_per_run == int(samples[j]['Integrator']['maxIter']))
    assert(dt == float(samples[0]['Integrator']['fixedDeltaTime']))

json.dump(config, open('0.json', 'w'), indent=4)
for i in range(1, args.num_runs):
    for j in range(0, len(samples)):
       samples[j]['Integrator']['startIter'] = i * iters_per_run
       samples[j]['Integrator']['startTime'] = i * iters_per_run * dt
       samples[j]['Integrator']['maxIter'] = (i+1) * iters_per_run
       samples[j]['Flow']['initCase'] = 'Restart'
       samples[j]['Flow']['restartDir'] = '%s/%s/sample%s/fluid_iter%010d' % (args.out_dir, i-1, j, i * iters_per_run)
       samples[j]['Particles']['initCase'] = 'Restart'
       samples[j]['Particles']['restartDir'] = '%s/%s/sample%s/particles_iter%010d' % (args.out_dir, i-1, j, i * iters_per_run)
    json.dump(config, open(str(i) + '.json', 'w'), indent=4)
    mkdirp('%s/%s' % (args.out_dir, i-1))

switch = '-i' if len(samples) == 1 else '-m'
print 'Testcases created in current directory; when ready, run:'
print 'JOBID=`QUEUE=batch $SOLEIL_DIR/src/soleil.sh %s 0.json -o %s/0`; echo $JOBID' % (switch, args.out_dir)
for i in range(1, args.num_runs):
    print 'JOBID=`AFTER=$JOBID QUEUE=batch $SOLEIL_DIR/src/soleil.sh %s %s.json -o %s/%s`; echo $JOBID' % (switch, i, args.out_dir, i)
