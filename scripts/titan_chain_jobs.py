#!/usr/bin/env python

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
assert(config['IO']['wrtRestart'])
assert(float(config['Integrator']['cfl']) < 0.0)
iters_per_run = int(config['Integrator']['maxIter'])
dt = float(config['Integrator']['fixedDeltaTime'])

json.dump(config, open('0.json', 'w'), indent=4)
for i in range(1, args.num_runs):
    config['Integrator']['restartIter'] = i * iters_per_run
    config['Integrator']['restartTime'] = i * dt
    config['Integrator']['maxIter'] = (i+1) * iters_per_run
    config['Flow']['initCase'] = 'Restart'
    config['Flow']['restartDir'] = '%s/%s/sample0/fluid_iter%010d' % (args.out_dir, i-1, i * iters_per_run)
    config['Particles']['initCase'] = 'Restart'
    config['Particles']['restartDir'] = '%s/%s/sample0/particles_iter%010d' % (args.out_dir, i-1, i * iters_per_run)
    json.dump(config, open(str(i) + '.json', 'w'), indent=4)
    mkdirp('%s/%s' % (args.out_dir, i-1))

print 'QUEUE=batch $SOLEIL_DIR/src/soleil.sh -i 0.json -o %s/0' % args.out_dir
for i in range(1, args.num_runs):
    print 'AFTER=<prev-job> QUEUE=batch $SOLEIL_DIR/src/soleil.sh -i %s.json -o %s/%s' % (i, args.out_dir, i)
