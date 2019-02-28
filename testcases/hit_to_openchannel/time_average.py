#!/usr/bin/env python2

import argparse
import json
import math

parser = argparse.ArgumentParser()
parser.add_argument('case_json', type=argparse.FileType('r'))
parser.add_argument('probe_file', type=argparse.FileType('r'))
args = parser.parse_args()

mc = json.load(args.case_json)
L = 0.16
U_0 = float(mc['configs'][1]['Flow']['initParams'][2])
delta_t_c = float(mc['configs'][0]['Integrator']['fixedDeltaTime'])
ftt_iters = int(math.ceil(L / U_0 / delta_t_c))
max_iter = int(mc['configs'][0]['Integrator']['maxIter'])

AvgFluidT = 0.0
AvgParticleT = 0.0
AvgCellOfParticleT = 0.0

next(args.probe_file)
for line in args.probe_file:
    line = line[:-1]
    toks = line.split()
    assert(len(toks) == 4)
    curr_iter = int(toks[0])
    if curr_iter < ftt_iters:
        continue
    AvgFluidT += float(toks[1])
    AvgParticleT += float(toks[2])
    AvgCellOfParticleT += float(toks[3])

AvgFluidT /= (max_iter - ftt_iters)
AvgParticleT /= (max_iter - ftt_iters)
AvgCellOfParticleT /= (max_iter - ftt_iters)
        
print '%s\t%s\t%s' % (AvgFluidT, AvgParticleT, AvgCellOfParticleT)
