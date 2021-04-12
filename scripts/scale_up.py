#!/usr/bin/env python2

import argparse
import json

parser = argparse.ArgumentParser()
parser.add_argument('base_json', type=argparse.FileType('r'))
parser.add_argument('-n', '--num-times', type=int, default=8)
parser.add_argument('-s', '--strong-scale', action='store_true')
args = parser.parse_args()

# Read base config
config = json.load(args.base_json)
x_tiles = int(config['Mapping']['tiles'][0])
y_tiles = int(config['Mapping']['tiles'][1])
z_tiles = int(config['Mapping']['tiles'][2])
x_tpr = int(config['Mapping']['tilesPerRank'][0])
y_tpr = int(config['Mapping']['tilesPerRank'][1])
z_tpr = int(config['Mapping']['tilesPerRank'][2])

# Double the dimension with the fewest tiles first
if x_tpr == min(x_tpr, y_tpr, z_tpr):
    if y_tpr <= z_tpr:
        doubling_order = [0,1,2]
    else:
        doubling_order = [0,2,1]
elif y_tpr == min(x_tpr, y_tpr, z_tpr):
    if x_tpr <= z_tpr:
        doubling_order = [1,0,2]
    else:
        doubling_order = [1,2,0]
else: # z_tpr == min(x_tpr, y_tpr, z_tpr)
    if x_tpr <= y_tpr:
        doubling_order = [2,0,1]
    else:
        doubling_order = [2,1,0]

# Scale down to 1 node
config['Mapping']['tiles'][0] = x_tpr
config['Mapping']['tiles'][1] = y_tpr
config['Mapping']['tiles'][2] = z_tpr
config['Grid']['xNum'] = int(config['Grid']['xNum']) / x_tiles * x_tpr
config['Grid']['yNum'] = int(config['Grid']['yNum']) / y_tiles * y_tpr
config['Grid']['zNum'] = int(config['Grid']['zNum']) / z_tiles * z_tpr
config['Grid']['xWidth'] = float(config['Grid']['xWidth']) / x_tiles * x_tpr
config['Grid']['yWidth'] = float(config['Grid']['yWidth']) / y_tiles * y_tpr
config['Grid']['zWidth'] = float(config['Grid']['zWidth']) / z_tiles * z_tpr
config['Particles']['initNum'] = (int(config['Particles']['initNum'])
                                  / x_tiles / y_tiles / z_tiles
                                  * x_tpr   * y_tpr   * z_tpr)
config['Particles']['maxNum'] = (int(config['Particles']['maxNum'])
                                 / x_tiles / y_tiles / z_tiles
                                 * x_tpr   * y_tpr   * z_tpr)
if config['Radiation']['type'] == 'DOM':
    config['Radiation']['xNum'] = int(config['Radiation']['xNum']) / x_tiles * x_tpr
    config['Radiation']['yNum'] = int(config['Radiation']['yNum']) / y_tiles * y_tpr
    config['Radiation']['zNum'] = int(config['Radiation']['zNum']) / z_tiles * z_tpr
with open('1.json', 'w') as fout:
    json.dump(config, fout, indent=4)

# Scale up
for i in range(0,args.num_times):
    dim = doubling_order[i % 3]
    config['Mapping']['tiles'][dim] *= 2
    if not args.strong_scale:
        if dim == 0:
            config['Grid']['xNum'] *= 2
            config['Grid']['xWidth'] *= 2
            if config['Radiation']['type'] == 'DOM':
                config['Radiation']['xNum'] *= 2
        elif dim == 1:
            config['Grid']['yNum'] *= 2
            config['Grid']['yWidth'] *= 2
            if config['Radiation']['type'] == 'DOM':
                config['Radiation']['yNum'] *= 2
        else: # dim == 2
            config['Grid']['zNum'] *= 2
            config['Grid']['zWidth'] *= 2
            if config['Radiation']['type'] == 'DOM':
                config['Radiation']['zNum'] *= 2
        config['Particles']['initNum'] *= 2
        config['Particles']['maxNum'] *= 2
    with open(str(2**(i+1)) + '.json', 'w') as fout:
        json.dump(config, fout, indent=4)
