#!/usr/bin/env python2

import argparse
import json

parser = argparse.ArgumentParser()
parser.add_argument('base_json', type=argparse.FileType('r'))
parser.add_argument('-n', '--num_times', type=int, default=8)
args = parser.parse_args()

# Read base config
config = json.load(args.base_json)
assert int(config['Mapping']['tiles'][0]) / int(config['Mapping']['tilesPerRank'][0]) == 1
assert int(config['Mapping']['tiles'][1]) / int(config['Mapping']['tilesPerRank'][1]) == 1
assert int(config['Mapping']['tiles'][2]) / int(config['Mapping']['tilesPerRank'][2]) == 1

# Scale up
with open('1.json', 'w') as fout:
    json.dump(config, fout, indent=4)
for i in range(0,args.num_times):
    config['Mapping']['tiles'][0] *= 2
    config['Grid']['xNum'] *= 2
    config['Grid']['xWidth'] *= 2
    if config['Radiation']['type'] == 'DOM':
        config['Radiation']['xNum'] *= 2
    config['Particles']['initNum'] *= 2
    config['Particles']['maxNum'] *= 2
    with open(str(2**(i+1)) + '.json', 'w') as fout:
        json.dump(config, fout, indent=4)
