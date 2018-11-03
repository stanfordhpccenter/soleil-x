#!/usr/bin/env python

import json

NUM_TIMES=8

with open('1.json') as fin:
    config = json.load(fin)
    assert config['Mapping']['tiles'][0] == 1
    assert config['Mapping']['tiles'][1] == 1
    assert config['Mapping']['tiles'][2] == 1
    for i in range(0,NUM_TIMES):
        if i % 3 == 0:
            config['Mapping']['tiles'][2] *= 2
            config['Grid']['zNum'] *= 2
            config['Grid']['zWidth'] *= 2
            config['Radiation']['zNum'] *= 2
        elif i % 3 == 1:
            config['Mapping']['tiles'][1] *= 2
            config['Grid']['yNum'] *= 2
            config['Grid']['yWidth'] *= 2
            config['Radiation']['yNum'] *= 2
        else: # i % 3 == 2
            config['Mapping']['tiles'][0] *= 2
            config['Grid']['xNum'] *= 2
            config['Grid']['xWidth'] *= 2
            config['Radiation']['xNum'] *= 2
        config['Particles']['initNum'] *= 2
        config['Particles']['maxNum'] *= 2
        with open(str(2**(i+1)) + '.json', 'w') as fout:
            json.dump(config, fout, indent=4)
