#!/usr/bin/env python2

import argparse
import json
import sys

parser = argparse.ArgumentParser()
parser.add_argument('hf_json', type=argparse.FileType('r'))
parser.add_argument('flow_x', type=int)
parser.add_argument('flow_y', type=int)
parser.add_argument('flow_z', type=int)
parser.add_argument('parcel_size', type=int)
parser.add_argument('dom_x', type=int)
parser.add_argument('dom_y', type=int)
parser.add_argument('dom_z', type=int)
parser.add_argument('quads', type=int)
parser.add_argument('use_dom', type=bool)
parser.add_argument('rk_order', type=int)
parser.add_argument('ftts', type=int)
args = parser.parse_args()

assert(args.flow_x % 4 == 0)
assert(args.flow_x >= args.dom_x and args.flow_x % args.dom_x == 0)
assert(args.flow_y >= args.dom_y and args.flow_y % args.dom_y == 0)
assert(args.flow_z >= args.dom_z and args.flow_z % args.dom_z == 0)
assert(args.ftts >= 2, 'At least one transient and one averaging FTT required')

# Parse json template
mc = json.load(args.hf_json)

# Set HIT section values
mc['configs'][0]['Grid']['xNum'] = args.flow_x / 4
mc['configs'][0]['Grid']['yNum'] = args.flow_y
mc['configs'][0]['Grid']['zNum'] = args.flow_z
mc['configs'][0]['Integrator']['rkOrder'] = args.rk_order
mc['configs'][0]['Particles']['parcelSize'] = args.parcel_size

# Set channel section values
mc['configs'][1]['Grid']['xNum'] = args.flow_x
mc['configs'][1]['Grid']['yNum'] = args.flow_y
mc['configs'][1]['Grid']['zNum'] = args.flow_z
mc['configs'][1]['Integrator']['rkOrder'] = args.rk_order
mc['configs'][1]['Particles']['parcelSize'] = args.parcel_size
if args.use_dom:
    mc['configs'][1]['Radiation']['xNum'] = args.dom_x
    mc['configs'][1]['Radiation']['yNum'] = args.dom_y
    mc['configs'][1]['Radiation']['zNum'] = args.dom_z
    mc['configs'][1]['Radiation']['angles'] = args.quads
else:
    mc['configs'][1]['Radiation'] = {}
    mc['configs'][1]['Radiation']['type'] = 'Algebraic'
    mc['configs'][1]['Radiation']['intensity'] = 'TBD'
    mc['configs'][1]['Radiation']['absorptivity'] = 'TBD'
mc['flowThroughTimes'] = args.ftts

# Update grid-related values
mc['configs'][1]['IO']['probes'][0]['fromCell'][0] = args.flow_x
mc['configs'][1]['IO']['probes'][0]['fromCell'][1] = 0
mc['configs'][1]['IO']['probes'][0]['fromCell'][2] = 0
mc['configs'][1]['IO']['probes'][0]['uptoCell'][0] = args.flow_x
mc['configs'][1]['IO']['probes'][0]['uptoCell'][1] = args.flow_y - 1
mc['configs'][1]['IO']['probes'][0]['uptoCell'][2] = args.flow_z - 1
mc['copySrc']['uptoCell'][1] = args.flow_y - 1
mc['copySrc']['uptoCell'][2] = args.flow_z - 1
mc['copyTgt']['uptoCell'][1] = args.flow_y - 1
mc['copyTgt']['uptoCell'][2] = args.flow_z - 1

# Dump final json config
json.dump(mc, sys.stdout, indent=4)
