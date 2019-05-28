#!/usr/bin/env python2

import argparse
import itertools

parser = argparse.ArgumentParser()
parser.add_argument('flow_x')
parser.add_argument('flow_y')
parser.add_argument('flow_z')
parser.add_argument('parcel_size')
parser.add_argument('dom_x')
parser.add_argument('dom_y')
parser.add_argument('dom_z')
parser.add_argument('quads')
parser.add_argument('use_dom')
parser.add_argument('rk_order')
parser.add_argument('ftts')
args = parser.parse_args()

print 'flow_x flow_y flow_z parcel_size dom_x dom_y dom_z quads use_dom rk_order ftts'
for (flow_x,
     flow_y,
     flow_z,
     parcel_size,
     use_dom,
     rk_order,
     ftts) in itertools.product(
         args.flow_x.split(','),
         args.flow_y.split(','),
         args.flow_z.split(','),
         args.parcel_size.split(','),
         args.use_dom.split(','),
         args.rk_order.split(','),
         args.ftts.split(',')):
    if use_dom == 'false':
        print '%6s %6s %6s %11s %5s %5s %5s %5s %7s %8s %4s' % (flow_x, flow_y, flow_z, parcel_size, -1, -1, -1, -1, use_dom, rk_order, ftts)
        continue
    assert use_dom == 'true'
    for (dom_x,
         dom_y,
         dom_z,
         quads) in itertools.product(
             args.dom_x.split(','),
             args.dom_y.split(','),
             args.dom_z.split(','),
             args.quads.split(',')):
        print '%6s %6s %6s %11s %5s %5s %5s %5s %7s %8s %4s' % (flow_x, flow_y, flow_z, parcel_size, dom_x, dom_y, dom_z, quads, use_dom, rk_order, ftts)
