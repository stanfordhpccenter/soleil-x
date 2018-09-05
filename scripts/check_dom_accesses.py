#!/usr/bin/env python

# Takes a log of accesses made by a DOM sweep and verifies that they were made
# in the correct order. Works with any tiling.
# Make sure you do the following:
# * run on one node
# * run on a single CPU (not OpenMP)
# * print the accesses in the format '(x,y,z)', one per line
# * print the values for one quadrant only
# * print the values for a single sweep
# * don't print multiple times per angle
# * run with -lg:inorder

import argparse
import numpy
import re

parser = argparse.ArgumentParser()
parser.add_argument('Nx', type=int)
parser.add_argument('Ny', type=int)
parser.add_argument('Nz', type=int)
parser.add_argument('quadrant', type=int)
parser.add_argument('accesses_log',  type=argparse.FileType('r'))
args = parser.parse_args()

accessed = numpy.zeros((args.Nx,args.Ny,args.Nz), numpy.int8)

count = 0
for line in args.accesses_log:
    line = line[:-1]
    m = re.match(r'\(([0-9]+),([0-9]+),([0-9]+)\)', line)
    assert m is not None, line
    x = int(m.group(1))
    if ((args.quadrant - 1) >> 2) % 2 == 1:
        x = args.Nx - x - 1
    y = int(m.group(2))
    if ((args.quadrant - 1) >> 1) % 2 == 1:
        y = args.Ny - y - 1
    z = int(m.group(3))
    if ((args.quadrant - 1) >> 0) % 2 == 1:
        z = args.Nz - z - 1
    assert accessed[(x,y,z)] == 0, line
    assert x == 0 or accessed[(x-1,y,  z  )] > 0, line
    assert y == 0 or accessed[(x,  y-1,z  )] > 0, line
    assert z == 0 or accessed[(x,  y,  z-1)] > 0, line
    accessed[(x,y,z)] = 1
    count += 1

assert count == args.Nx * args.Ny * args.Nz
