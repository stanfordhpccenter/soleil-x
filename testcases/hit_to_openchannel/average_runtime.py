#!/usr/bin/env python2

import argparse

parser = argparse.ArgumentParser()
parser.add_argument('-n', type=int)
parser.add_argument('console_file', type=argparse.FileType('r'), nargs='+')
args = parser.parse_args()
if args.n is None:
    args.n = len(args.console_file)

total = 0
for f in args.console_file:
    for line in f:
        pass
    total += float(line.split()[2])
print total / args.n
