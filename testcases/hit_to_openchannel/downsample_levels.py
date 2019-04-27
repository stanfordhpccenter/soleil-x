#!/usr/bin/env python2

import argparse

parser = argparse.ArgumentParser()
parser.add_argument('max_ftt', type=int)
parser.add_argument('levels_dat', type=argparse.FileType('r'))
args = parser.parse_args()

past_header = False
for line in args.levels_dat:
    toks = line.split()
    if not past_header:
        past_header = True
    else:
        orig_ftts = int(toks[-1])
        assert orig_ftts > args.max_ftt
        toks[-1] = str(args.max_ftt)
    print ' '.join(toks)
