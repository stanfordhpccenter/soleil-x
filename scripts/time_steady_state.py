#!/usr/bin/env python2

import argparse
import collections

parser = argparse.ArgumentParser()
parser.add_argument('-a', '--after', type=int, default=5)
parser.add_argument('-b', '--before', type=int, default=5)
parser.add_argument('console', type=argparse.FileType('r'))
args = parser.parse_args()

lineno = 0
t_start = None
tail = collections.deque()
for line in args.console:
    lineno += 1
    if lineno - 2 < args.before:
        continue
    t = float(line.split()[2])
    if t_start is None:
        t_start = t
    if len(tail) == args.after + 1:
        tail.popleft()
    tail.append(t)
assert t_start is not None
assert len(tail) == args.after + 1
t_end = tail[0]
print t_end - t_start
