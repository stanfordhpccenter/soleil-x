#!/usr/bin/env python2

import argparse
import collections

parser = argparse.ArgumentParser()
parser.add_argument('-a', '--after', type=int, default=5)
parser.add_argument('-b', '--before', type=int, default=5)
parser.add_argument('console', type=argparse.FileType('r'))
args = parser.parse_args()

i_start = None
t_start = None
tail = collections.deque()
next(args.console)
for line in args.console:
    toks = line.split()
    i = int(toks[0])
    assert i >= 0
    if i < args.before:
        continue
    t = float(toks[2])
    if i_start is None:
        i_start = i
        t_start = t
    else:
        assert i > i_start and t >= t_start
    if len(tail) == args.after + 1:
        tail.popleft()
    tail.append((i,t))
assert t_start is not None
assert len(tail) == args.after + 1
(i_end,t_end) = tail[0]
delta_i = i_end - i_start
delta_t = t_end - t_start
print '%d iterations in %f seconds: %f iters/s' % (delta_i, delta_t, delta_i/delta_t)
