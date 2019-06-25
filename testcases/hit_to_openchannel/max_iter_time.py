#!/usr/bin/env python2

import fileinput

prev_ts = None
max_dt = 0
for line in fileinput.input():
    line = line[:-1]
    curr_ts = line.split()[2]
    if curr_ts == 'Time':
        continue
    curr_ts = float(curr_ts)
    if prev_ts is not None:
        max_dt = max(curr_ts - prev_ts, max_dt)
    prev_ts = curr_ts
print max_dt    
