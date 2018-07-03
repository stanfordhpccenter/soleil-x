#!/usr/bin/env python

import argparse

parser = argparse.ArgumentParser()
parser.add_argument('f1')
parser.add_argument('f2')
args = parser.parse_args()

def parse_angle_file(f):
    angles = []
    num_angles = None
    xi_idx = 0
    eta_idx = 0
    mu_idx = 0
    w_idx = 0
    for line in f:
        line = line[:-1]
        if num_angles is None:
            num_angles = int(line)
            for i in range(num_angles):
                angles.append([None,None,None,None])
        elif xi_idx < num_angles:
            angles[xi_idx][0] = line
            xi_idx += 1
        elif eta_idx < num_angles:
            angles[eta_idx][1] = line
            eta_idx += 1
        elif mu_idx < num_angles:
            angles[mu_idx][2] = line
            mu_idx += 1
        elif w_idx < num_angles:
            angles[w_idx][3] = line
            w_idx += 1
        else:
            assert False
    assert xi_idx == num_angles
    assert eta_idx == num_angles
    assert mu_idx == num_angles
    assert w_idx == num_angles
    return angles

with open(args.f1) as f1:
    a1 = parse_angle_file(f1)
    a1.sort()
with open(args.f2) as f2:
    a2 = parse_angle_file(f2)
    a2.sort()
if len(a1) != len(a2):
    print 'DIFFERENT LENGTH'
elif not all([x == y for (x,y) in zip(a1,a2)]):
    print 'DIFFERENT CONTENTS'
else:
    print 'SAME'
