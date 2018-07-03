#!/usr/bin/env python

import fileinput
import random

def is_rev_sorted(l):
    return all(l[i] >= l[i+1] for i in range(len(l)-1))

xi = []
eta = []
mu = []
w = []
num_angles = None

for line in fileinput.input():
    line = line[:-1]
    if num_angles is None:
        num_angles = int(line)
    elif len(xi) < num_angles:
        xi.append(line)
    elif len(eta) < num_angles:
        eta.append(line)
    elif len(mu) < num_angles:
        mu.append(line)
    elif len(w) < num_angles:
        w.append(line)
    else:
        assert False

assert len(xi) == num_angles
assert len(eta) == num_angles
assert len(mu) == num_angles
assert len(w) == num_angles

def possible_indices(m):
    indices = []
    if float(xi[m]) >= 0.0:
        if float(eta[m]) >= 0.0:
            if float(mu[m]) >= 0.0:
                indices.append(0)
            if float(mu[m]) <= 0.0:
                indices.append(1)
        if float(eta[m]) <= 0.0:
            if float(mu[m]) >= 0.0:
                indices.append(2)
            if float(mu[m]) <= 0.0:
                indices.append(3)
    if float(xi[m]) <= 0.0:
        if float(eta[m]) >= 0.0:
            if float(mu[m]) >= 0.0:
                indices.append(4)
            if float(mu[m]) <= 0.0:
                indices.append(5)
        if float(eta[m]) <= 0.0:
            if float(mu[m]) >= 0.0:
                indices.append(6)
            if float(mu[m]) <= 0.0:
                indices.append(7)
    assert len(indices) >= 1
    return indices

quadrants = None
while True:
    quadrants = [[],[],[],[],[],[],[],[]]
    def filter_shortest(indices):
        min_len = min([len(quadrants[i]) for i in indices])
        return [i for i in indices if len(quadrants[i]) == min_len]
    for m in range(num_angles):
        indices = possible_indices(m)
        # TODO: Don't have to do this from scratch on every try
        if len(indices) == 1:
            quadrants[indices[0]].append(m)
    for m in range(num_angles):
        indices = possible_indices(m)
        if len(indices) > 1:
            indices = filter_shortest(indices)
            quadrants[random.choice(indices)].append(m)
    quad_sizes = [len(q) for q in quadrants]
    if max(quad_sizes) <= min(quad_sizes) + 1 and is_rev_sorted(quad_sizes):
        break

def print_round_robin(field):
    max_quad_len = max([len(q) for q in quadrants])
    for i in range(max_quad_len):
        for q in quadrants:
            if i < len(q):
                m = q[i]
                print field[m]

print num_angles
print_round_robin(xi)
print_round_robin(eta)
print_round_robin(mu)
print_round_robin(w)
