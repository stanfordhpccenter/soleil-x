#!/usr/bin/env python2

import math
import numpy as np
import random

np.set_printoptions(linewidth=np.nan)

A = np.genfromtxt('all2.csv', dtype='S20,<i8,<i8,<i8,<i8,<i8,<i8,<i8,<i8,<i8,<i8,<f8,<f8', delimiter=',', names=True)
X = np.vstack((A['flow_x'], A['flow_y'], A['flow_z'], A['parcel_size'], A['dom_x'], A['dom_y'], A['dom_z'], A['quads'], A['rk_order'], A['ftts'])).T
r = A['rho']
M = {}
for i in range(len(X)):
    M[tuple(X[i])] = r[i]

def neighbors(c):
    def swap(v, i, choice1, choice2):
        assert v[i] == choice1 or v[i] == choice2
        v[i] = choice1 if v[i] == choice2 else choice2
    def swapped(i, choice1, choice2):
        v = list(c)
        swap(v, i, choice1, choice2)
        return tuple(v)
    res = []
    res.append(swapped(0, 128, 64))
    res.append(swapped(1, 32, 16))
    res.append(swapped(2, 32, 16))
    res.append(swapped(3, 10, 100))
    if (c[4] == 0  and c[5] == 0 and c[6] == 0 and c[7] == 0  or
        c[4] == 32 and c[5] == 8 and c[6] == 8 and c[7] == 50):
        v = list(c)
        swap(v, 4, 0, 32)
        swap(v, 5, 0, 8)
        swap(v, 6, 0, 8)
        swap(v, 7, 0, 50)
        res.append(tuple(v))
    if c[4] != 0:
        res.append(swapped(4, 64, 32))
        res.append(swapped(5, 16, 8))
        res.append(swapped(6, 16, 8))
        res.append(swapped(7, 86, 50))
    res.append(swapped(8, 4, 3))
    res.append(swapped(9, 3, 2))
    return res

# =============================================================================
# Hill climbing
# =============================================================================

# def climb(c):
#     visited = set()
#     visited.add(c)
#     # print 'LF =', c, 'rho =', M[c], 'visited =', len(visited)
#     while True:
#         ns = [n for n in neighbors(c) if n not in visited]
#         if len(ns) == 0:
#             break
#         visited.update(ns)
#         best = ns[ np.argmax([M[n] for n in ns]) ]
#         if M[best] <= M[c]:
#             break
#         c = best
#         # print 'LF =', c, 'rho =', M[c], 'visited =', len(visited)
#     # print
#     return c, len(visited)

# random.seed()
# # c = (128,32,32,10,64,16,16,86,4,3)
# max_corr = max(r)
# full_work = len(X)
# num_runs = 100
# corr_loss = []
# work_saved = []
# for dummy in range(num_runs):
#     c, work = climb(random.choice(M.keys()))
#     corr_loss.append(max_corr - M[c])
#     work_saved.append((full_work-work)/float(full_work))
# print np.mean(corr_loss), np.mean(work_saved)

# =============================================================================
# MCMC
# =============================================================================

def a(b, new, prev):
    return min(1.0, math.exp(b*(new-prev)))

def search(b, max_tries, search_space_size):
    visited = set()
    tradeoff = []
    c = random.choice(M.keys())
    best = M[c]
    visited.add(c)
    tradeoff.append(0)
    # print 'step =', -1, 'best =', best, 'visited =', len(visited)
    for step in range(max_tries):
        n = random.choice(neighbors(c))
        visited.add(n)
        if M[n] > best:
            for i in range(len(tradeoff),len(visited)):
                tradeoff.append(best)
            best = M[n]
            # print 'step =', step, 'best =', best, 'visited =', len(visited)
        if random.random() < a(b, M[n], M[c]):
            c = n
    for i in range(len(tradeoff),search_space_size+1):
        tradeoff.append(best)
    # print
    return tradeoff

random.seed()
max_corr = max(r)
full_work = len(X)
num_runs = 100
acc = np.zeros(full_work + 1)
for dummy in range(num_runs):
    acc += search(1.0, 10000, full_work)
avg_tradeoff = acc / float(num_runs)
print 'work_saved\tcorr_loss'
for (work,corr) in enumerate(avg_tradeoff):
    work_saved = (full_work-work)/float(full_work)
    corr_loss = max_corr - corr
    print '{:f}\t{:f}'.format(work_saved, corr_loss)
