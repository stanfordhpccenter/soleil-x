#!/usr/bin/env python3

import argparse
import numpy as np

parser = argparse.ArgumentParser()
parser.add_argument('tiles', type=int, nargs=3)
parser.add_argument('-v', '--verbose', action='store_true')
args = parser.parse_args()
X = args.tiles[0]
Y = args.tiles[1]
Z = args.tiles[2]

quads = [1,2,3,4,5,6,7,8]

state = {q:np.zeros((X,Y,Z), dtype=np.int8) for q in quads}
state[1][  0,  0,  0] = 1
state[2][  0,  0,Z-1] = 1
state[3][  0,Y-1,  0] = 1
state[4][  0,Y-1,Z-1] = 1
state[5][X-1,  0,  0] = 1
state[6][X-1,  0,Z-1] = 1
state[7][X-1,Y-1,  0] = 1
state[8][X-1,Y-1,Z-1] = 1

curr = np.zeros((X,Y,Z), dtype=np.int8)

t = 0
while True:
    done = True
    for x in range(X):
        for y in range(Y):
            for z in range(Z):
                def downstream_diagonals(q):
                    if state[q][x,y,z] != 1:
                        return -1
                    return (((X-x-1) if q in [1,2,3,4] else x) +
                            ((Y-y-1) if q in [1,2,5,6] else y) +
                            ((Z-z-1) if q in [1,3,5,7] else z))
                q = max(quads, key=downstream_diagonals)
                if state[q][x,y,z] != 1:
                    continue
                done = False
                state[q][x,y,z] = 2
                curr[x,y,z] = q
    if done:
        break
    t += 1
    if args.verbose:
        print()
        print('== %s' % t)
    for x in range(X):
        if args.verbose:
            print()
        for y in range(Y):
            for z in range(Z):
                q = curr[x,y,z]
                if args.verbose:
                    print('.' if q == 0 else str(q) , end='')
                if q == 0:
                    continue
                xx = (x+1) if q in [1,2,3,4] else (x-1)
                yy = (y+1) if q in [1,2,5,6] else (y-1)
                zz = (z+1) if q in [1,3,5,7] else (z-1)
                if 0 <= xx and xx < X:
                    state[q][xx, y, z] = max(1, state[q][xx, y, z])
                if 0 <= yy and yy < Y:
                    state[q][ x,yy, z] = max(1, state[q][ x,yy, z])
                if 0 <= zz and zz < Z:
                    state[q][ x, y,zz] = max(1, state[q][ x, y,zz])
                curr[x,y,z] = 0
            if args.verbose:
                print()
print(t)
