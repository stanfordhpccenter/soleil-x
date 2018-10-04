#!/usr/bin/env python2

# TODO: Adapt this script to run on output of -fpretty

import fileinput
import re

privs = {} # map(task:string,privileges:string)
curr_task = None
in_privs = False
in_main = False

for line in fileinput.input():
    line = line[:-1]
    if line.startswith('task') or line.startswith('local task'):
        assert(curr_task is None and not in_privs)
        task_name_idx = 1 if line.startswith('task') else 2
        curr_task = line.split()[task_name_idx].split('(')[0]
        assert(curr_task not in privs)
        privs[curr_task] = ''
    elif line.startswith('where'):
        assert(curr_task is not None and not in_privs)
        in_privs = True
    elif line.startswith('do'):
        assert(curr_task is not None and in_privs)
        in_privs = False
    elif line.startswith('end'):
        assert(not in_privs)
        curr_task = None
    elif in_privs:
        privs[curr_task] += ('\n' + line)
    elif 'MAIN' in line:
        in_main = True
    elif in_main:
        for task in privs:
            if (task+'(') in line and privs[task] != '':
                print task + privs[task]
