#!/usr/bin/env python

import fileinput
import re

fluid_reads = {} # map(task:string,list(field:string))
fluid_writes = {} # map(task:string,list(field:string))
curr_task = None
in_privs = False
in_main = False

for line in fileinput.input():
    line = line[:-1]
    if line.startswith('task') or line.startswith('local task'):
        assert(curr_task is None and not in_privs)
        task_name_idx = 1 if line.startswith('task') else 2
        curr_task = line.split()[task_name_idx].split('(')[0]
        assert(curr_task not in fluid_reads)
        fluid_reads[curr_task] = []
        assert(curr_task not in fluid_writes)
        fluid_writes[curr_task] = []
    if line.startswith('where'):
        assert(curr_task is not None and not in_privs)
        in_privs = True
    if line.startswith('do'):
        assert(curr_task is not None and in_privs)
        in_privs = False
    if line.startswith('end'):
        assert(not in_privs)
        curr_task = None
    if in_privs and 'Fluid' in line:
        flds = []
        for m in re.finditer(r'Fluid\.{([^}]*)}', line):
            for f in m.group(1).split(','):
                flds.append(f.strip())
        for m in re.finditer(r'Fluid\.(\w+)', line):
            flds.append(m.group(1))
        if 'reads' in line:
            fluid_reads[curr_task].extend(flds)
        if 'writes' in line:
            fluid_writes[curr_task].extend(flds)
    if 'MAIN' in line:
        in_main = True
    if in_main:
        for task in fluid_reads:
            if (task+'(') not in line:
                continue
            if len(fluid_reads[task]) == 0 and len(fluid_writes[task]) == 0:
                continue
            print task
            if len(fluid_reads[task]) > 0:
                print '  R: ' + ', '.join(fluid_reads[task])
            if len(fluid_writes[task]) > 0:
                print '  W: ' + ', '.join(fluid_writes[task])
