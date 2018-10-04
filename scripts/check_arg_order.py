#!/usr/bin/env python2

import fileinput
from itertools import count

formal_args = {}
curr_task = None
curr_args = []
in_sig = False
open_parens = 0

for (line,lineno) in zip(fileinput.input(),count(1)):
    line = line[:-1]
    toks = line.split()
    if line.startswith('task') or line.startswith('local task'):
        assert(curr_task is None)
        task_name_idx = 1 if line.startswith('task') else 2
        curr_task = toks[task_name_idx].split('(')[0]
        toks[task_name_idx] = toks[task_name_idx].split('(')[1]
        toks = toks[task_name_idx:]
        in_sig = True
    else:
        for task in formal_args:
            if (task + '(') in line:
                curr_task = task
                in_sig = False
    if curr_task is None:
        continue
    if in_sig:
        for i in xrange(0,len(toks)):
            if toks[i] == ':' and i > 0:
                curr_args.append(toks[i-1])
    else:
        before_after_call = line.split('(')
        if len(before_after_call) > 2:
            print 'Line %s: Nested call' % lineno
            curr_task = None
            curr_args = []
            open_parens = 0
            continue
        params = [p.strip().strip(')')
                  for p in before_after_call[-1].split(',')]
        if params[-1] == '':
            params = params[:-1]
        curr_args.extend(params)
    open_parens += line.count('(')
    open_parens -= line.count(')')
    if open_parens < 0:
        print lineno
        assert(False)
    if open_parens != 0:
        continue
    if in_sig:
        formal_args[curr_task] = curr_args
    else:
        task = curr_task
        formals = formal_args[curr_task]
        actuals = curr_args
    curr_task = None
    curr_args = []
    if in_sig:
        continue
    if len(formals) != len(actuals):
        print('Line %s: Task %s: %s formals vs %s actuals' %
              (lineno, task, len(formals), len(actuals)))
    else:
        for (f,a) in zip(formals, actuals):
            if f != a:
                print('Line %s: Task %s: argument mismatch: %s vs %s' %
                      (lineno, task, f, a))

if open_parens != 0:
    print('Unbalanced parens')
