#!/usr/bin/env python

from itertools import count
import fileinput
import re

task_name = None
task_start = -1
num_tasks = 0
demands = []
annots = []

def record_annot(a):
    if a in line:
        annots.append(a)

in_openmp = False
openmp_start = -1
openmp_depth = -1

def check_task():
    def warn(msg):
        print 'Warning: Line %s: %s' % (task_start, msg)
    if '__inline' in demands or task_name == 'work' or task_name == 'main':
        return
    if '__parallel' not in demands and 'MANUALLY PARALLELIZED' not in annots:
        warn('Task not auto-parallelized')
    if '__cuda' not in demands and 'NO CUDA' not in annots:
        warn('Task has no GPU variant')
    if '__openmp' not in demands and 'NO OPENMP' not in annots:
        warn('Task has no OpenMP loop')

def check_in_openmp(line, lineno):
    def warn(msg):
        print 'Warning: Line %s: %s' % (lineno, msg)
    if 'rand(' in line or 'drand48(' in line or 'drand48_r(' in line:
        warn('Calling random function inside OpenMP loop')
    if 'C.' in line and '__cuda' in demands:
        warn('Calling external function inside potentially CUDA-ized loop')
    if 'config' in line and '__cuda' in demands:
        warn('Accessing struct inside potentially CUDA-ized loop')

def check_line(line, lineno):
    def warn(msg):
        print 'Warning: Line %s: %s' % (lineno, msg)
    if 'rand(' in line or 'drand48(' in line:
        warn('Calling non-reentrant random function')

for (line, lineno) in zip(fileinput.input(), count(start=1)) :
    line = line[:-1]
    # Record linter annotations
    record_annot('MANUALLY PARALLELIZED')
    record_annot('NO CUDA')
    record_annot('NO OPENMP')
    # Record __demand annotations
    for m in re.findall(r'__demand\(([^)]*)\)', line):
        demands.extend([a.strip() for a in m.split(',')])
    # Process start of openmp loop
    if line.strip().startswith('__demand(__openmp)'):
        assert(task_name is not None and not in_openmp)
        in_openmp = True
        openmp_start = lineno
        openmp_depth = line.find('_')
    # Process end of openmp loop
    elif (line.strip().startswith('end')
          and in_openmp
          and line.find('e') == openmp_depth):
        in_openmp = False
        openmp_start = -1
        openmp_depth = -1
    # Process start of task
    elif line.startswith('task') or line.startswith('local task'):
        assert(task_name is None and not in_openmp)
        task_name_idx = 1 if line.startswith('task') else 2
        task_name = line.split()[task_name_idx].split('(')[0]
        task_start = lineno
    # Process end of task
    elif line.startswith('end') and task_name is not None:
        assert(not in_openmp)
        check_task()
        task_name = None
        task_start = -1
        num_tasks += 1
        demands = []
        annots = []
    # Check statements inside OpenMP loops
    elif in_openmp:
        check_in_openmp(line, lineno)
    # Check statements
    check_line(line, lineno)

print 'Total: %s task(s)' % num_tasks
