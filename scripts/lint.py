#!/usr/bin/env python2

from itertools import count
import fileinput
import re

task_name = None
task_start = -1
task_depth = -1
num_tasks = 0
demands = []
annots = []
region_args = 0

def record_annot(a):
    if a in line:
        annots.append(a)

in_openmp = False
openmp_start = -1
openmp_depth = -1

def check_task():
    def warn(msg):
        print 'Warning: Line %s: %s' % (task_start, msg)
    if '__inline' in demands:
        return
    if task_name.startswith('work') or task_name == 'main':
        if '__inner' not in demands:
            warn('Top-level task not declared inner')
        return
    if '__leaf' not in demands and 'NOT LEAF' not in annots:
        warn('Task not declared leaf')
    if '__parallel' in demands and region_args != 1:
        warn('Auto-parallelized tasks should (usually) have exactly 1 region argument')
    if '__parallel' not in demands and 'MANUALLY PARALLELIZED' not in annots:
        warn('Task not auto-parallelized')
    if '__cuda' not in demands and 'NO CUDA' not in annots:
        warn('Task has no GPU variant')
    if '__openmp' not in demands and 'NO OPENMP' not in annots:
        warn('Task has no OpenMP loop')

# TODO: Warn if launching a non-inlined task
def check_in_openmp(line, lineno):
    line = line.split('--')[0]
    def warn(msg):
        print 'Warning: Line %s: %s' % (lineno, msg)
    if 'rand(' in line or 'drand48(' in line or 'drand48_r(' in line:
        warn('Calling random function inside OpenMP loop')
    if ('C.' in line or 'regentlib' in line) and '__cuda' in demands:
        warn('Calling external function inside potentially CUDA-ized loop')
    if 'config' in line and '__cuda' in demands:
        warn('Accessing struct inside potentially CUDA-ized loop')

def check_line(line, lineno):
    line = line.split('--')[0]
    def warn(msg):
        print 'Warning: Line %s: %s' % (lineno, msg)
    if 'rand(' in line or 'drand48(' in line:
        warn('Calling non-reentrant random function')

for (line, lineno) in zip(fileinput.input(), count(start=1)):
    line = line[:-1]
    # Record linter annotations
    record_annot('MANUALLY PARALLELIZED')
    record_annot('NO CUDA')
    record_annot('NO OPENMP')
    record_annot('NOT LEAF')
    # Record region arguments
    for m in re.finditer(r':\s*region', line):
        region_args += 1
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
    elif (line.strip().startswith('task') or
          line.strip().startswith('local task')):
        assert(task_name is None and not in_openmp)
        if line.strip().startswith('task'):
            task_name_idx = 1
            task_depth = line.find('t')
        else:
            task_name_idx = 2
            task_depth = line.find('l')
        task_name = line.strip().split()[task_name_idx].split('(')[0]
        task_start = lineno
    # Process end of task
    elif (line.strip().startswith('end')
          and task_name is not None
          and line.find('e') == task_depth):
        assert(not in_openmp)
        check_task()
        task_name = None
        task_start = -1
        task_depth = -1
        num_tasks += 1
        demands = []
        annots = []
        region_args = 0
    # Check statements inside OpenMP loops
    elif in_openmp:
        check_in_openmp(line, lineno)
    # Check statements
    check_line(line, lineno)

print 'Total: %s task(s)' % num_tasks
