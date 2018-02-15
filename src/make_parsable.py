#!/usr/bin/env python

import fileinput
import os
import re

HDF_LIBNAME = os.environ['HDF_LIBNAME']
OBJNAME = os.environ['OBJNAME']
USE_HDF = os.environ['USE_HDF'] != '0'

STDLIB_FUNS = ['acos', 'asin', 'atan', 'cbrt', 'ceil', 'exit', 'fabs',
               'fclose', 'fmod', 'fopen', 'fread', 'free', 'fscanf',
               'fseek', 'ftell', 'malloc', 'memcpy', 'pow', 'printf',
               'rand', 'snprintf', 'sqrt', 'strcmp', 'strncpy', 'tan']

# Add required imports
print 'import "regent"'
print 'local C = terralib.includecstring[['
print '#include <math.h>'
print '#include <stdlib.h>'
print '#include <stdio.h>'
print '#include <string.h>'
print ']]'
print 'local JSON = terralib.includec("json.h")'

for line in fileinput.input():
    line = line[:-1]
    # Remove some debug numbers
    line = re.sub(r'region#[0-9]+', r'region', line)
    line = re.sub(r'ispace#[0-9]+', r'ispace', line)
    # Remove type annotations where they can be inferred
    line = re.sub(r'for ([\w$#]+) : .* in ', r'for \1 in ', line)
    line = re.sub(r'var ([\w$#]+) : [^:=]* =', r'var \1 =', line)
    line = re.sub(r' : partition\(.*', '', line)
    # Make variable names valid identifiers
    # Regent vars: $abc, $abc#123, $123
    # Terra vars: abc, $abc, abc$123, $abc$123
    line = re.sub(r'\$(\w+)#([0-9]+)', r'\1__\2', line)
    line = re.sub(r'\$?(\w+)\$([0-9]+)', r'\1__\2', line)
    line = re.sub(r'\$([0-9]+)', r'__\1', line)
    line = re.sub(r'\$(\w+)', r'\1', line)
    # Wrap type casts
    line = re.sub(r'(float\[[0-9]+\])\(', r'[\1](', line)
    line = re.sub(r'(double\[[0-9]+\])\(', r'[\1](', line)
    line = re.sub(r'(&\w+)\(', r'[\1](', line)
    # inf -> math.huge
    line = re.sub(r'([^\w])inf', r'\1math.huge', line)
    # Add proper namespaces
    line = re.sub(r'([^\w.])(%s)\(' % '|'.join(STDLIB_FUNS), r'\1C.\2(', line)
    line = line.replace('_IO_FILE', 'C._IO_FILE')
    line = line.replace('H5', 'HDF5.H5')
    line = line.replace('legion_', 'regentlib.c.legion_')
    line = line.replace('std.', 'regentlib.')
    line = line.replace('base.', 'regentlib.')
    line = re.sub(r'(_?json_)', r'JSON.\1', line)
    line = line.replace('parseConfig(', 'SCHEMA.parseConfig(')
    # Remove unparsable terra annotations
    line = re.sub(r'extern global (.*) : \w+', r'\1', line)
    # Fix some terra mis-prints
    line = line.replace('if not', 'if not ')
    line = re.sub(r'\(_0,[ _0-9,]* = ([^)]*)\)', r'({\1})', line)
    # (@x). -> x.
    line = re.sub(r'\(@(\w+)\)\.', r'\1.', line)
    # Remove unnecessary casts
    line = re.sub(r'uint32\(([0-9]+)\)', r'\1', line)
    line = re.sub(r'int32\(([0-9]+)\)', r'\1', line)
    line = re.sub(r'double\(([0-9]+)\)', r'\1.0', line)
    # Print filtered line (unless it's a comment)
    if not line.startswith('--'):
        print line

# Add final compile command
LIBS = '"-ljsonparser","-lm"' + ((',"-l%s"' % HDF_LIBNAME) if USE_HDF else '')
print 'regentlib.saveobj(main, "%s", "executable", nil, {%s})' % (OBJNAME, LIBS)
