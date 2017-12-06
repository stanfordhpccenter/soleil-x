#!/usr/bin/env python

# NOTE: Run Admiral with DEBUG=1, set '-fdebug 1', and feed the stdout
# (not stderr) output into this script.

import fileinput
import os
import re

HDF_LIBNAME = os.environ['HDF_LIBNAME'] if 'HDF_LIBNAME' in os.environ else 'hdf5'
HDF_HEADER = os.environ['HDF_HEADER'] if 'HDF_HEADER' in os.environ else 'hdf5.h'

# Add required imports
print 'import "regent"'
print 'local HDF5 = terralib.includec("%s")' % HDF_HEADER
print 'local JSON = terralib.includec("json.h")'
print 'local C = terralib.includecstring[['
print '#include <math.h>'
print '#include <stdlib.h>'
print '#include <stdio.h>'
print '#include <string.h>'
print ']]'

for line in fileinput.input():
    line = line[:-1]
    # Remove type annotations where they can be inferred
    line = re.sub(r'for ([\w$#]+) : .* in ', r'for \1 in ', line)
    line = re.sub(r'var ([\w$#]+) : [^:=]* =', r'var \1 =', line)
    # Make variable names valid identifiers
    line = re.sub(r'\$\w+#([0-9]+)', r'v\1', line)
    line = re.sub(r'\$(\w+)\$([0-9]+)', r'v\1__\2', line)
    line = re.sub(r'\$(\w+)', r'v\1', line)
    # Remove remaining debug numbers
    line = re.sub(r'region#[0-9]+', r'region', line)
    line = re.sub(r'ispace#[0-9]+', r'ispace', line)
    # Wrap type casts
    line = re.sub(r'(float\[[0-9]+\])\(', r'[\1](', line)
    line = re.sub(r'(double\[[0-9]+\])\(', r'[\1](', line)
    line = re.sub(r'(&\w+)\(', r'[\1](', line)
    # inf -> math.huge
    line = re.sub(r'([^\w])inf', r'\1math.huge', line)
    # Add proper namespaces
    line = re.sub(r'([^\w])(printf|snprintf|strcmp|malloc|free|memcpy|exit)\(', r'\1C.\2(', line)
    line = re.sub(r'([^\w])(acos|asin|atan|cbrt|tan|pow|fmod|ceil)\(', r'\1C.\2(', line)
    line = re.sub(r'([^\w])(fopen|fseek|ftell|fread|fclose|rand)\(', r'\1C.\2(', line)
    line = line.replace('_IO_FILE', 'C._IO_FILE')
    line = line.replace('H5', 'HDF5.H5')
    line = line.replace('legion_', 'regentlib.c.legion_')
    line = line.replace('std.', 'regentlib.')
    line = re.sub(r'(_?json_)', r'JSON.\1', line)
    # Remove unparsable terra annotations
    line = re.sub(r'extern global (.*) : \w+', r'\1', line)
    # Fix some terra mis-prints
    line = line.replace('if not', 'if not ')
    line = re.sub(r'\(_0,[ _0-9,]* = ([^)]*)\)', r'({\1})', line)
    # (@x). -> x.
    line = re.sub(r'\(@(\w+)\)\.', r'\1.', line)
    # Print filtered line
    print line

# Add final compile command
print 'regentlib.saveobj(main, "a.out", "executable", nil, {"-ljsonparser","-lm","-l%s"})' % HDF_LIBNAME
