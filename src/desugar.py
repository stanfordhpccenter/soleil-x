#!/usr/bin/env python

import fileinput

for line in fileinput.input():
    line = line[:-1]
    line = line.replace('@ESCAPE', '[(function() local __quotes = terralib.newlist()')
    line = line.replace('@EPACSE', 'return __quotes end)()];')
    line = line.replace('@EMIT', '__quotes:insert(rquote')
    line = line.replace('@TIME', 'end)')
    print line
