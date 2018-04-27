#!/usr/bin/env python

# Written by Elliott Slaughter
# See https://upc-bugs.lbl.gov/bugzilla/show_bug.cgi?id=3447

from __future__ import print_function

import os
import re
import subprocess
import sys

# Retreive the list of HCAs on the current node.
devinfo = subprocess.check_output(['ibv_devinfo'])
hcas = re.findall(r'hca_id:\s+([A-Za-z0-9_]+)', devinfo)
assert(len(hcas) > 0)

# Pick an HCA by round-robin over ranks.
rank = int(os.environ['JSM_NAMESPACE_LOCAL_RANK'])
hca = hcas[rank % len(hcas)]

env = dict(list(os.environ.items()) + [('GASNET_IBV_PORTS', hca)])
assert(len(sys.argv) >= 2)
sys.exit(subprocess.call(sys.argv[1:], env=env))
