#!/usr/bin/env python2

import argparse
import csv
import numpy as np
import glob
import os
import sys

parser = argparse.ArgumentParser()
parser.add_argument('-n', '--num-cases', type=int, default=32)
parser.add_argument('launch_dir', nargs='+')
args = parser.parse_args()

columns = []
for d in args.launch_dir:
    csv_files = glob.glob(os.path.join(d, '*.csv'))
    if len(csv_files) == 0:
        columns.append([''] * args.num_cases)
        continue
    latest_csv = max(csv_files)
    with open(latest_csv) as f:
        reader = csv.reader(f, dialect=csv.excel_tab)
        next(reader)
        column = [row[0] for row in reader]
        assert len(column) == args.num_cases
        columns.append(column)

arr = np.transpose(np.array(columns))
writer = csv.writer(sys.stdout, dialect=csv.excel_tab)
writer.writerows(arr.tolist())
