#!/usr/bin/env python2

import argparse
import csv
import glob
import os
import sys

parser = argparse.ArgumentParser()
parser.add_argument('-n', '--num-cases', type=int, default=32)
parser.add_argument('-c', '--column', type=int, default=0)
parser.add_argument('launch_dir', nargs='+')
args = parser.parse_args()

rows = []
for d in args.launch_dir:
    csv_files = glob.glob(os.path.join(d, '*.csv'))
    if len(csv_files) == 0:
        rows.append([''] * args.num_cases)
        continue
    latest_csv = max(csv_files)
    with open(latest_csv) as f:
        reader = csv.reader(f, dialect=csv.excel_tab)
        next(reader)
        column = [row[args.column] for row in reader]
        assert len(column) == args.num_cases
        rows.append(column)

writer = csv.writer(sys.stdout, dialect=csv.excel_tab)
writer.writerows(rows)
