#!/bin/bash -eu

"$SOLEIL_DIR"/scripts/compare_console.py sample0/console.txt ref/sample0/console.txt
"$SOLEIL_DIR"/scripts/compare_console.py sample1/console.txt ref/sample1/console.txt
