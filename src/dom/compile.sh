#!/bin/bash -eu

DOM_DIR="$(cd "$(dirname "$(perl -MCwd -le 'print Cwd::abs_path(shift)' "${BASH_SOURCE[0]}")")" && pwd)"

"$LEGION_DIR"/language/regent.py "$DOM_DIR"/host.rg
