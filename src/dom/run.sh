#!/bin/bash -eu

DOM_DIR="$(cd "$(dirname "$(perl -MCwd -le 'print Cwd::abs_path(shift)' "${BASH_SOURCE[0]}")")" && pwd)"

LD_LIBRARY_PATH=.:"$LEGION_DIR"/bindings/regent/ "$DOM_DIR"/a.out
