#!/bin/bash -eu

###############################################################################
# Inputs
###############################################################################

# Which queue/partition to use (if a scheduler is available)
export QUEUE="${QUEUE:-}"

# Which job to wait for before starting (if a scheduler is available)
export AFTER="${AFTER:-}"

# Whether to use GPUs (if available)
export USE_CUDA="${USE_CUDA:-1}"

# Whether to emit Legion profiler logs
export PROFILE="${PROFILE:-0}"

# Whether to print a backtrace on crash (interferes with signal handling)
export GASNET_BACKTRACE="${GASNET_BACKTRACE:-1}"

# Whether to freeze Legion on crash
export LEGION_FREEZE_ON_ERROR="${LEGION_FREEZE_ON_ERROR:-0}"

###############################################################################
# Helper functions
###############################################################################

function quit {
    echo "$1" >&2
    exit 1
}

function read_json {
    python2 -c "
import json
def wallTime(sample):
  return int(sample['Mapping']['wallTime'])
def numRanks(sample):
  tiles = sample['Mapping']['tiles']
  tilesPerRank = sample['Mapping']['tilesPerRank']
  xRanks = int(tiles[0]) / int(tilesPerRank[0])
  yRanks = int(tiles[1]) / int(tilesPerRank[1])
  zRanks = int(tiles[2]) / int(tilesPerRank[2])
  return xRanks * yRanks * zRanks
f = json.load(open('$1'))
if '$2' == 'single':
  print wallTime(f), numRanks(f)
elif '$2' == 'dual':
  print max(wallTime(f['configs'][0]), wallTime(f['configs'][1])), \
        numRanks(f['configs'][0]) + numRanks(f['configs'][1])
else:
  assert(false)"
}

###############################################################################
# Derived options
###############################################################################

export EXECUTABLE="$SOLEIL_DIR"/src/soleil.exec

# Total wall-clock time is the maximum across all samples.
# Total number of ranks is the sum of all sample rank requirements.
MINUTES=0
NUM_RANKS=0
function parse_config {
    read -r _MINUTES _NUM_RANKS <<<"$(read_json "$@")"
    MINUTES=$(( MINUTES > _MINUTES ? MINUTES : _MINUTES ))
    NUM_RANKS=$(( NUM_RANKS + _NUM_RANKS ))
}
for (( i = 1; i <= $#; i++ )); do
    j=$((i+1))
    if [[ "${!i}" == "-i" ]] && (( $i < $# )); then
        parse_config "${!j}" "single"
    elif [[ "${!i}" == "-m" ]] && (( $i < $# )); then
        parse_config "${!j}" "dual"
    fi
done
if (( NUM_RANKS < 1 )); then
    quit "No configuration files provided"
fi
export MINUTES
export NUM_RANKS

###############################################################################

source "$SOLEIL_DIR"/src/run.sh
