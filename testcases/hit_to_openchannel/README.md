This directory contains the scripts, inputs and processed output from the automated LF search study performed on Lassen in early 2019.

Contents
========

* `all_combinations.py`: enumerates all valid coarsening parameter combinations within a search space
* `average_runtime.py`: computes average running time of multiple samples of same level
* `chain_jobs.py`: breaks down a long-running simulation into multiple cluster jobs
* `hf-*.json`: HF configuration, partitioned according to the limits of different clusters
* `job_status.sh`: queries job status, post-processes output of completed runs
* `make_cases.sh`: instantiates a concrete case configuration for each choice of uncertain inputs
* `make_levels.sh`: produces an LF base case for each coarsening parameter combination
* `merge_files.sh`: merges results from the consituent jobs of a long-running simulation
* `time_average_all.sh`: averages quantities at the outlet across iterations
* `uncertainties_*.dat`: preselected samplings of uncertain inputs
* `results/`: processed output from Lassen run
  * `*.csv`: lists each LF's coarsening parameter values, correlation with HF and running time
  * `model.py`: script to evaluate classifier-based acceleration method
  * `search.py`: script to evaluate search-based acceleration methods

Instructions (Lassen)
=====================

Create all levels:

```
mkdir hf
cp hf-sierra.json hf/base.json
./all_combinations.py 128,64 32,16 32,16 10,100 64,32 16,8 16,8 86,50 true,false 4,3 3,2 > levels.dat
./make_levels.sh levels.dat hf-sierra.json # creates directories lf{0..1087}/
```

Perform HF runs (must be split into multiple jobs):

```
cd hf
mkdir -p "$SCRATCH"/hf/aggr
../make_cases.sh ../uncertainties_32.dat base.json # creates files hf/case{0..31}.json
../chain_jobs.py "$SCRATCH"/hf case{0..31}.json 150000
# follow the printed instructions to launch all jobs (assuming 16 jobs in this example)
# wait for runs to finish
for I in {0..31}; do ../merge_files.sh "$SCRATCH"/hf/job{0..15}/sample$((2*I+1))/probe0.csv > "$SCRATCH"/hf/aggr/sample$((2*I+1))/probe0.csv; done
../time_average_all.sh 32 . "$SCRATCH"/hf/aggr > averages.csv
../average_runtime.py -n 64 "$SCRATCH"/hf/job*/sample*/console.txt
```

Perform LF runs (short enough to fit in one job):

```
cd lf"$I"
../make_cases.sh ../uncertainties_32.dat base.json # creates files "lf$I"/case{0..31}.json
RANKS_PER_NODE=4 "$SOLEIL_DIR"/src/soleil.sh $(echo -m\ case{0..31}.json)
# wait for runs to finish
RUNTIME=1 AVERAGE=1 ../job_status.sh . # time averages printed to csv file
```

Instructions (other machines)
=============================

The same scripts should work on other machines, but you need to edit (at least) the following:

* set `RANKS_PER_NODE` to the number of GPUs per node
* supply an appropriate value for `iters_per_job` to `chain_jobs.py`
* set `wallTime` in the base HF configuration according to the cluster's job time limit
* set `tiles` and `collocateSections` in the base HF configuration according to the GPU memory limits
