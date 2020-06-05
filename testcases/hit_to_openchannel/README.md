This directory contains the scripts, input and output from the automated low-fidelity model selection case study performed on Lassen in early 2019

Experiment description
======================

This case study was run to evaluate a prototype implementation for a search-based approach that, given a high-fidelity (HF) simulation, finds an optimal coarsening of the HF's physics to use as a low-fidelity (LF) model in a Control Variates-based UQ study.

The method being evaluated operates as follows:

* Build a case configuration for the HF (we used the case described in [1]).
* Designate a single metric as the quantity of interest (we used the same one as the original paper, space-averaged over the outlet plane and time-averaged starting from 1 flow-through time until the end of the simulation).
* Designate some of the simulation's inputs as uncertain (we chose 6 of the simulation's physical parameters).
* Designate a distribution for the uncertain inputs (we assumed each uncertain input is uniformly and independently distributed within a ±5% range around a nominal value).
* Sample a number of inputs from this distribution (32 samples in our case).
* Evaluate the HF on all inputs.
* Define a search space of LF models (we exposed 11 simulation parameters as coarsening knobs, with 2 choices for each, for a total of 1088 LF combinations).
* Evaluate each LF on all inputs.
* Compare LFs in terms of their running time and their correlation with the HF (sample correlation between the HF and LF outputs on the same set of inputs).

Human experts were given the same problem and their answers were compared against the optimal configurations discovered by the search. Different acceleration ideas were evaluated using the data from the full run.

More details about this method and case study can be found in the accompanying paper (currently under submission).

Directory contents
==================

* `all_combinations.py`: enumerates all valid coarsening parameter combinations within a search space
* `average_runtime.py`: computes average running time of multiple samples of the same level
* `chain_jobs.py`: breaks down a long-running simulation into multiple cluster jobs
* `hf-*.json`: HF configuration, partitioned according to the limits of different clusters
* `job_status.sh`: queries job status, post-processes output of completed runs
* `make_cases.sh`: instantiates a concrete case configuration for each choice of uncertain inputs
* `make_levels.sh`: produces an LF base case for each coarsening parameter combination
* `merge_files.sh`: merges results from the constituent jobs of a long-running simulation
* `time_average_all.sh`: averages quantities at the outlet across iterations
* `uncertainties_*.dat`: preselected samplings of uncertain inputs
* `results/`: processed output from Lassen run
  * `*.csv`: lists each LF's coarsening parameter values, correlation with HF and running time
  * `model.py`: script to evaluate classifier-based acceleration method
  * `search.py`: script to evaluate search-based acceleration methods

Experiment configuration
========================

This study was performed on the Lassen supercomputer. At the time of the study, Lassen's hardware configuration was as follows:

* CPUs: 44 IBM Power9 3.5GHz cores per node
* GPUs: 4 NVIDIA Tesla V100 GPUs per node
* CPU memory:: 256GB per node
* GPU memory: 64GB per node
* Interconnect: EDR InfiniBand
* Operating system: RHEL

Each sample was ran on a dedicated GPU.

For this study, Soleil-X was built and run using the following software:

* GCC 7.3.1
* LLVM 3.8
* CUDA 9.2.148
* GASNet 1.30.0
* HDF5 1.10.1
* Legion: nightly version from early 2019
* Terra: nightly version from early 2019

Instructions to reproduce (Lassen)
==================================

Create directories for HF and LF runs:

```
mkdir hf
cp hf-sierra.json hf/base.json
./all_combinations.py 128,64 32,16 32,16 10,100 64,32 16,8 16,8 86,50 true,false 4,3 3,2 > levels.dat
./make_levels.sh levels.dat hf-sierra.json # creates directories lf{0..1087}/
```

Perform HF runs (must be split into multiple jobs):

```
cd hf
for I in {0..31}; do mkdir -p "$SCRATCH"/hf/aggr/sample$((2*I+1)); done
../make_cases.sh ../uncertainties_32.dat base.json # creates files hf/case{0..31}.json
../chain_jobs.py "$SCRATCH"/hf case{0..31}.json 150000
# follow the printed instructions to launch all jobs (16 jobs in this case)
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

The different LFs can be compared in terms of average running time across the 32 samples (printed by the `job_status.sh` script), and in terms of correlation with the HF (computed by running the sample correlation formula on the csv files produced above, using any statistical package).

Different search acceleration ideas can be evaluated from the results of the full run:

```
cd results/
./model.py
./search.py
```

Instructions to reproduce (other machines)
==========================================

The same scripts should work on other machines, but you will need to edit (at least) the following:

* set `RANKS_PER_NODE` to the number of GPUs per node
* supply an appropriate value for `iters_per_job` to `chain_jobs.py`
* set `wallTime` in the base HF configuration according to the cluster's job time limit
* set `tiles` and `collocateSections` in the base HF configuration according to the GPU memory limits

References
==========

[1] M. Rahmani, G. Geraci, G. Iaccarino, and A. Mani, "Effects of Particle Polydispersity on Radiative Heat Transfer in Particle-Laden Turbulent Flows", International Journal of Multiphase Flow, vol. 104, pp. 42–59, 2018
