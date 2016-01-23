-----------------------------------------------------------
  Soleil-X version 0.0.1
-----------------------------------------------------------

Soleil-X is a turbulence/particle/radiation solver written in the Liszt-Ebb DSL for execution with the Legion runtime.

---------------------------------------------------
  Soleil-X Installation
---------------------------------------------------

Soleil-X is an application written in the Ebb domain specific language (part of the Liszt project at Stanford University), which means that there is no explicit compilation required by the user (Liszt takes care of this using a just-in-time compilation process). Think of the Soleil-X application as a script you would prepare for a dynamic programming language like Python. You can find the source code for Soleil-X in the file 'soleil-x.t' under the src/ directory.

To execute Soleil-X application, at a minimum, you must download and install the Liszt-Ebb DSL. Ebb is open-source and available on GitHub [here](https://github.com/gilbo/liszt-ebb). See the instructions in its README for the basic install. More information on the Ebb language can be found [here](http://ebblang.org).

Research to connect the Liszt-Ebb DSL with the [Legion Programming System](https://github.com/StanfordLegion/legion) is underway. Legion is a runtime for achieving high performance on distributed, heterogeneous hardware. More information on Legion can be found [here](http://legion.stanford.edu).

----------------------------------------------------------
  Executing Soleil-X
----------------------------------------------------------

After installing Liszt-Ebb and obtaining a copy of the Soleil-X source repository, one can execute the solver. Apart from the source script, you will need to provide a configuration file as input to the solver that describes your particular problem conditions. A thoroughly commented template configuration file is available to serve as an example (params_template.lua) in the root directory of the source. In addition, there are several existing configuration files available in the testcases/ directory.

Get started quickly by doing the following:

1. Prepare your own config file or use one of the examples. We will use the Taylor Green Vortex case that is provided.

2. Move to the directory containing the config file.

    ```
    $ cd soleil-x/testcases/taylor_green_vortex
    ```

3. Execute the script through Ebb and provide the config file as a command line arg (with the '-f' option). You will see console output as the case runs, and any output will be written to the current working directory, such as restart files or files to visualize the flow and particle solutions with Tecplot.

    ```
    $ /path/to/ebb ../../src/soleil-x.t -f taylor_green_vortex.lua
    ```

That's it!

----------------------------------------------------------
  Partitioning Options (UNSTABLE - DEVELOPERS ONLY)
----------------------------------------------------------

There are new options to control the partitioning of the underlying grid and particle relations in Liszt. They take the form of optional command line args after the config file specification. The grid partitions are defined in each direction independently (x,y,z) along with a separate value for partitioning particles:

```
$ /path/to/ebb /path/to/soleil-x.t -f your_config.lua -x 2 -y 4 -z 5 -p 1
```

where here we are requesting 2 partitions along the x direction, 4 in the y, and 5 in the z for the grid (fluid), and no partitioning (1) for the particles.

----------------------------------------------------------
  Soleil-X Developers
----------------------------------------------------------

Soleil-X is under active development as part of the Predictive Science Academic Alliance Program II (PSAAP II) supported by the Department of Energy National Nuclear Security Administration under Award Number DE-NA0002373-1.

The initial Soleil-X source code was developed by:

   - Dr. Thomas D. Economon (economon@stanford.edu)
   - Dr. Ivan Bermejo-Moreno

with much ongoing support and help from:

   - Prof. Pat Hanrahan's group at Stanford University (Liszt).
   - Prof. Alex Aiken's group at Stanford University (Legion).
