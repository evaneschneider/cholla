CHOLLA
============
A 3D GPU-based hydrodynamics code (Schneider & Robertson, ApJS, 2015).

Getting started
----------------
This is a branch of the *Cholla* code with dynamic 3D gravity.
It includes a variety of 3D Poisson solvers selected using the `POISSON_SOLVER` environment variable to define a compile-time macro.
The following settings of `POISSON_SOLVER` are currently supported.
- `-DPFFT` uses the *PFFT* library to compute 3D FFTs on host processors, distributed with MPI.
This option requires at least two MPI tasks.
- `-DCUFFT` uses the *CuFFT* library to compute 3D FFTs on a single GPU.
This options can use only one MPI task.
- `-DPARIS` uses the *Paris* Poisson solvers, provided in the `src/gravity` directory.
This solvers perform FFTs on GPUs, distributed with MPI.
The intent is to extend and tune the Paris solvers to run efficiently on exascale computers.
- `-DPFFT -DPARIS` or `-DCUFFT -DPARIS` uses the *PFFT* or *CuFFT* solver, respectively, and compares the result of each Poisson solve against the result of the *Paris* periodic solver.
At the beginning of the run, it also compares each solver against an analytic solution. The comparisons are single-line printouts of the L1, L2, and L-infinity norms.
- `-DSOR -DPARIS` uses an iterative successive-over-relaxation solver with nonzero nonperiodic boundary conditions, and it compares the result of each Poisson solve against the result of the *Paris* zero-valued-boundary-condition solver. The Cholla solver that calls Paris in this case subtracts an analtyic density before calling Paris and adds an analytic potential afterward. These analytic fields match the nonzero boundary conditions of the problem. At the beginning of the run, Cholla compares the SOR solver and the Paris solver against an analytic solution. Again, the comparisons are single-line printouts of the L1, L2, and L-infinity norms.

*Cholla* is designed to 
be run using GPUs, and can be run in serial mode using one GPU
or with MPI.

A Makefile is included in this directory, and sample build scripts are in `builds/scripts`. After downloading the code, you should
be able to configure it for your machine by creating a build script based on one of the following examples.

| Script | Description |
| --- | --- |
| `amake-*.sh` | Build for AMD GPUs using MPI with host-memory message buffers. |
| `aomake-*.sh` | Build for AMD GPUs using MPI with GPU-memory message buffers. |
| `nmake-*.sh` | Build for NVidia GPUs using MPI with GPU-memory message buffers. |
| `smake-*.sh` | Build for NVidia GPUs on OLCF Summit using MPI with GPU-memory message buffers. |

You can run the sample scripts directly from the root Cholla directory, and they will install a Cholla binary in the `bin` directory. For example,
```
builds/scripts/amake-hydro.sh
```
installs `bin/cholla.hydro-amd`.


Configuration Notes
------------
Most of the configuration options available in *Cholla* are selected by commenting/uncommenting
the appropriate line in the Makefile or by setting environment variables in the build script.
Examples of configurations that require edits to the Makefile include single vs
double precision, output format, the reconstruction method, Riemann solver, integrator, and cooling. The entire code must be recompiled any time you change the configuration. For more information on the various options, see the "[Makefile](https://github.com/cholla-hydro/cholla/wiki/Makefile)" page of the wiki.


Running Cholla
--------------
To run the code after it is compiled, you must supply an input file with parameters and a problem that matches a function
in the `initial_conditions` file. For example, to run a 1D Sod Shock Tube test, you would type

```./cholla tests/1D/Sod.txt```

in the directory with the `cholla` binary. Some output will be generated in the terminal, and output files will be written in the directory specified in the input parameter file.

To run *Cholla* in parallel mode, the `CHOLLA_MPI` macro must be defined at compile time, either in the Makefile or by the build script. Then you can run
using a command like the following.

```srun -n4 ./cholla tests/1D/Sod.txt```

Each process will be assigned a GPU. *Cholla* cannot be run with more processes than available GPUs,
so MPI mode is most useful on a cluster (or for testing parallel behavior with a single process).

The `tests/scripts` directory contains `*run-*.sh` examples for each of the `*make-*.sh` examples. They create output directories in the `run` directory. For example
```
test/scripts/arun-hydro.sh
```
runs and generates output in `run/out.hyrdo-amd.<date>.<process ID>`.

More information about compiling and running *Cholla* can be found in the wiki associated with this repository.

Other Notes
--------------

*Cholla* can be run without GPUs by commenting out `-DCUDA` in the Makefile, but this configuration is not recommended. Because *Cholla*
was designed with GPUs in mind, the CPU performance is lackluster at best. In addition, some 
of the configuration options are not available in the non-GPU mode (and warnings are not always included).

When running tests in fewer than 3 dimensions, *Cholla* assumes that the X direction will be used first, then
the Y, then Z. This is to say, in 1D `nx` must always be greater than 1, and in 2D `nx` and `ny` must be greater than 1.

In practice, we have found the Van Leer integrator to be the most stable. *Cholla* is set to run with a default CFL coefficient of 0.3, but this can be changed within the grid initialization function.

Currently, *Cholla* is not designed to be able to run very large (>1e6) 1 dimensional problems. If this is a functionality you are
interested in, please let us know.
