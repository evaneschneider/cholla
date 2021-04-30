#!/bin/bash

module load xl cuda hdf5

export GPU_MPI="-DGPU_MPI"

export CHOLLA_ENVSET=1

make clean
make TYPE=gpu_hydro -j 8
