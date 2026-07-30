#!/bin/bash

#-- This script needs to be sourced in the terminal, e.g.
#   source ./setup.crc.gcc.sh

module load cuda/11.8  gcc/10.2.0  openmpi/4.0.5  hdf5/1.12.1

echo "mpicxx --version is: "
mpicxx --version

# export MPI_GPU="-DMPI_GPU"
export F_OFFLOAD="-fopenmp"
export CHOLLA_ENVSET=1
