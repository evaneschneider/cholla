/*! \file hydro_cuda.h
 *  \brief Declarations of functions used in all cuda integration algorithms. */

#ifdef CUDA
#ifndef HYDRO_CUDA_H
#define HYDRO_CUDA_H

#include"global.h"


__global__ void Update_Conserved_Variables_2D(Real *dev_conserved, Real *dev_F_x, Real *dev_F_y, int nx, int ny,
                                              int n_ghost, Real dx, Real dy, Real dt, Real gamma);

__global__ void Calc_dt_2D(Real *dev_conserved, int nx, int ny, int n_ghost, Real dx, Real dy, Real *dti_array, Real gamma);

__global__ void Sync_Energies_2D(Real *dev_conserved, int nx, int ny, int n_ghost, Real gamma);


#endif //HYDRO_CUDA_H
#endif //CUDA

