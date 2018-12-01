/*! \file grid3D.cpp
 *  \brief Definitions of the Grid3D class */

#include <stdlib.h>
#include <math.h>
#include <string.h>
#ifdef HDF5
#include <hdf5.h>
#endif
#include "global.h"
#include "grid3D.h"
#include "CTU_1D.h"
#include "CTU_2D.h"
#include "CTU_3D.h"
#include "CTU_1D_cuda.h"
#include "CTU_2D_cuda.h"
#include "CTU_3D_cuda.h"
#include "VL_1D_cuda.h"
#include "VL_2D_cuda.h"
#include "VL_3D_cuda.h"
#include "io.h"
#include "error_handling.h"
#ifdef MPI_CHOLLA
#include <mpi.h>
#ifdef HDF5
#include <H5FDmpio.h>
#endif
#include "mpi_routines.h"
#endif
#include <stdio.h>
#include "flux_correction.h"
#ifdef CLOUDY_COOL
#include "cooling_wrapper.h"
#endif


/*! \fn Grid3D(void)
 *  \brief Constructor for the Grid. */
Grid3D::Grid3D(void)
{
  // set initialization flag to 0
  flag_init = 0;

  // set number of ghost cells
  #ifdef PCM
  H.n_ghost = 2;
  #endif //PCM
  #ifdef PLMP
  H.n_ghost = 3;
  #endif //PLMP
  #ifdef PLMC
  H.n_ghost = 3;
  #endif //PLMC
  #ifdef PPMP
  H.n_ghost = 4;
  #endif //PPMP
  #ifdef PPMC
  H.n_ghost=4;
  #endif //PPMC

}

/*! \fn void Get_Position(long i, long j, long k, Real *xpos, Real *ypos, Real *zpos)
 *  \brief Get the cell-centered position based on cell index */ 
void Grid3D::Get_Position(long i, long j, long k, Real *x_pos, Real *y_pos, Real *z_pos)
{

#ifndef   MPI_CHOLLA

  *x_pos = H.xbound + H.dx*(i-H.n_ghost) + 0.5*H.dx;
  *y_pos = H.ybound + H.dy*(j-H.n_ghost) + 0.5*H.dy;
  *z_pos = H.zbound + H.dz*(k-H.n_ghost) + 0.5*H.dz;

#else   /*MPI_CHOLLA*/

  /* position relative to local xyz bounds */
  *x_pos = H.xblocal + H.dx*(i-H.n_ghost) + 0.5*H.dx;
  *y_pos = H.yblocal + H.dy*(j-H.n_ghost) + 0.5*H.dy;
  *z_pos = H.zblocal + H.dz*(k-H.n_ghost) + 0.5*H.dz;
  

#endif  /*MPI_CHOLLA*/

}


/*! \fn void Initialize(int nx_in, int ny_in, int nz_in)
 *  \brief Initialize the grid. */
void Grid3D::Initialize(struct parameters *P)
{
  // number of fields to track (default 5 is # of conserved variables)
  H.n_fields = 5;

  // if including passive scalars increase the number of fields
  #ifdef SCALAR
  H.n_fields += NSCALARS;
  #endif

  // if using dual energy formalism must track internal energy - always the last field!
  #ifdef DE
  H.n_fields++;
  #endif  

  int nx_in = P->nx;
  int ny_in = P->ny;
  int nz_in = P->nz;

  // Set the CFL coefficient (a global variable)
  C_cfl = 0.2;

  // Set the output timestep
  H.out_step = P->gridstep;
  #ifdef PROJECTION
  H.out_step = fmin(P->projstep, H.out_step);
  #endif
  #ifdef SLICES
  H.out_step = fmin(P->slicestep, H.out_step);
  #endif

#ifndef MPI_CHOLLA

  // set grid dimensions
  H.nx = nx_in+2*H.n_ghost;
  H.nx_real = nx_in;
  if (ny_in == 1) H.ny = 1;
  else H.ny = ny_in+2*H.n_ghost;
  H.ny_real = ny_in;
  if (nz_in == 1) H.nz = 1;
  else H.nz = nz_in+2*H.n_ghost;
  H.nz_real = nz_in;

  // set total number of cells
  H.n_cells = H.nx * H.ny * H.nz;

#else  /*MPI_CHOLLA*/

  /* perform domain decomposition
   * and set grid dimensions      
   * and allocate comm buffers */ 
  DomainDecomposition(P, &H, nx_in, ny_in, nz_in); 
  
#endif /*MPI_CHOLLA*/

  // failsafe
  if(H.n_cells<=0)
  {
    chprintf("Error initializing grid: H.n_cells = %d\n", H.n_cells);
    chexit(-1);
  }

  // check for initilization
  if(flag_init)
  {
    chprintf("Already initialized. Please reset.\n");
    return;
  }
  else
  {
    // mark that we are initializing
    flag_init = 1;
  }

  // Set the flag that tells Update_Grid which buffer to read from
  gflag = 0;

  // Set header variables for time within the simulation
  H.t = 0.0;
  // and the number of timesteps taken
  H.n_step = 0;
  // and the wall time
  H.t_wall = 0.0;
  // and inialize the timestep
  H.dt = 0.0;


  // allocate memory
  AllocateMemory();


#ifdef ROTATED_PROJECTION
  //x-dir pixels in projection 
  R.nx = P->nxr;
  //z-dir pixels in projection
  R.nz = P->nzr;
  //minimum x location to project
  R.nx_min = 0;
  //minimum z location to project
  R.nz_min = 0;
  //maximum x location to project
  R.nx_max = R.nx;
  //maximum z location to project
  R.nz_max = R.nz;
  //rotation angle about z direction
  R.delta = M_PI*(P->delta/180.); //convert to radians
  //rotation angle about x direction
  R.theta = M_PI*(P->theta/180.); //convert to radians
  //rotation angle about y direction
  R.phi = M_PI*(P->phi/180.); //convert to radians
  //x-dir physical size of projection
  R.Lx = P->Lx;
  //z-dir physical size of projection
  R.Lz = P->Lz;
  //initialize a counter for rotated outputs
  R.i_delta = 0;
  //number of rotated outputs in a complete revolution
  R.n_delta = P->n_delta;
  //rate of rotation between outputs, for an actual simulation
  R.ddelta_dt = P->ddelta_dt;
  //are we not rotating about z(0)?
  //are we outputting multiple rotations(1)? or rotating during a simulation(2)?
  R.flag_delta = P->flag_delta;
#endif /*ROTATED_PROJECTION*/

}


/*! \fn void AllocateMemory(void)
 *  \brief Allocate memory for the arrays. */
void Grid3D::AllocateMemory(void)
{


  // allocate memory for the conserved variable arrays
  // allocate all the memory to density, to insure contiguous memory
  buffer0 = (Real *) malloc(H.n_fields*H.n_cells*sizeof(Real));
  buffer1 = (Real *) malloc(H.n_fields*H.n_cells*sizeof(Real));

  // point conserved variables to the appropriate locations in buffer
  C.density  = &(buffer0[0]);
  C.momentum_x = &(buffer0[H.n_cells]);
  C.momentum_y = &(buffer0[2*H.n_cells]);
  C.momentum_z = &(buffer0[3*H.n_cells]);
  C.Energy   = &(buffer0[4*H.n_cells]);
  #ifdef SCALAR
  C.scalar  = &(buffer0[5*H.n_cells]);
  #endif
  #ifdef DE
  C.GasEnergy = &(buffer0[(H.n_fields-1)*H.n_cells]);
  #endif

  // initialize array
  for (int i=0; i<H.n_fields*H.n_cells; i++)
  {
    buffer0[i] = 0.0;
    buffer1[i] = 0.0;
  }

  #ifdef CLOUDY_COOL
  //printf("Warning: Cloudy cooling isn't currently working. No cooling will be applied.\n");
  Load_Cuda_Textures();
  #endif

  Set_Cluster_Locations();

}


/*! \fn void set_dt(Real dti)
 *  \brief Set the timestep. */
 void Grid3D::set_dt(Real dti)
{
  Real max_dti;

  if (H.n_step == 0) {
    max_dti = calc_dti_CPU();
  }
  else {
    #ifndef CUDA
    max_dti = calc_dti_CPU();
    #endif /*NO_CUDA*/
    #ifdef CUDA
    max_dti = dti;
    #endif /*CUDA*/
  }

  #ifdef MPI_CHOLLA
  max_dti = ReduceRealMax(max_dti);
  #endif /*MPI_CHOLLA*/
  
  /*
  if (H.n_step > 1) {
    H.dt = fmin(2*H.dt, C_cfl / max_dti);
  }
  else 
    H.dt = C_cfl / max_dti;
  */
  //chprintf("Within set_dt: %f %f %f\n", C_cfl, H.dt, max_dti);
  H.dt = C_cfl / max_dti;

}


/*! \fn Real calc_dti_CPU()
 *  \brief Calculate the maximum inverse timestep, according to the CFL condition (Toro 6.17). */
Real Grid3D::calc_dti_CPU()
{
  int i, j, k, id;
  Real d_inv, vx, vy, vz, P, cs;
  Real max_vx, max_vy, max_vz;
  Real max_dti = 0.0;
  max_vx = max_vy = max_vz = 0.0;

  // 1D
  if (H.nx > 1 && H.ny == 1 && H.nz == 1) {
    //Find the maximum wave speed in the grid
    for (i=H.n_ghost; i<H.nx-H.n_ghost; i++) {
      id = i;
      d_inv = 1.0 / C.density[id];
      vx = d_inv * C.momentum_x[id];
      vy = d_inv * C.momentum_y[id];
      vz = d_inv * C.momentum_z[id];
      P = fmax((C.Energy[id] - 0.5*C.density[id]*(vx*vx + vy*vy + vz*vz) )*(gama-1.0), TINY_NUMBER);
      cs = sqrt(d_inv * gama * P);
      // compute maximum cfl velocity
      max_vx = fmax(max_vx, fabs(vx) + cs);
    }
    // compute max inverse of dt
    max_dti = max_vx / H.dx;
  }
  // 2D
  else if (H.nx > 1 && H.ny > 1 && H.nz == 1) {
    // Find the maximum wave speed in the grid
    for (i=H.n_ghost; i<H.nx-H.n_ghost; i++) {
      for (j=H.n_ghost; j<H.ny-H.n_ghost; j++) {
        id = i + j*H.nx;
        d_inv = 1.0 / C.density[id];
        vx = d_inv * C.momentum_x[id];
        vy = d_inv * C.momentum_y[id];
        vz = d_inv * C.momentum_z[id];
        P = fmax((C.Energy[id] - 0.5*C.density[id]*(vx*vx + vy*vy + vz*vz) )*(gama-1.0), TINY_NUMBER);
        cs = sqrt(d_inv * gama * P);
        // compute maximum cfl velocity
        max_vx = fmax(max_vx, fabs(vx) + cs);
        max_vy = fmax(max_vy, fabs(vy) + cs);
      }
    }
    // compute max inverse of dt
    max_dti = max_vx / H.dx;
    max_dti = fmax(max_dti, max_vy / H.dy);
  }
  // 3D
  else if (H.nx > 1 && H.ny > 1 && H.nz > 1) {
    // Find the maximum wave speed in the grid
    for (i=0; i<H.nx-H.n_ghost; i++) {
      for (j=0; j<H.ny-H.n_ghost; j++) {
        for (k=0; k<H.nz-H.n_ghost; k++) {
          id = i + j*H.nx + k*H.nx*H.ny;
          d_inv = 1.0 / C.density[id];
          vx = d_inv * C.momentum_x[id];
          vy = d_inv * C.momentum_y[id];
          vz = d_inv * C.momentum_z[id];
          P = fmax((C.Energy[id] - 0.5*C.density[id]*(vx*vx + vy*vy + vz*vz) )*(gama-1.0), TINY_NUMBER);
          cs = sqrt(d_inv * gama * P);
          // compute maximum cfl velocity
          max_vx = fmax(max_vx, fabs(vx) + cs);
          max_vy = fmax(max_vy, fabs(vy) + cs);
          max_vz = fmax(max_vz, fabs(vz) + cs);
        }
      }
    }
    // compute max inverse of dt
    max_dti = max_vx / H.dx;
    max_dti = fmax(max_dti, max_vy / H.dy);
    max_dti = fmax(max_dti, max_vz / H.dy);
  } 
  else {
    chprintf("Invalid grid dimensions. Failed to compute dt.\n");
    chexit(-1);
  }

  return max_dti;

}



/*! \fn void Update_Grid(void)
 *  \brief Update the conserved quantities in each cell. */
Real Grid3D::Update_Grid(void)
{
  Real *g0, *g1;
  if (gflag == 0) {
    g0 = &(buffer0[0]);
    g1 = &(buffer1[0]);
  }
  else {
    g0 = &(buffer1[0]);
    g1 = &(buffer0[0]);
  }

  Real max_dti = 0;
  int x_off, y_off, z_off;

  // set x, y, & z offsets of local CPU volume to pass to GPU
  // so global position on the grid is known
  x_off = y_off = z_off = 0;
  #ifdef MPI_CHOLLA
  x_off = nx_local_start;
  y_off = ny_local_start;
  z_off = nz_local_start;
  #endif

  // Pass the structure of conserved variables to the CTU update functions
  // The function returns the updated variables
  if (H.nx > 1 && H.ny == 1 && H.nz == 1) //1D
  {
    #ifndef CUDA
    #ifndef VL
    CTU_Algorithm_1D(&(C.density[0]), H.nx, H.n_ghost, H.dx, H.dt);
    #endif //not_VL
    #ifdef VL
    chprintf("VL algorithm not implemented in non-cuda version.");
    chexit(-1);
    #endif //VL
    #endif //not_CUDA

    #ifdef CUDA
    #ifndef VL
    max_dti = CTU_Algorithm_1D_CUDA(g0, g1, H.nx, x_off, H.n_ghost, H.dx, H.xbound, H.dt, H.n_fields);
    #endif //not_VL
    #ifdef VL
    max_dti = VL_Algorithm_1D_CUDA(g0, g1, H.nx, x_off, H.n_ghost, H.dx, H.xbound, H.dt, H.n_fields);
    #endif //VL
    #endif //CUDA
  }
  else if (H.nx > 1 && H.ny > 1 && H.nz == 1) //2D
  {
    #ifndef CUDA
    #ifndef VL
    CTU_Algorithm_2D(&(C.density[0]), H.nx, H.ny, H.n_ghost, H.dx, H.dy, H.dt);
    #endif //not_VL
    #ifdef VL
    chprintf("VL algorithm not implemented in non-cuda version.");
    chexit(-1);    
    #endif //VL
    #endif //not_CUDA

    #ifdef CUDA
    #ifndef VL
    max_dti = CTU_Algorithm_2D_CUDA(g0, g1, H.nx, H.ny, x_off, y_off, H.n_ghost, H.dx, H.dy, H.xbound, H.ybound, H.dt, H.n_fields);
    #endif //not_VL
    #ifdef VL
    max_dti = VL_Algorithm_2D_CUDA(g0, g1, H.nx, H.ny, x_off, y_off, H.n_ghost, H.dx, H.dy, H.xbound, H.ybound, H.dt, H.n_fields);
    #endif //VL
    #endif //CUDA
  }
  else if (H.nx > 1 && H.ny > 1 && H.nz > 1) //3D
  {
    #ifndef CUDA
    #ifndef VL
    CTU_Algorithm_3D(&(C.density[0]), H.nx, H.ny, H.nz, H.n_ghost, H.dx, H.dy, H.dz, H.dt);
    #endif //not_VL
    #ifdef VL
    chprintf("VL algorithm not implemented in non-cuda version.");
    chexit(-1);    
    #endif //VL
    #endif //not_CUDA

    #ifdef CUDA
    #ifndef VL
    max_dti = CTU_Algorithm_3D_CUDA(g0, g1, H.nx, H.ny, H.nz, x_off, y_off, z_off, H.n_ghost, H.dx, H.dy, H.dz, H.xbound, H.ybound, H.zbound, H.dt, H.n_fields);
    #endif //not_VL
    #ifdef VL
    max_dti = VL_Algorithm_3D_CUDA(g0, g1, H.nx, H.ny, H.nz, x_off, y_off, z_off, H.n_ghost, H.dx, H.dy, H.dz, H.xbound, H.ybound, H.zbound, H.dt, H.n_fields);
    #endif //VL
    #endif    
  }
  else
  {
    chprintf("Error: Grid dimensions nx: %d  ny: %d  nz: %d  not supported.\n", H.nx, H.ny, H.nz);
    chexit(-1);
  }
  // at this point g0 has the old data, g1 has the new data
  // point the grid variables at the new data
  C.density  = &g1[0];
  C.momentum_x = &g1[H.n_cells];
  C.momentum_y = &g1[2*H.n_cells];
  C.momentum_z = &g1[3*H.n_cells];
  C.Energy   = &g1[4*H.n_cells];
  #ifdef SCALAR
  C.scalar = &g1[5*H.n_cells];
  #endif
  #ifdef DE
  C.GasEnergy = &g1[(H.n_fields-1)*H.n_cells];
  #endif

  // reset the grid flag to swap buffers
  gflag = (gflag+1)%2;

  return max_dti;

}


void Grid3D::Fix_Cells(void)
{
  int i, j, k, id;
  Real d, mx, my, mz, P, E;
  Real n, T, mu;
  Real c;
  mu = 0.6;

  for (k=H.n_ghost; k<H.nz-H.n_ghost; k++) {
    for (j=H.n_ghost; j<H.ny-H.n_ghost; j++) {
      for (i=H.n_ghost; i<H.nx-H.n_ghost; i++) {

        id = i + j*H.nx + k*H.nx*H.ny;

        d = C.density[id];
        E = C.Energy[id];
        mx = C.momentum_x[id];
        my = C.momentum_y[id];
        mz = C.momentum_z[id];
        P = (E - (0.5/d)*(mx*mx+ my*my+ mz*mz))*(gama-1.0);
        n = d*DENSITY_UNIT/(mu*MP);
        #ifdef DE
        T = C.GasEnergy[id]*(gama-1.0)*PRESSURE_UNIT/(n*KB); 
        #else
        T = P*PRESSURE_UNIT/(n*KB);
        #endif
        c = C.scalar[id]/d;

        // if there is a problem, replace the cell value with surrounding cell average
        if (d < 0.0 || d != d || P < 0.0 || P != P|| E < 0.0 || E != E|| T > 1.0e9) {
        //if (d < 0.0 || d != d || P < 0.0 || P != P|| E < 0.0 || E != E || T > 1.0e10) {

          printf("%3d %3d %3d BC: d: %e  E:%e  P:%e  n:%e  T:%e  c:%e\n", i+nx_local_start, j+ny_local_start, k+nz_local_start, d, E, P, n, T, c);

          int idn;
          int N = 0;
          int S = 0;
          Real d_av, vx_av, vy_av, vz_av, P_av;
          d_av = vx_av = vy_av = vz_av = P_av = 0.0;
          #ifdef SCALAR
          Real scalar[NSCALARS], scalar_av[NSCALARS];
          for (int n=0; n<NSCALARS; n++) {
            scalar_av[n] = 0.0;
          }
          #endif

          for (int kk=k-1; kk<=k+1; kk++) {
          for (int jj=j-1; jj<=j+1; jj++) {
          for (int ii=i-1; ii<=i+1; ii++) {

            idn = ii + jj*H.nx + kk*H.nx*H.ny; 
            d = C.density[idn];
            mx = C.momentum_x[idn];
            my = C.momentum_y[idn];
            mz = C.momentum_z[idn];
            P  = (C.Energy[idn] - (0.5/d)*(mx*mx + my*my + mz*mz))*(gama-1.0);
            #ifdef SCALAR
            for (int n = 0; n<NSCALARS; n++) {
              scalar[n]  = C.scalar[idn+n*H.n_cells];
            }
            #endif
            if (d > 0.0 && P > 0.0) {
              d_av += d;
              vx_av += mx;
              vy_av += my;
              vz_av += mz;
              P_av += P/(gama-1.0);
              #ifdef SCALAR
              for (int n=0; n<NSCALARS; n++) {
                scalar_av[n] += scalar[n];
              }
              #endif
              N++;
            }

          }
          }
          }

          P_av = P_av / N;
          vx_av = vx_av/d_av;
          vy_av = vy_av/d_av;
          vz_av = vz_av/d_av;
          #ifdef SCALAR
          for (int n=0; n<NSCALARS; n++) {
            scalar_av[n] = scalar_av[n]/d_av;
          }
          #endif
          d_av = d_av/N;

          // replace cell values with new averaged values
          C.density[id] = d_av;
          C.momentum_x[id] = d_av*vx_av;
          C.momentum_y[id] = d_av*vy_av;
          C.momentum_z[id] = d_av*vz_av;
          C.Energy[id] = P_av/(gama-1.0) + 0.5*d_av*(vx_av*vx_av + vy_av*vy_av + vz_av*vz_av);
          #ifdef DE
          C.GasEnergy[id] = P_av/(gama-1.0);
          #endif
          #ifdef SCALAR
          for (int n=0; n<NSCALARS; n++) {
            C.scalar[id+n*H.n_cells] = d_av*scalar_av[n];
          }
          #endif

          d = d_av;
          E = P_av/(gama-1.0) + 0.5*d_av*(vx_av*vx_av + vy_av*vy_av + vz_av*vz_av);
          P = P_av;
          n = d*DENSITY_UNIT/(mu*MP);
          T = P_av*PRESSURE_UNIT/(n*KB);
          c = scalar_av[0];

          printf("%3d %3d %3d FC: d: %e  E:%e  P:%e  n:%e  T:%e  c:%e\n", i+nx_local_start, j+ny_local_start, k+nz_local_start, d, E, P, n, T, c);
          if (d < 0.0 || d != d || P < 0.0 || P != P|| E < 0.0 || E != E|| T > 1.0e9) {
            printf("%3d %3d %3d Flux correction failed: d: %e  E:%e  P:%e  n:%e  T:%e  c:%e\n", i+nx_local_start, j+ny_local_start, k+nz_local_start, d, E, P, n, T, c);
          }

        }
      }
    }
  }

}


/*! \fn void Reset(void)
 *  \brief Reset the Grid3D class. */
void Grid3D::Reset(void)
{
  // free the memory
  FreeMemory();

  // reset the initialization flag
  flag_init = 0;

}


/*! \fn void FreeMemory(void)
 *  \brief Free the memory allocated by the Grid3D class. */
void Grid3D::FreeMemory(void)
{
  // free the conserved variable arrays
  free(buffer0);
  free(buffer1);

  #ifdef COOLING_GPU
  #ifdef CLOUDY_COOL
  Free_Cuda_Textures();
  #endif
  #endif

}
