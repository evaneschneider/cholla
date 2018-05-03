/*! \file hydro_cuda.cu
 *  \brief Definitions of functions used in all cuda integration algorithms. */
#ifdef CUDA

#include<stdio.h>
#include<math.h>
#include<cuda.h>
#include"global.h"
#include"global_cuda.h"
#include"hydro_cuda.h"


__global__ void Update_Conserved_Variables_1D(Real *dev_conserved, Real *dev_F, int n_cells, int x_off, int n_ghost, Real dx, Real xbound, Real dt, Real gamma, int n_fields)
{
  int id;
  #if defined(DE) || defined(STATIC_GRAV)
  Real d, d_inv, vx;  
  #endif
  #ifdef DE
  Real vx_imo, vx_ipo, vy, vz, P;
  #endif
  #ifdef STATIC_GRAV
  Real gx, d_n, d_inv_n, vx_n;
  gx = 0.0;
  #endif
  
  Real dtodx = dt/dx;

  // get a global thread ID
  id = threadIdx.x + blockIdx.x * blockDim.x;


  // threads corresponding to real cells do the calculation
  if (id > n_ghost - 1 && id < n_cells-n_ghost)
  {
    #if defined(DE) || defined(STATIC_GRAV)
    d  =  dev_conserved[            id];
    d_inv = 1.0 / d;
    vx =  dev_conserved[1*n_cells + id] * d_inv;
    #endif
    #ifdef DE
    vy =  dev_conserved[2*n_cells + id] * d_inv;
    vz =  dev_conserved[3*n_cells + id] * d_inv;
    P  = (dev_conserved[4*n_cells + id] - 0.5*d*(vx*vx + vy*vy + vz*vz)) * (gamma - 1.0);
    vx_imo = dev_conserved[1*n_cells + id-1]/dev_conserved[id-1];
    vx_ipo = dev_conserved[1*n_cells + id+1]/dev_conserved[id+1];
    #endif
  
    // update the conserved variable array
    dev_conserved[            id] += dtodx * (dev_F[            id-1] - dev_F[            id]);
    dev_conserved[  n_cells + id] += dtodx * (dev_F[  n_cells + id-1] - dev_F[  n_cells + id]);
    dev_conserved[2*n_cells + id] += dtodx * (dev_F[2*n_cells + id-1] - dev_F[2*n_cells + id]);
    dev_conserved[3*n_cells + id] += dtodx * (dev_F[3*n_cells + id-1] - dev_F[3*n_cells + id]);
    dev_conserved[4*n_cells + id] += dtodx * (dev_F[4*n_cells + id-1] - dev_F[4*n_cells + id]);
    #ifdef SCALAR
    for (int i=0; i<NSCALARS; i++) {
      dev_conserved[(5+i)*n_cells + id] += dtodx * (dev_F[(5+i)*n_cells + id-1] - dev_F[(5+i)*n_cells + id]);
    }
    #endif
    #ifdef DE
    dev_conserved[(n_fields-1)*n_cells + id] += dtodx * (dev_F[(n_fields-1)*n_cells + id-1] - dev_F[(n_fields-1)*n_cells + id])
                                  +  dtodx * P * 0.5 * (vx_imo - vx_ipo);
    #endif
    #ifdef STATIC_GRAV // add gravitational source terms, time averaged from n to n+1
    calc_g_1D(id, x_off, n_ghost, dx, xbound, &gx);    
    d_n  =  dev_conserved[            id];
    d_inv_n = 1.0 / d_n;
    vx_n =  dev_conserved[1*n_cells + id] * d_inv_n;
    dev_conserved[  n_cells + id] += 0.5*dt*gx*(d + d_n);
    dev_conserved[4*n_cells + id] += 0.25*dt*gx*(d + d_n)*(vx + vx_n);
    #endif    
    //printf("%3d %e %e %e\n", id-n_ghost, fz, gcorr, gx);
    if (dev_conserved[id] != dev_conserved[id]) printf("%3d Thread crashed in final update. %f\n", id, dev_conserved[id]);
    /*
    d  =  dev_conserved[            id];
    d_inv = 1.0 / d;
    vx =  dev_conserved[1*n_cells + id] * d_inv;
    vy =  dev_conserved[2*n_cells + id] * d_inv;
    vz =  dev_conserved[3*n_cells + id] * d_inv;
    P  = (dev_conserved[4*n_cells + id] - 0.5*d*(vx*vx + vy*vy + vz*vz)) * (gamma - 1.0);
    if (P < 0.0) printf("%d Negative pressure after final update.\n", id);
    */
  }


}


__global__ void Update_Conserved_Variables_2D(Real *dev_conserved, Real *dev_F_x, Real *dev_F_y, int nx, int ny, int x_off, int y_off, int n_ghost, Real dx, Real dy, Real xbound, Real ybound, Real dt, Real gamma, int n_fields)
{
  int id, xid, yid, n_cells;
  int imo, jmo;

  #if defined (DE) || defined(STATIC_GRAV)
  Real d, d_inv, vx, vy;
  #endif
  #ifdef DE
  Real vx_imo, vx_ipo, vy_jmo, vy_jpo, vz, P;
  int ipo, jpo;
  #endif

  #ifdef STATIC_GRAV
  Real gx, gy, d_n, d_inv_n, vx_n, vy_n;
  //gx = 0.0;
  //gy = 0.0;
  #endif

  Real dtodx = dt/dx;
  Real dtody = dt/dy;

  n_cells = nx*ny;

  // get a global thread ID
  int blockId = blockIdx.x + blockIdx.y*gridDim.x;
  id = threadIdx.x + blockId * blockDim.x;
  yid = id / nx;
  xid = id - yid*nx;
  imo = xid-1 + yid*nx;
  jmo = xid + (yid-1)*nx;

  // threads corresponding to real cells do the calculation
  if (xid > n_ghost-1 && xid < nx-n_ghost && yid > n_ghost-1 && yid < ny-n_ghost)
  {
    #if defined (DE) || defined (STATIC_GRAV)
    d  =  dev_conserved[            id];
    d_inv = 1.0 / d;
    vx =  dev_conserved[1*n_cells + id] * d_inv;
    vy =  dev_conserved[2*n_cells + id] * d_inv;
    #endif
    #ifdef DE
    vz =  dev_conserved[3*n_cells + id] * d_inv;
    P  = (dev_conserved[4*n_cells + id] - 0.5*d*(vx*vx + vy*vy + vz*vz)) * (gamma - 1.0);
    ipo = xid+1 + yid*nx;
    jpo = xid + (yid+1)*nx;
    vx_imo = dev_conserved[1*n_cells + imo] / dev_conserved[imo]; 
    vx_ipo = dev_conserved[1*n_cells + ipo] / dev_conserved[ipo]; 
    vy_jmo = dev_conserved[2*n_cells + jmo] / dev_conserved[jmo]; 
    vy_jpo = dev_conserved[2*n_cells + jpo] / dev_conserved[jpo]; 
    #endif
    // update the conserved variable array
    dev_conserved[            id] += dtodx * (dev_F_x[            imo] - dev_F_x[            id])
                                  +  dtody * (dev_F_y[            jmo] - dev_F_y[            id]);
    dev_conserved[  n_cells + id] += dtodx * (dev_F_x[  n_cells + imo] - dev_F_x[  n_cells + id]) 
                                  +  dtody * (dev_F_y[  n_cells + jmo] - dev_F_y[  n_cells + id]);
    dev_conserved[2*n_cells + id] += dtodx * (dev_F_x[2*n_cells + imo] - dev_F_x[2*n_cells + id]) 
                                  +  dtody * (dev_F_y[2*n_cells + jmo] - dev_F_y[2*n_cells + id]); 
    dev_conserved[3*n_cells + id] += dtodx * (dev_F_x[3*n_cells + imo] - dev_F_x[3*n_cells + id])
                                  +  dtody * (dev_F_y[3*n_cells + jmo] - dev_F_y[3*n_cells + id]);
    dev_conserved[4*n_cells + id] += dtodx * (dev_F_x[4*n_cells + imo] - dev_F_x[4*n_cells + id])
                                  +  dtody * (dev_F_y[4*n_cells + jmo] - dev_F_y[4*n_cells + id]);
    #ifdef SCALAR
    for (int i=0; i<NSCALARS; i++) {
      dev_conserved[(5+i)*n_cells + id] += dtodx * (dev_F_x[(5+i)*n_cells + imo] - dev_F_x[(5+i)*n_cells + id])
                                        +  dtody * (dev_F_y[(5+i)*n_cells + jmo] - dev_F_y[(5+i)*n_cells + id]);
    }
//    if (xid == 4) {
//      printf("%3d %f %f %f\n", yid, dev_conserved[5*n_cells + id]), dtodx*(dev_F_x[(5+i)*n_cells + imo] - dev_F_x[(5+i)*n_cells + id]), dtody * (dev_F_y[(5+i)*n_cells + jmo] - dev_F_y[(5+i)*n_cells + id]);
//    }
    #endif
    // add gravitational source terms, time averaged from n to n+1                                 
    #ifdef DE
    dev_conserved[(n_fields-1)*n_cells + id] += dtodx * (dev_F_x[(n_fields-1)*n_cells + imo] - dev_F_x[(n_fields-1)*n_cells + id])
                                  +  dtody * (dev_F_y[(n_fields-1)*n_cells + jmo] - dev_F_y[(n_fields-1)*n_cells + id])
                                  +  0.5*P*(dtodx*(vx_imo-vx_ipo) + dtody*(vy_jmo-vy_jpo));
    //if (dev_conserved[5*n_cells + id] < 0.0) printf("%3d %3d Negative internal energy after final update.\n", xid, yid);
    #endif
    #ifdef STATIC_GRAV 
    // calculate the gravitational acceleration as a function of x & y position
    calc_g_2D(xid, yid, x_off, y_off, n_ghost, dx, dy, xbound, ybound, &gx, &gy);
    d_n  =  dev_conserved[            id];
    d_inv_n = 1.0 / d_n;
    vx_n =  dev_conserved[1*n_cells + id] * d_inv_n;
    vy_n =  dev_conserved[2*n_cells + id] * d_inv_n;
    dev_conserved[  n_cells + id] += 0.5*dt*gx*(d + d_n);
    dev_conserved[2*n_cells + id] += 0.5*dt*gy*(d + d_n);
    dev_conserved[4*n_cells + id] += 0.25*dt*gx*(d + d_n)*(vx + vx_n)
                                  +  0.25*dt*gy*(d + d_n)*(vy + vy_n);
    //dev_conserved[4*n_cells + id] += 0.5*dt*gx*(d*vx + d_n*vx_n)
    //                              +  0.5*dt*gy*(d*vy + d_n*vy_n);
    #endif
    if (dev_conserved[id] < 0.0 || dev_conserved[id] != dev_conserved[id]) {
      printf("%3d %3d Thread crashed in final update. %f %f %f\n", xid, yid, dtodx*(dev_F_x[imo]-dev_F_x[id]), dtody*(dev_F_y[jmo]-dev_F_y[id]), dev_conserved[id]);
    }   
    /*
    d  =  dev_conserved[            id];
    d_inv = 1.0 / d;
    vx =  dev_conserved[1*n_cells + id] * d_inv;
    vy =  dev_conserved[2*n_cells + id] * d_inv;
    vz =  dev_conserved[3*n_cells + id] * d_inv;
    P  = (dev_conserved[4*n_cells + id] - 0.5*d*(vx*vx + vy*vy + vz*vz)) * (gamma - 1.0);
    if (P < 0.0)
      printf("%3d %3d Negative pressure after final update. %f %f %f %f\n", xid, yid, dev_conserved[4*n_cells + id], 0.5*d*vx*vx, 0.5*d*vy*vy, P);    
    */
  }

}



__global__ void Update_Conserved_Variables_3D(Real *dev_conserved, Real *dev_F_x, Real *dev_F_y,  Real *dev_F_z,
                                              int nx, int ny, int nz, int x_off, int y_off, int z_off, int n_ghost, 
                                              Real dx, Real dy, Real dz, Real xbound, Real ybound, Real zbound, Real dt,
                                              Real gamma, int n_fields)
{
  int id, xid, yid, zid, n_cells;
  int imo, jmo, kmo;
  #if defined (DE) || defined(STATIC_GRAV)
  Real d, d_inv, vx, vy, vz;
  #endif
  #ifdef DE
  Real vx_imo, vx_ipo, vy_jmo, vy_jpo, vz_kmo, vz_kpo, P;
  int ipo, jpo, kpo;
  #endif

  #ifdef STATIC_GRAV
  Real gx, gy, gz, d_n, d_inv_n, vx_n, vy_n, vz_n;
  gx = 0.0;
  gy = 0.0;
  gz = 0.0;
  #endif

  Real dtodx = dt/dx;
  Real dtody = dt/dy;
  Real dtodz = dt/dz;
  n_cells = nx*ny*nz;

  // get a global thread ID
  id = threadIdx.x + blockIdx.x * blockDim.x;
  zid = id / (nx*ny);
  yid = (id - zid*nx*ny) / nx;
  xid = id - zid*nx*ny - yid*nx;
  imo = xid-1 + yid*nx + zid*nx*ny;
  jmo = xid + (yid-1)*nx + zid*nx*ny;
  kmo = xid + yid*nx + (zid-1)*nx*ny;

  // threads corresponding to real cells do the calculation
  if (xid > n_ghost-1 && xid < nx-n_ghost && yid > n_ghost-1 && yid < ny-n_ghost && zid > n_ghost-1 && zid < nz-n_ghost)
  {
    #if defined (DE) || defined(STATIC_GRAV)
    d  =  dev_conserved[            id];
    d_inv = 1.0 / d;
    vx =  dev_conserved[1*n_cells + id] * d_inv;
    vy =  dev_conserved[2*n_cells + id] * d_inv;
    vz =  dev_conserved[3*n_cells + id] * d_inv;
    #endif
    #ifdef DE
    P  = (dev_conserved[4*n_cells + id] - 0.5*d*(vx*vx + vy*vy + vz*vz)) * (gamma - 1.0);
    //if (d < 0.0 || d != d) printf("Negative density before final update.\n");
    //if (P < 0.0) printf("%d Negative pressure before final update.\n", id);
    ipo = xid+1 + yid*nx + zid*nx*ny;
    jpo = xid + (yid+1)*nx + zid*nx*ny;
    kpo = xid + yid*nx + (zid+1)*nx*ny;
    vx_imo = dev_conserved[1*n_cells + imo] / dev_conserved[imo]; 
    vx_ipo = dev_conserved[1*n_cells + ipo] / dev_conserved[ipo]; 
    vy_jmo = dev_conserved[2*n_cells + jmo] / dev_conserved[jmo]; 
    vy_jpo = dev_conserved[2*n_cells + jpo] / dev_conserved[jpo]; 
    vz_kmo = dev_conserved[3*n_cells + kmo] / dev_conserved[kmo]; 
    vz_kpo = dev_conserved[3*n_cells + kpo] / dev_conserved[kpo]; 
    #endif

    // update the conserved variable array
    dev_conserved[            id] += dtodx * (dev_F_x[            imo] - dev_F_x[            id])
                                  +  dtody * (dev_F_y[            jmo] - dev_F_y[            id])
                                  +  dtodz * (dev_F_z[            kmo] - dev_F_z[            id]);
    dev_conserved[  n_cells + id] += dtodx * (dev_F_x[  n_cells + imo] - dev_F_x[  n_cells + id])
                                  +  dtody * (dev_F_y[  n_cells + jmo] - dev_F_y[  n_cells + id])
                                  +  dtodz * (dev_F_z[  n_cells + kmo] - dev_F_z[  n_cells + id]);
    dev_conserved[2*n_cells + id] += dtodx * (dev_F_x[2*n_cells + imo] - dev_F_x[2*n_cells + id])
                                  +  dtody * (dev_F_y[2*n_cells + jmo] - dev_F_y[2*n_cells + id])
                                  +  dtodz * (dev_F_z[2*n_cells + kmo] - dev_F_z[2*n_cells + id]);
    dev_conserved[3*n_cells + id] += dtodx * (dev_F_x[3*n_cells + imo] - dev_F_x[3*n_cells + id])
                                  +  dtody * (dev_F_y[3*n_cells + jmo] - dev_F_y[3*n_cells + id])
                                  +  dtodz * (dev_F_z[3*n_cells + kmo] - dev_F_z[3*n_cells + id]);
    //fz = dtodx * (dev_F_x[3*n_cells + imo] - dev_F_x[3*n_cells + id])
    //                              +  dtody * (dev_F_y[3*n_cells + jmo] - dev_F_y[3*n_cells + id])
    //                              +  dtodz * (dev_F_z[3*n_cells + kmo] - dev_F_z[3*n_cells + id]);
    dev_conserved[4*n_cells + id] += dtodx * (dev_F_x[4*n_cells + imo] - dev_F_x[4*n_cells + id])
                                  +  dtody * (dev_F_y[4*n_cells + jmo] - dev_F_y[4*n_cells + id])
                                  +  dtodz * (dev_F_z[4*n_cells + kmo] - dev_F_z[4*n_cells + id]);
    #ifdef SCALAR
    for (int i=0; i<NSCALARS; i++) {
      dev_conserved[(5+i)*n_cells + id] += dtodx * (dev_F_x[(5+i)*n_cells + imo] - dev_F_x[(5+i)*n_cells + id])
                                    +  dtody * (dev_F_y[(5+i)*n_cells + jmo] - dev_F_y[(5+i)*n_cells + id])
                                    +  dtodz * (dev_F_z[(5+i)*n_cells + kmo] - dev_F_z[(5+i)*n_cells + id]);
    }                              
    #endif
    #ifdef DE
    dev_conserved[(n_fields-1)*n_cells + id] += dtodx * (dev_F_x[(n_fields-1)*n_cells + imo] - dev_F_x[(n_fields-1)*n_cells + id])
                                  +  dtody * (dev_F_y[(n_fields-1)*n_cells + jmo] - dev_F_y[(n_fields-1)*n_cells + id])
                                  +  dtodz * (dev_F_z[(n_fields-1)*n_cells + kmo] - dev_F_z[(n_fields-1)*n_cells + id])
                                  +  0.5*P*(dtodx*(vx_imo-vx_ipo) + dtody*(vy_jmo-vy_jpo) + dtodz*(vz_kmo-vz_kpo));
    #endif
    // add a density floor of n = 1.0e-3
    //dev_conserved[id] = fmax(dev_conserved[id], 0.6*MP*1.0e-3/DENSITY_UNIT);
    #ifdef STATIC_GRAV 
    calc_g_3D_CUDA(xid, yid, zid, x_off, y_off, z_off, n_ghost, dx, dy, dz, xbound, ybound, zbound, &gx, &gy, &gz);
    d_n  =  dev_conserved[            id];
    d_inv_n = 1.0 / d_n;
    vx_n =  dev_conserved[1*n_cells + id] * d_inv_n;
    vy_n =  dev_conserved[2*n_cells + id] * d_inv_n;
    vz_n =  dev_conserved[3*n_cells + id] * d_inv_n;
    dev_conserved[  n_cells + id] += 0.5*dt*gx*(d + d_n);
    dev_conserved[2*n_cells + id] += 0.5*dt*gy*(d + d_n);
    dev_conserved[3*n_cells + id] += 0.5*dt*gz*(d + d_n);
    dev_conserved[4*n_cells + id] += 0.25*dt*gx*(d + d_n)*(vx + vx_n)
                                  +  0.25*dt*gy*(d + d_n)*(vy + vy_n)
                                  +  0.25*dt*gz*(d + d_n)*(vz + vz_n);
    #endif    
    
    if (dev_conserved[id] < 0.0 || dev_conserved[id] != dev_conserved[id] || dev_conserved[4*n_cells + id] < 0.0 || dev_conserved[4*n_cells+id] != dev_conserved[4*n_cells+id]) {
      printf("%3d %3d %3d Thread crashed in final update. %e %e %e %e %e\n", xid+x_off, yid+y_off, zid+z_off, dev_conserved[id], dtodx*(dev_F_x[imo]-dev_F_x[id]), dtody*(dev_F_y[jmo]-dev_F_y[id]), dtodz*(dev_F_z[kmo]-dev_F_z[id]), dev_conserved[4*n_cells+id]);
      //printf("%3d %3d %3d Thread crashed in final update. %e %e\n", xid+x_off, yid+y_off, zid+z_off, dev_conserved[id], dev_conserved[4*n_cells+id]);
      //Real ge = dev_conserved[5*n_cells + id];
      //Real T = ge * (gamma-1.0)*SP_ENERGY_UNIT*0.6*MP/(d_n*KB);
      //printf("Internal energy: %e  Temperature: %e\n", ge, T);
    }
    
    /*
    d  =  dev_conserved[            id];
    d_inv = 1.0 / d;
    vx =  dev_conserved[1*n_cells + id] * d_inv;
    vy =  dev_conserved[2*n_cells + id] * d_inv;
    vz =  dev_conserved[3*n_cells + id] * d_inv;
    P  = (dev_conserved[4*n_cells + id] - 0.5*d*(vx*vx + vy*vy + vz*vz)) * (gamma - 1.0);
    if (P < 0.0) printf("%3d %3d %3d Negative pressure after final update. %f %f %f %f %f\n", xid, yid, zid, dev_conserved[4*n_cells + id], 0.5*d*vx*vx, 0.5*d*vy*vy, 0.5*d*vz*vz, P);
    */
  }

}


__global__ void Sync_Energies_1D(Real *dev_conserved, int n_cells, int n_ghost, Real gamma, int n_fields)
{
  int id;
  Real d, d_inv, vx, vy, vz, P, E;
  Real ge1, ge2, Emax;
  int im1, ip1;

  // get a global thread ID
  id = threadIdx.x + blockIdx.x * blockDim.x;
  
  im1 = max(id-1, n_ghost);
  ip1 = min(id+1, n_cells-n_ghost-1);

  // threads corresponding to real cells do the calculation
  if (id > n_ghost - 1 && id < n_cells-n_ghost)
  {
    // every thread collects the conserved variables it needs from global memory
    d  =  dev_conserved[            id];
    d_inv = 1.0 / d;
    vx =  dev_conserved[1*n_cells + id] * d_inv;
    vy =  dev_conserved[2*n_cells + id] * d_inv;
    vz =  dev_conserved[3*n_cells + id] * d_inv;
    E  =  dev_conserved[4*n_cells + id];
    P  = (dev_conserved[4*n_cells + id] - 0.5*d*(vx*vx + vy*vy + vz*vz)) * (gamma - 1.0);
    // separately tracked internal energy 
    ge1 = dev_conserved[(n_fields-1)*n_cells + id];
    // internal energy calculated from total energy
    ge2 = dev_conserved[4*n_cells + id] - 0.5*d*(vx*vx + vy*vy + vz*vz);
    // if the ratio of conservatively calculated internal energy to total energy
    // is greater than 1/1000, use the conservatively calculated internal energy
    // to do the internal energy update
    if (ge2/E > 0.001) {
      dev_conserved[(n_fields-1)*n_cells + id] = ge2;
      ge1 = ge2;
    }     
    // find the max nearby total energy 
    Emax = fmax(dev_conserved[4*n_cells + im1], E);
    Emax = fmax(dev_conserved[4*n_cells + ip1], Emax);
    // if the ratio of conservatively calculated internal energy to max nearby total energy
    // is greater than 1/10, continue to use the conservatively calculated internal energy 
    if (ge2/Emax > 0.1) {
      dev_conserved[(n_fields-1)*n_cells + id] = ge2;
    }
    // sync the total energy with the internal energy 
    else {
      dev_conserved[4*n_cells + id] += ge1 - ge2;
    }
    /*
    // if the conservatively calculated internal energy is greater than the estimate of the truncation error,
    // use the internal energy computed from the total energy to do the update
    //find the max nearby velocity difference (estimate of truncation error) 
    vmax = fmax(fabs(vx-dev_conserved[1*n_cells + im1]/dev_conserved[im1]), fabs(dev_conserved[1*n_cells + ip1]/dev_conserved[ip1]-vx));
    //printf("%3d %f %f %f %f\n", id, ge1, ge2, vmax, 0.25*d*vmax*vmax);
    if (ge2 > 0.25*d*vmax*vmax) {
      dev_conserved[5*n_cells + id] = ge2;
      ge1 = ge2;
    }
    //else printf("%d Using ge1 %f %f %f %f\n", id, ge1, ge2, vmax, 0.25*d*vmax*vmax);
    */
    // update the total energy
     
    // recalculate the pressure 
    P = (dev_conserved[4*n_cells + id] - 0.5*d*(vx*vx + vy*vy + vz*vz)) * (gamma - 1.0);    
    if (P < 0.0) printf("%d Negative pressure after internal energy sync. %f %f \n", id, ge1, ge2);    
  }

}


__global__ void Sync_Energies_2D(Real *dev_conserved, int nx, int ny, int n_ghost, Real gamma, int n_fields)
{
  int id, xid, yid, n_cells;
  Real d, d_inv, vx, vy, vz, P, E;
  Real ge1, ge2, Emax;
  int imo, ipo, jmo, jpo;
  n_cells = nx*ny;

  // get a global thread ID
  int blockId = blockIdx.x + blockIdx.y*gridDim.x;
  id = threadIdx.x + blockId * blockDim.x;
  yid = id / nx;
  xid = id - yid*nx;

  imo = max(xid-1, n_ghost);
  imo = imo + yid*nx;
  ipo = min(xid+1, nx-n_ghost-1);
  ipo = ipo + yid*nx;
  jmo = max(yid-1, n_ghost);
  jmo = xid + jmo*nx;
  jpo = min(yid+1, ny-n_ghost-1);
  jpo = xid + jpo*nx;

  // threads corresponding to real cells do the calculation
  if (xid > n_ghost-1 && xid < nx-n_ghost && yid > n_ghost-1 && yid < ny-n_ghost)
  {
    // every thread collects the conserved variables it needs from global memory
    d  =  dev_conserved[            id];
    d_inv = 1.0 / d;
    vx =  dev_conserved[1*n_cells + id] * d_inv;
    vy =  dev_conserved[2*n_cells + id] * d_inv;
    vz =  dev_conserved[3*n_cells + id] * d_inv;
    E  =  dev_conserved[4*n_cells + id];
    P  = (E - 0.5*d*(vx*vx + vy*vy + vz*vz)) * (gamma - 1.0);
    // separately tracked internal energy 
    ge1 =  dev_conserved[(n_fields-1)*n_cells + id];
    // internal energy calculated from total energy
    ge2 = dev_conserved[4*n_cells + id] - 0.5*d*(vx*vx + vy*vy + vz*vz);
    // if the ratio of conservatively calculated internal energy to total energy
    // is greater than 1/1000, use the conservatively calculated internal energy
    // to do the internal energy update
    if (ge2/E > 0.001) {
      dev_conserved[(n_fields-1)*n_cells + id] = ge2;
      ge1 = ge2;
    }     
    //find the max nearby total energy 
    Emax = fmax(dev_conserved[4*n_cells + imo], E);
    Emax = fmax(Emax, dev_conserved[4*n_cells + ipo]);
    Emax = fmax(Emax, dev_conserved[4*n_cells + jmo]);
    Emax = fmax(Emax, dev_conserved[4*n_cells + jpo]);
    // if the ratio of conservatively calculated internal energy to max nearby total energy
    // is greater than 1/10, continue to use the conservatively calculated internal energy 
    if (ge2/Emax > 0.1) {
      dev_conserved[(n_fields-1)*n_cells + id] = ge2;
    }
    // sync the total energy with the internal energy 
    else {
      dev_conserved[4*n_cells + id] += ge1 - ge2;
    }
    // recalculate the pressure 
    P = (dev_conserved[4*n_cells + id] - 0.5*d*(vx*vx + vy*vy + vz*vz)) * (gamma - 1.0);    
    if (P < 0.0) printf("%d Negative pressure after internal energy sync. %f %f \n", id, ge1, ge2);    
  }
}




__global__ void Sync_Energies_3D(Real *dev_conserved, int nx, int ny, int nz, int n_ghost, Real gamma, int n_fields)
{
  int id, xid, yid, zid, n_cells;
  Real d, d_inv, vx, vy, vz, E;
  Real ge1, ge2, Emax;
  int imo, ipo, jmo, jpo, kmo, kpo;
  n_cells = nx*ny*nz;

  // get a global thread ID
  id = threadIdx.x + blockIdx.x * blockDim.x;
  zid = id / (nx*ny);
  yid = (id - zid*nx*ny) / nx;
  xid = id - zid*nx*ny - yid*nx;

  imo = max(xid-1, n_ghost);
  imo = imo + yid*nx + zid*nx*ny;
  ipo = min(xid+1, nx-n_ghost-1);
  ipo = ipo + yid*nx + zid*nx*ny;
  jmo = max(yid-1, n_ghost);
  jmo = xid + jmo*nx + zid*nx*ny;
  jpo = min(yid+1, ny-n_ghost-1);
  jpo = xid + jpo*nx + zid*nx*ny;
  kmo = max(zid-1, n_ghost);
  kmo = xid + yid*nx + kmo*nx*ny;
  kpo = min(zid+1, nz-n_ghost-1);
  kpo = xid + yid*nx + kpo*nx*ny;

  // threads corresponding to real cells do the calculation
  if (xid > n_ghost-1 && xid < nx-n_ghost && yid > n_ghost-1 && yid < ny-n_ghost && zid > n_ghost-1 && zid < nz-n_ghost)
  {
    // every thread collects the conserved variables it needs from global memory
    d  =  dev_conserved[            id];
    d_inv = 1.0 / d;
    vx =  dev_conserved[1*n_cells + id] * d_inv;
    vy =  dev_conserved[2*n_cells + id] * d_inv;
    vz =  dev_conserved[3*n_cells + id] * d_inv;
    E  =  dev_conserved[4*n_cells + id];
    // don't do the energy sync if this thread has crashed
    if (E < 0.0 || E != E) return;
    // separately tracked internal energy 
    ge1 =  dev_conserved[(n_fields-1)*n_cells + id];
    // internal energy calculated from total energy
    ge2 = dev_conserved[4*n_cells + id] - 0.5*d*(vx*vx + vy*vy + vz*vz);
    /*
    if (ge2 < 0.0) {
      Real T = (ge1/d) * (gamma-1.0)*SP_ENERGY_UNIT*0.6*MP/KB;
      printf("%3d %3d %3d Temperature: %e\n", xid, yid, zid, T);
    }
    */
    // if the ratio of conservatively calculated internal energy to total energy
    // is greater than 1/1000, use the conservatively calculated internal energy
    // to do the internal energy update
    if (ge2 > 0.0 && E > 0.0 && ge2/E > 0.001) {
      dev_conserved[(n_fields-1)*n_cells + id] = ge2;
      ge1 = ge2;
    }
    //find the max nearby total energy 
    Emax = fmax(dev_conserved[4*n_cells + imo], E);
    Emax = fmax(Emax, dev_conserved[4*n_cells + ipo]);
    Emax = fmax(Emax, dev_conserved[4*n_cells + jmo]);
    Emax = fmax(Emax, dev_conserved[4*n_cells + jpo]);
    Emax = fmax(Emax, dev_conserved[4*n_cells + kmo]);
    Emax = fmax(Emax, dev_conserved[4*n_cells + kpo]);
    // if the ratio of conservatively calculated internal energy to max nearby total energy
    // is greater than 1/10, continue to use the conservatively calculated internal energy 
    if (ge2/Emax > 0.1 && ge2 > 0.0 && Emax > 0.0) {
      dev_conserved[(n_fields-1)*n_cells + id] = ge2;
    }
    // sync the total energy with the internal energy 
    else {
      if (ge1 > 0.0) dev_conserved[4*n_cells + id] += ge1 - ge2;
      else dev_conserved[(n_fields-1)*n_cells+id] = ge2;
    }
    // recalculate the pressure 
    //Real P = (dev_conserved[4*n_cells + id] - 0.5*d*(vx*vx + vy*vy + vz*vz)) * (gamma - 1.0);
    //if (P < 0.0) printf("%3d %3d %3d Negative pressure after internal energy sync. %f %f %f\n", xid, yid, zid, P/(gamma-1.0), ge1, ge2);    
  }
}



__global__ void Calc_dt_1D(Real *dev_conserved, int n_cells, int n_ghost, Real dx, Real *dti_array, Real gamma)
{
  __shared__ Real max_dti[TPB];

  Real d, d_inv, vx, vy, vz, P, cs;
  int id, tid;

  // get a global thread ID
  id = threadIdx.x + blockIdx.x * blockDim.x;
  // and a thread id within the block
  tid = threadIdx.x;

  // set shared memory to 0
  max_dti[tid] = 0;
  __syncthreads();


  // threads corresponding to real cells do the calculation
  if (id > n_ghost - 1 && id < n_cells-n_ghost)
  {
    // start timestep calculation here
    // every thread collects the conserved variables it needs from global memory
    d  =  dev_conserved[            id];
    d_inv = 1.0 / d;
    vx =  dev_conserved[1*n_cells + id] * d_inv;
    vy =  dev_conserved[2*n_cells + id] * d_inv;
    vz =  dev_conserved[3*n_cells + id] * d_inv;
    P  = (dev_conserved[4*n_cells + id] - 0.5*d*(vx*vx + vy*vy + vz*vz)) * (gamma - 1.0);
    P  = fmax(P, (Real) TINY_NUMBER);
    // find the max wavespeed in that cell, use it to calculate the inverse timestep
    cs = sqrt(d_inv * gamma * P);
    max_dti[tid] = (fabs(vx)+cs)/dx;
  }
  __syncthreads();
  
  // do the reduction in shared memory (find the max inverse timestep in the block)
  for (unsigned int s=1; s<blockDim.x; s*=2) {
    if (tid % (2*s) == 0) {
      max_dti[tid] = fmax(max_dti[tid], max_dti[tid + s]);
    }
    __syncthreads();
  }

  // write the result for this block to global memory
  if (tid == 0) dti_array[blockIdx.x] = max_dti[0];


}



__global__ void Calc_dt_2D(Real *dev_conserved, int nx, int ny, int n_ghost, Real dx, Real dy, Real *dti_array, Real gamma)
{
  __shared__ Real max_dti[TPB];

  Real d, d_inv, vx, vy, vz, P, cs;
  int id, tid, xid, yid, n_cells;
  n_cells = nx*ny;

  // get a global thread ID
  int blockId = blockIdx.x + blockIdx.y*gridDim.x;
  id = threadIdx.x + blockId * blockDim.x;
  yid = id / nx;
  xid = id - yid*nx;
  // and a thread id within the block
  tid = threadIdx.x;

  // set shared memory to 0
  max_dti[tid] = 0;
  __syncthreads();

  // threads corresponding to real cells do the calculation
  if (xid > n_ghost-1 && xid < nx-n_ghost && yid > n_ghost-1 && yid < ny-n_ghost)
  {
    // every thread collects the conserved variables it needs from global memory
    d  =  dev_conserved[            id];
    d_inv = 1.0 / d;
    vx =  dev_conserved[1*n_cells + id] * d_inv;
    vy =  dev_conserved[2*n_cells + id] * d_inv;
    vz =  dev_conserved[3*n_cells + id] * d_inv;
    P  = (dev_conserved[4*n_cells + id] - 0.5*d*(vx*vx + vy*vy + vz*vz)) * (gamma - 1.0);
    P  = fmax(P, (Real) 1.0e-20);
    // find the max wavespeed in that cell, use it to calculate the inverse timestep
    cs = sqrt(d_inv * gamma * P);
    max_dti[tid] = fmax((fabs(vx)+cs)/dx, (fabs(vy)+cs)/dy);
  }
  __syncthreads();
  
  // do the reduction in shared memory (find the max inverse timestep in the block)
  for (unsigned int s=1; s<blockDim.x; s*=2) {
    if (tid % (2*s) == 0) {
      max_dti[tid] = fmax(max_dti[tid], max_dti[tid + s]);
    }
    __syncthreads();
  }

  // write the result for this block to global memory
  if (tid == 0) dti_array[blockId] = max_dti[0];

}


__global__ void Calc_dt_3D(Real *dev_conserved, int nx, int ny, int nz, int n_ghost, Real dx, Real dy, Real dz, Real *dti_array, Real gamma)
{
  __shared__ Real max_dti[TPB];

  Real d, d_inv, vx, vy, vz, E, P, cs;
  int id, xid, yid, zid, n_cells;
  int tid;

  n_cells = nx*ny*nz;

  // get a global thread ID
  id = threadIdx.x + blockIdx.x * blockDim.x;
  zid = id / (nx*ny);
  yid = (id - zid*nx*ny) / nx;
  xid = id - zid*nx*ny - yid*nx;
  // and a thread id within the block  
  tid = threadIdx.x;

  // set shared memory to 0
  max_dti[tid] = 0;
  __syncthreads();

  // threads corresponding to real cells do the calculation
  if (xid > n_ghost-1 && xid < nx-n_ghost && yid > n_ghost-1 && yid < ny-n_ghost && zid > n_ghost-1 && zid < nz-n_ghost)
  {
    // every thread collects the conserved variables it needs from global memory
    d  =  dev_conserved[            id];
    d_inv = 1.0 / d;
    vx =  dev_conserved[1*n_cells + id] * d_inv;
    vy =  dev_conserved[2*n_cells + id] * d_inv;
    vz =  dev_conserved[3*n_cells + id] * d_inv;
    E  = dev_conserved[4*n_cells + id];
    P  = (E - 0.5*d*(vx*vx + vy*vy + vz*vz)) * (gamma - 1.0);
    cs = sqrt(d_inv * gamma * P);
    max_dti[tid] = fmax((fabs(vx)+cs)/dx, (fabs(vy)+cs)/dy);
    max_dti[tid] = fmax(max_dti[tid], (fabs(vz)+cs)/dz);
    max_dti[tid] = fmax(max_dti[tid], 0);
    Real n = d*DENSITY_UNIT / (1.27*MP);
    Real T = P*PRESSURE_UNIT/(n*KB);
    // if this thread has crashed or will be reset, don't include it in the timestep calculation
    if (d < 0 || d != d || P < 0 || P != P || E < 0 || E != E) max_dti[tid] = 0;
    //max_dti[tid] = (fabs(vx)+cs)/dx + (fabs(vy)+cs)/dy + (fabs(vz)+cs)/dz;
    /*
    if (0.3/max_dti[tid] < 1.0) {
      Real n = d*DENSITY_UNIT/(0.6*MP);
      Real T = P*PRESSURE_UNIT/(n*KB);
      printf("%3d %3d %3d n: %e  T: %e\n", xid, yid, zid, n, T);
    }
    */
  }
  __syncthreads();
  
  // do the reduction in shared memory (find the max inverse timestep in the block)
  for (unsigned int s=1; s<blockDim.x; s*=2) {
    if (tid % (2*s) == 0) {
      max_dti[tid] = fmax(max_dti[tid], max_dti[tid + s]);
    }
    __syncthreads();
  }

  // write the result for this block to global memory
  if (tid == 0) dti_array[blockIdx.x] = max_dti[0];

}


__device__ void calc_g_1D(int xid, int x_off, int n_ghost, Real dx, Real xbound, Real *gx)
{
  Real x_pos, r_disk, r_halo;
  x_pos = (x_off + xid - n_ghost + 0.5)*dx + xbound;

  // for disk components, calculate polar r
  //r_disk = 0.220970869121;
  //r_disk = 6.85009694274;
  r_disk = 13.9211647546;
  //r_disk = 20.9922325665;
  // for halo, calculate spherical r
  r_halo = sqrt(x_pos*x_pos + r_disk*r_disk);

  // set properties of halo and disk (these must match initial conditions)
  Real a_disk_z, a_halo, M_vir, M_d, R_vir, R_d, z_d, R_h, M_h, c_vir, phi_0_h, x;
  M_vir = 1.0e12; // viral mass of MW in M_sun
  M_d = 6.5e10; // mass of disk in M_sun
  M_h = M_vir - M_d; // halo mass in M_sun
  R_vir = 261; // viral radius in kpc
  c_vir = 20.0; // halo concentration
  R_h = R_vir / c_vir; // halo scale length in kpc
  R_d = 3.5; // disk scale length in kpc
  z_d = 3.5/5.0; // disk scale height in kpc
  phi_0_h = GN * M_h / (log(1.0+c_vir) - c_vir / (1.0+c_vir));
  x = r_halo / R_h;
  
  // calculate acceleration due to NFW halo & Miyamoto-Nagai disk
  a_halo = - phi_0_h * (log(1+x) - x/(1+x)) / (r_halo*r_halo);
  a_disk_z = - GN * M_d * x_pos * (R_d + sqrt(x_pos*x_pos + z_d*z_d)) / ( pow(r_disk*r_disk + pow(R_d + sqrt(x_pos*x_pos + z_d*z_d), 2), 1.5) * sqrt(x_pos*x_pos + z_d*z_d) );

  // total acceleration is the sum of the halo + disk components
  *gx = (x_pos/r_halo)*a_halo + a_disk_z;

  return;

}


__device__ void calc_g_2D(int xid, int yid, int x_off, int y_off, int n_ghost, Real dx, Real dy, Real xbound, Real ybound, Real *gx, Real *gy)
{
  Real x_pos, y_pos, r, phi;
  // use the subgrid offset and global boundaries to calculate absolute positions on the grid
  x_pos = (x_off + xid - n_ghost + 0.5)*dx + xbound;
  y_pos = (y_off + yid - n_ghost + 0.5)*dy + ybound;

  // for Gresho, also need r & phi
  r = sqrt(x_pos*x_pos + y_pos*y_pos);
  phi = atan2(y_pos, x_pos);

/*
  // set acceleration to balance v_phi in Gresho problem
  if (r < 0.2) {
    *gx = -cos(phi)*25.0*r;
    *gy = -sin(phi)*25.0*r;
  }
  else if (r >= 0.2 && r < 0.4) {
    *gx = -cos(phi)*(4.0 - 20.0*r + 25.0*r*r)/r;
    *gy = -sin(phi)*(4.0 - 20.0*r + 25.0*r*r)/r;
  }
  else {
    *gx = 0.0;
    *gy = 0.0;
  }
*/
/*
  // set gravitational acceleration for Keplarian potential
  Real M;
  M = 1*Msun;
  *gx = -cos(phi)*GN*M/(r*r);
  *gy = -sin(phi)*GN*M/(r*r);
*/
  // set gravitational acceleration for Kuzmin disk + NFW halo
  Real a_d, a_h, a, M_vir, M_d, R_vir, R_d, R_s, M_h, c_vir, x;
  M_vir = 1.0e12; // viral mass of MW in M_sun
  M_d = 6.5e10; // mass of disk in M_sun (assume all gas)
  M_h = M_vir - M_d; // halo mass in M_sun
  R_vir = 261; // viral radius in kpc
  c_vir = 20; // halo concentration
  R_s = R_vir / c_vir; // halo scale length in kpc
  R_d = 3.5; // disk scale length in kpc
  
  // calculate acceleration
  x = r / R_s;
  a_d = GN * M_d * r * pow(r*r + R_d*R_d, -1.5);
  a_h = GN * M_h * (log(1+x)- x / (1+x)) / ((log(1+c_vir) - c_vir / (1+c_vir)) * r*r);
  a = a_d + a_h;

  *gx = -cos(phi)*a;
  *gy = -sin(phi)*a;

  return;
}


__device__ void calc_g_3D_CUDA(int xid, int yid, int zid, int x_off, int y_off, int z_off, int n_ghost, Real dx, Real dy, Real dz, Real xbound, Real ybound, Real zbound, Real *gx, Real *gy, Real *gz)
{
  Real x_pos, y_pos, z_pos, r_disk, r_halo;
  // use the subgrid offset and global boundaries to calculate absolute positions on the grid
  x_pos = (x_off + xid - n_ghost + 0.5)*dx + xbound;
  y_pos = (y_off + yid - n_ghost + 0.5)*dy + ybound;
  z_pos = (z_off + zid - n_ghost + 0.5)*dz + zbound;

  // for disk components, calculate polar r
  r_disk = sqrt(x_pos*x_pos + y_pos*y_pos);
  // for halo, calculate spherical r
  r_halo = sqrt(x_pos*x_pos + y_pos*y_pos + z_pos*z_pos);

  // set properties of halo and disk (these must match initial conditions)
  Real a_disk_r, a_disk_z, a_halo, a_halo_r, a_halo_z;
  Real M_vir, M_d, R_vir, R_d, z_d, R_h, M_h, c_vir, phi_0_h, x;
  // MW model
  M_vir = 1.0e12; // viral mass of in M_sun
  M_d = 6.5e10; // viral mass of in M_sun
  R_d = 3.5; // disk scale length in kpc
  z_d = 3.5/5.0; // disk scale height in kpc
  R_vir = 261.; // virial radius in kpc
  c_vir = 20.0; // halo concentration
  // M82 model
  //M_vir = 5.0e10; // viral mass of in M_sun
  //M_d = 1.0e10; // mass of disk in M_sun
  //R_d = 0.8; // disk scale length in kpc
  //z_d = 0.15; // disk scale height in kpc
  //R_vir = R_d/0.015; // viral radius in kpc
  //c_vir = 10.0; // halo concentration

  M_h = M_vir - M_d; // halo mass in M_sun
  R_h = R_vir / c_vir; // halo scale length in kpc
  phi_0_h = GN * M_h / (log(1.0+c_vir) - c_vir / (1.0+c_vir));
  x = r_halo / R_h;
  
  // calculate acceleration due to NFW halo & Miyamoto-Nagai disk
  a_halo = - phi_0_h * (log(1+x) - x/(1+x)) / (r_halo*r_halo);
  a_halo_r = a_halo*(r_disk/r_halo);
  a_halo_z = a_halo*(z_pos/r_halo);
  a_disk_r = - GN * M_d * r_disk * pow(r_disk*r_disk+ pow(R_d + sqrt(z_pos*z_pos + z_d*z_d),2), -1.5);
  a_disk_z = - GN * M_d * z_pos * (R_d + sqrt(z_pos*z_pos + z_d*z_d)) / ( pow(r_disk*r_disk + pow(R_d + sqrt(z_pos*z_pos + z_d*z_d), 2), 1.5) * sqrt(z_pos*z_pos + z_d*z_d) );

  // total acceleration is the sum of the halo + disk components
  //*gx = 0.0;
  //*gy = 0.0;
  //*gz = 0.0;
  *gx = (x_pos/r_disk)*(a_disk_r+a_halo_r);
  *gy = (y_pos/r_disk)*(a_disk_r+a_halo_r);
  //*gx = x_pos*a_halo/r_halo + x_pos*a_disk/r_disk;
  //*gy = y_pos*a_halo/r_halo + y_pos*a_disk/r_disk;
  *gz = a_disk_z+a_halo_z;
  //*gz = (z_pos/r_halo)*a_halo;
  //*gz = (z_pos/r_halo)*a_halo + a_disk_z;

  return;
}



#endif //CUDA
