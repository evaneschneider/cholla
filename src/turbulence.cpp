/*! \file turbulence.cpp
 *  \brief Definitions of grid functions used to generate turbulence. */

#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include"global.h"
#include"grid3D.h"

/*! \fn Apply_Forcing(void)
 *  \brief Apply a forcing field to continuously generate turbulence. */
void Grid3D::Apply_Forcing(void)
{
  int i, j, k, id, fid;
  int n_cells = H.nx_real * H.ny_real * H.nz_real;  
  int istart, jstart, kstart, iend, jend, kend;
  Real x_pos, y_pos, z_pos;
  Real *A1, *A2, *A3, *A4, *B1, *B2, *B3, *B4;
  Real *vxp, *vyp, *vzp;
  Real *p;
  Real mxp, myp, mzp;
  Real vx_av, vy_av, vz_av, v_av;
  Real vxsq, vysq, vzsq;
  Real mxsq, mysq, mzsq;
  Real mx, my, mz;
  Real M, d_tot;  
  A1 = rng_direction(3);
  A2 = rng_direction(3);
  A3 = rng_direction(3);
  A4 = rng_direction(3);
  B1 = rng_direction(3);
  B2 = rng_direction(3);
  B3 = rng_direction(3);
  B4 = rng_direction(3);
  p = rng_direction(3);
  vxp = (Real *) malloc(n_cells*sizeof(Real));
  vyp = (Real *) malloc(n_cells*sizeof(Real));
  vzp = (Real *) malloc(n_cells*sizeof(Real));
  
  //printf("%f %f %f\n", A1[0], A1[1], A1[2]);
  //printf("%f %f %f\n", B1[0], B1[1], B1[2]);
  //printf("%f %f %f\n", A2[0], A2[1], A2[2]);
  //printf("%f %f %f\n", B2[0], B2[1], B2[2]);
  //set the desired Mach number
  M = 2.0;
  d_tot = 0;
  istart = H.n_ghost;
  iend   = H.nx-H.n_ghost;
  if (H.ny > 1) {
    jstart = H.n_ghost;
    jend   = H.ny-H.n_ghost;
  }
  else {
    jstart = 0;
    jend   = H.ny;
  }
  if (H.nz > 1) {
    kstart = H.n_ghost;
    kend   = H.nz-H.n_ghost;
  }
  else {
    kstart = 0;
    kend   = H.nz;
  }
  mxp = myp = mzp = 0.0;
  mx = my = mz = 0.0;
  vx_av = vy_av = vz_av = 0.0;
  for(k=kstart; k<kend; k++) {
    for(j=jstart; j<jend; j++) {
      for(i=istart; i<iend; i++) {
        //get cell index
        id = i + j*H.nx + k*H.nx*H.ny;
        fid = (i-H.n_ghost) + (j-H.n_ghost)*H.nx_real + (k-H.n_ghost)*H.nx_real*H.ny_real;
        // get cell-centered position
        Get_Position(i, j, k, &x_pos, &y_pos, &z_pos);
        
        // calculate velocity perturbations
        vxp[fid] = B1[0]*sin( (2*PI/sqrt(A1[0]*A1[0] + A1[1]*A1[1] + A1[2]*A1[2])) * (A1[0]*(x_pos+p[0]) + A1[1]*(y_pos+p[1]) + A1[2]*(z_pos+p[2])));
        vyp[fid] = B1[1]*sin( (2*PI/sqrt(A1[0]*A1[0] + A1[1]*A1[1] + A1[2]*A1[2])) * (A1[0]*(x_pos+p[0]) + A1[1]*(y_pos+p[1]) + A1[2]*(z_pos+p[2])));
        vzp[fid] = B1[2]*sin( (2*PI/sqrt(A1[0]*A1[0] + A1[1]*A1[1] + A1[2]*A1[2])) * (A1[0]*(x_pos+p[0]) + A1[1]*(y_pos+p[1]) + A1[2]*(z_pos+p[2])));
        vxp[fid] += B2[0]*sin( (4*PI/sqrt(A2[0]*A2[0] + A2[1]*A2[1] + A2[2]*A2[2])) * (A2[0]*(x_pos+p[0]) + A2[1]*(y_pos+p[1]) + A2[2]*(z_pos+p[2])));
        vyp[fid] += B2[1]*sin( (4*PI/sqrt(A2[0]*A2[0] + A2[1]*A2[1] + A2[2]*A2[2])) * (A2[0]*(x_pos+p[0]) + A2[1]*(y_pos+p[1]) + A2[2]*(z_pos+p[2])));
        vzp[fid] += B2[2]*sin( (4*PI/sqrt(A2[0]*A2[0] + A2[1]*A2[1] + A2[2]*A2[2])) * (A2[0]*(x_pos+p[0]) + A2[1]*(y_pos+p[1]) + A2[2]*(z_pos+p[2])));
        vxp[fid] += B3[0]*cos( (2*PI/sqrt(A3[0]*A3[0] + A3[1]*A3[1] + A3[2]*A3[2])) * (A3[0]*(x_pos+p[0]) + A3[1]*(y_pos+p[1]) + A3[2]*(z_pos+p[2])));
        vyp[fid] += B3[1]*cos( (2*PI/sqrt(A3[0]*A3[0] + A3[1]*A3[1] + A3[2]*A3[2])) * (A3[0]*(x_pos+p[0]) + A3[1]*(y_pos+p[1]) + A3[2]*(z_pos+p[2])));
        vzp[fid] += B3[2]*cos( (2*PI/sqrt(A3[0]*A3[0] + A3[1]*A3[1] + A3[2]*A3[2])) * (A3[0]*(x_pos+p[0]) + A3[1]*(y_pos+p[1]) + A3[2]*(z_pos+p[2])));
        vxp[fid] += B4[0]*cos( (4*PI/sqrt(A4[0]*A4[0] + A4[1]*A4[1] + A4[2]*A4[2])) * (A4[0]*(x_pos+p[0]) + A4[1]*(y_pos+p[1]) + A4[2]*(z_pos+p[2])));
        vyp[fid] += B4[1]*cos( (4*PI/sqrt(A4[0]*A4[0] + A4[1]*A4[1] + A4[2]*A4[2])) * (A4[0]*(x_pos+p[0]) + A4[1]*(y_pos+p[1]) + A4[2]*(z_pos+p[2])));
        vzp[fid] += B4[2]*cos( (4*PI/sqrt(A4[0]*A4[0] + A4[1]*A4[1] + A4[2]*A4[2])) * (A4[0]*(x_pos+p[0]) + A4[1]*(y_pos+p[1]) + A4[2]*(z_pos+p[2])));
        
        // calculate momentum of forcing field
        mxp = C.density[id] * vxp[fid];
        myp = C.density[id] * vyp[fid];
        mzp = C.density[id] * vzp[fid];
        // track total momentum
        mx += mxp;
        my += myp;
        mz += mzp;
        d_tot += C.density[id];
      }
    }
  }
  // calculate density weighted average velocity of forcing field
  vx_av = mx / d_tot;
  vy_av = my / d_tot;
  vz_av = mz / d_tot;
  //printf("%f %f %f\n", vx_av, vy_av, vz_av);
  mxsq = mysq = mzsq = 0.0;
  mx = my = mz = 0.0;
  for(k=kstart; k<kend; k++) {
    for(j=jstart; j<jend; j++) {
      for(i=istart; i<iend; i++) {
        //get cell index
        id = i + j*H.nx + k*H.nx*H.ny;
        fid = (i-H.n_ghost) + (j-H.n_ghost)*H.nx_real + (k-H.n_ghost)*H.nx_real*H.ny_real;
        // get cell-centered position
        Get_Position(i, j, k, &x_pos, &y_pos, &z_pos);
        // subtract off average velocity to create a field with zero net momentum
        vxp[fid] -= vx_av; 
        vyp[fid] -= vy_av; 
        vzp[fid] -= vz_av; 
        // calculate momentum of forcing field
        mxp = C.density[id]*vxp[fid];
        myp = C.density[id]*vyp[fid];
        mzp = C.density[id]*vzp[fid];
        // track total momentum
        mx += C.momentum_x[id] + mxp;
        my += C.momentum_y[id] + myp;
        mz += C.momentum_z[id] + mzp;
        // calculate <v^2> for each direction
        mxsq += mxp*mxp/C.density[id];
        mysq += myp*myp/C.density[id];
        mzsq += mzp*mzp/C.density[id];
      }
    }
  }
  vx_av = sqrt(mxsq / d_tot);
  vy_av = sqrt(mysq / d_tot);
  vz_av = sqrt(mzsq / d_tot);
  v_av = sqrt(vx_av*vx_av + vy_av*vy_av + vz_av*vz_av);
  //printf("%f %f %f %f\n", vx_av, vy_av, vz_av, v_av);
  //printf("%f %f %f\n", mx, my, mz);   
  mx = my = mz = 0.0;
  mxp = myp = mzp = 0.0;
  mxsq = mysq = mzsq = 0.0;
  for(k=kstart; k<kend; k++) {
    for(j=jstart; j<jend; j++) {
      for(i=istart; i<iend; i++) {
        // get cell index
        id = i + j*H.nx + k*H.nx*H.ny;
        fid = (i-H.n_ghost) + (j-H.n_ghost)*H.nx_real + (k-H.n_ghost)*H.nx_real*H.ny_real;
        
        // rescale velocities to get desired Mach number
        // only apply a tenth of initial energy since forcing is 
        // applied every tenth of a crossing time
        vxp[fid] *= sqrt(0.1*M*M/3.0) / vx_av;
        vyp[fid] *= sqrt(0.1*M*M/3.0) / vy_av;
        vzp[fid] *= sqrt(0.1*M*M/3.0) / vz_av;
        // calculate momentum perturbations and apply
        mxp = C.density[id]*vxp[fid];
        myp = C.density[id]*vyp[fid];
        mzp = C.density[id]*vzp[fid];
        C.momentum_x[id] += mxp;
        C.momentum_y[id] += myp;
        C.momentum_z[id] += mzp;
        C.Energy[id] += 0.5*(mxp*mxp + myp*myp + mzp*mzp)/C.density[id];        
        // track total momentum
        mx += C.momentum_x[id];
        my += C.momentum_y[id];
        mz += C.momentum_z[id];
        // calcultate <v^2> for each direction
        mxsq += C.momentum_x[id]*C.momentum_x[id]/C.density[id];
        mysq += C.momentum_y[id]*C.momentum_y[id]/C.density[id];
        mzsq += C.momentum_z[id]*C.momentum_z[id]/C.density[id];        
        
      }
    }
  }
  vx_av = sqrt(mxsq / d_tot);
  vy_av = sqrt(mysq / d_tot);
  vz_av = sqrt(mzsq / d_tot);
  v_av = sqrt(vx_av*vx_av + vy_av*vy_av + vz_av*vz_av);
  printf("%f %f %f %f\n", vx_av, vy_av, vz_av, v_av);
  //printf("%f %f %f\n", mx, my, mz);  
  free(vxp);
  free(vyp);
  free(vzp);
  free(A1);
  free(A2);
  free(A3);
  free(A4);
  free(B1);
  free(B2);
  free(B3);
  free(B4);
  free(p);
}
