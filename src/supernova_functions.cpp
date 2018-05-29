/*! \file supernova_functions.cpp
 *  \brief Definitions of functions used to add SN energy to cells.
           Functions are members of the Grid3D class. */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <cmath>
#include <string.h>
#include <time.h>
#include "global.h"
#include "grid3D.h"
#include "mpi_routines.h"
#include "io.h"
#include "error_handling.h"
#include "ran.h"


//Add a single supernova with 10^51 ergs of thermal energy and 10 M_sun
Real Grid3D::Add_Supernova(void)
{
  int i, j, k, id;
  Real x_pos, y_pos, z_pos, r, R_s;
  Real M, E, V, rho, Ed;
  Real xl, xr, yl, yr, zl, zr, rl, rr;
  int incount, ii;
  Real weight, xpoint, ypoint, zpoint;
  //R_s = 2*H.dx; // supernova radius, pc
  R_s = 5; // supernova radius, pc
  M = 15.0; // mass input, in M_sun
  E = 1.0e51; // energy input, in erg
  Real x_sn, y_sn, z_sn; // central location of SN
  x_sn = y_sn = z_sn = 0;
  Real R_cl = 2.0; // size of cluster, pc (for random distribution of SN)
  srand (int(H.t)); // change location of SN
  x_sn = 2*R_cl*float(rand() % 10000)/10000.0 - R_cl; // pick a random z between -R_cl and R_cl
  y_sn = 2*R_cl*float(rand() % 10000)/10000.0 - R_cl; // pick a random z between -R_cl and R_cl
  z_sn = 2*R_cl*float(rand() % 10000)/10000.0 - R_cl; // pick a random z between -R_cl and R_cl
  printf("%f %f %f\n", x_sn, y_sn, z_sn);

  Real max_dti = 0;

  E = E/(MASS_UNIT*VELOCITY_UNIT*VELOCITY_UNIT); // convert to code units
  V = (4.0/3.0)*PI*R_s*R_s*R_s;
  //chprintf("%f %f %f %f\n", V, (4.0/3.0)*PI*0.3*0.3*0.3, M_dot, E_dot);
  rho = M / V;
  Ed = E / V;

  Real d_inv, vx, vy, vz, P, cs;
  Real max_vx, max_vy, max_vz;
  max_dti = max_vx = max_vy = max_vz = 0.0;
  Real M_dot_tot, E_dot_tot;
  M_dot_tot = E_dot_tot = 0.0;

  for (k=H.n_ghost; k<H.nz-H.n_ghost; k++) {
    for (j=H.n_ghost; j<H.ny-H.n_ghost; j++) {
      for (i=H.n_ghost; i<H.nx-H.n_ghost; i++) {

        id = i + j*H.nx + k*H.nx*H.ny;

        Get_Position(i, j, k, &x_pos, &y_pos, &z_pos);
        
        // calculate spherical radius
        xl = fabs(x_pos-x_sn)-0.5*H.dx;
        yl = fabs(y_pos-y_sn)-0.5*H.dy;
        zl = fabs(z_pos-z_sn)-0.5*H.dz;
        xr = fabs(x_pos-x_sn)+0.5*H.dx;
        yr = fabs(y_pos-y_sn)+0.5*H.dy;
        zr = fabs(z_pos-z_sn)+0.5*H.dz;
        rl = sqrt(xl*xl + yl*yl + zl*zl);
        rr = sqrt(xr*xr + yr*yr + zr*zr);
        //r = sqrt(x_pos*x_pos + y_pos*y_pos + z_pos*z_pos);

        // within SN radius, inject mass and thermal energy
        // entire cell is within sphere
        if (rr < R_s) {
          C.density[id] += rho;
          C.Energy[id] += Ed;
          #ifdef DE
          C.GasEnergy[id] += Ed;
          //Real n = C.density[id]*DENSITY_UNIT/(0.6*MP);
          //Real T = C.GasEnergy[id]*(gama-1.0)*PRESSURE_UNIT/(n*KB);
          //printf("%3d %3d %3d Starburst zone n: %e T:%e.\n", i, j, k, n, T);
          #endif
          #ifdef SCALAR
          C.scalar[id] += 1.0*rho;
          #endif
          //M_dot_tot += rho_dot*H.dx*H.dy*H.dz;
          //E_dot_tot += Ed_dot*H.dx*H.dy*H.dz;
          // recalculate the timestep for these cells
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
        // on the sphere
        //if (rl < R_s && rr > R_s && fabs(z_pos) < z_s) {
        if (rl < R_s && rr > R_s) {
          // quick Monte Carlo to determine weighting
          Ran quickran(50);
          incount = 0;
          for (ii=0; ii<1000; ii++) {
            // generate a random number between x_pos and dx
            xpoint = xl + H.dx*quickran.doub();
            // generate a random number between y_pos and dy
            ypoint = yl + H.dy*quickran.doub();
            // generate a random number between z_pos and dz
            zpoint = zl + H.dz*quickran.doub();
            // check to see whether the point is within the sphere 
            if (xpoint*xpoint + ypoint*ypoint + zpoint*zpoint < R_s*R_s) incount++;
            //if (xpoint*xpoint + ypoint*ypoint < R_s*R_s) incount++;
          }
          weight = incount / 1000.0;
          C.density[id] += rho * weight;
          C.Energy[id]  += Ed * weight;
          #ifdef DE
          C.GasEnergy[id] += Ed * weight;
          #endif
          #ifdef SCALAR
          C.scalar[id] += 1.0*rho*weight;
          #endif
          // recalculate the timestep for these cells
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
  }

  /*
  printf("procID: %d M_dot: %e E_dot: %e\n", procID, M_dot_tot, E_dot_tot);
  MPI_Barrier(MPI_COMM_WORLD);
  Real global_M_dot, global_E_dot;
  MPI_Reduce(&M_dot_tot, &global_M_dot, 1, MPI_CHREAL, MPI_SUM, 0, MPI_COMM_WORLD);
  MPI_Reduce(&E_dot_tot, &global_E_dot, 1, MPI_CHREAL, MPI_SUM, 0, MPI_COMM_WORLD);
  chprintf("Total M_dot: %e E_dot: %e \n", global_M_dot, global_E_dot); 
  */

  // compute max inverse of dt
  max_dti = max_vx / H.dx;
  max_dti = fmax(max_dti, max_vy / H.dy);
  max_dti = fmax(max_dti, max_vz / H.dy);

  return max_dti;

}


Real Grid3D::Add_Supernovae(void)
{
  int i, j, k, id;
  Real x_pos, y_pos, z_pos, r, R_s, z_s, R_c, f, t, t1, t2, t3;
  Real M1, M2, M3, E1, E2, E3, M_dot, E_dot, V, rho_dot, Ed_dot;
  Real r_sn, phi_sn, x_sn, y_sn, z_sn;
  Real xl, xr, yl, yr, zl, zr, rl, rr;
  int incount, ii;
  Real weight, xpoint, ypoint, zpoint;
  Real SFR, E_sn;
  int N_sn, nx_sn, ny_sn, nz_sn;
  int N_cluster;
  SFR = 20.0; // star formation rate, in M_sun / yr
  R_s = 0.75; // starburst radius, in kpc
  R_c = 0.15; // cluster radius, in kpc
  z_s = 0.105; // starburst height, in kpc
  M1 = 1.5e3; // mass input rate, in M_sun / kyr
  E1 = 1.5e42; // energy input rate, in erg/s
  M2 = 12.0e3;
  E2 = 5.4e42;
  M_dot = 0.0;
  E_dot = 0.0;

  Real max_dti = 0;

  N_cluster = 8;

  // start feedback after 5 Myr, ramp up for 5 Myr, high for 30 Myr, ramp down for 5 Myr
  Real tstart = 5.0;
  Real tramp = 5.0;
  Real thigh = 30.0;
  t = H.t/1000 - tstart;
  t1 = tramp;
  t2 = tramp+thigh;
  t3 = 2*tramp+thigh;
  if (t >= 0) {

  if (t >= 0 && t < t1) {
    M_dot = M1 + (1.0/tramp)*t*(M2-M1); 
    E_dot = E1 + (1.0/tramp)*t*(E2-E1);
  }
  if (t >= t1 && t < t2) {
    M_dot = M2;
    E_dot = E2;
  } 
  if (t >= t2 && t < t3) {
    M_dot = M1 + (1.0/tramp)*(t-t3)*(M1-M2);
    E_dot = E1 + (1.0/tramp)*(t-t3)*(E1-E2);
  }
  if (t >= t3) {
    M_dot = M1;
    E_dot = E1;
  }

  E_dot = E_dot*TIME_UNIT/(MASS_UNIT*VELOCITY_UNIT*VELOCITY_UNIT); // convert to code units
  V = N_cluster*(4.0/3.0)*PI*R_c*R_c*R_c;
  //chprintf("%f %f %f %f\n", V, (4.0/3.0)*PI*0.3*0.3*0.3, M_dot, E_dot);
  f = H.dx*H.dy*H.dz / V;
  rho_dot = f * M_dot / (H.dx*H.dy*H.dz);
  Ed_dot = f * E_dot / (H.dx*H.dy*H.dz);

  Real d_inv, vx, vy, vz, P, cs;
  Real max_vx, max_vy, max_vz;
  max_dti = max_vx = max_vy = max_vz = 0.0;
  Real M_dot_tot, E_dot_tot;
  M_dot_tot = E_dot_tot = 0.0;

  //if (H.n_step==0) srand (1); // initialize random seed
  srand (int(t)/15+1); // change location of clusters every 15 Myr 

  for (int nn=0; nn<N_cluster; nn++) {

    r_sn = R_s*float(rand() % 10000)/10000.0; // pick a random radius within R_s
    phi_sn = 2*PI*float(rand() % 10000)/10000.0; // pick a random phi between 0 and 2pi
    z_sn = 2*z_s*float(rand() % 10000)/10000.0 - z_s; // pick a random z between -z_s and z_s
    x_sn = r_sn*cos(phi_sn);
    y_sn = r_sn*sin(phi_sn);

    int xid_sn, yid_sn, zid_sn, nl_x, nl_y, nl_z;
    // identify the global id of the cell containing the cluster center
    xid_sn = int((x_sn + 0.5*H.xdglobal)/H.dx);
    yid_sn = int((y_sn + 0.5*H.ydglobal)/H.dx);
    zid_sn = int((z_sn + 0.5*H.zdglobal)/H.dx);
    // how many cells to loop through around the center
    nl_x = ceil(R_c/H.dx);
    nl_y = ceil(R_c/H.dy);
    nl_z = ceil(R_c/H.dz);
    //chprintf("x: %f y: %f z: %f xid: %d yid: %d zid: %d nx: %d ny: %d nz: %d\n", x_sn, y_sn, z_sn, xid_sn, yid_sn, zid_sn, nl_x, nl_y, nl_z);

    for (int kk=zid_sn-nl_z; kk<=zid_sn+nl_z; kk++) {
    for (int jj=yid_sn-nl_y; jj<=yid_sn+nl_y; jj++) {
    for (int ii=xid_sn-nl_x; ii<=xid_sn+nl_x; ii++) {

      // is this cell in your domain?
      #ifdef MPI_CHOLLA
      if (ii >= nx_local_start && ii < nx_local_start+nx_local && jj >= ny_local_start && jj < ny_local_start+ny_local && kk >= nz_local_start && kk < nz_local_start+nz_local) 
      {
      #endif
        i = ii + H.n_ghost;
        j = jj + H.n_ghost;
        k = kk + H.n_ghost;
        #ifdef MPI_CHOLLA
        i -= nx_local_start;
        j -= ny_local_start;
        k -= nz_local_start;
        #endif

        //printf("procID: %d  ig: %d  jg: %d  kg: %d  il: %d  jl: %d  kl: %d\n", procID, ii, jj, kk, i, j, k);

        // local domain cell id
        id = i + j*H.nx + k*H.nx*H.ny;
        // global position
        Get_Position(i, j, k, &x_pos, &y_pos, &z_pos);
        
        // calculate radius from the cluster center
        xl = fabs(x_pos-x_sn)-0.5*H.dx;
        yl = fabs(y_pos-y_sn)-0.5*H.dy;
        zl = fabs(z_pos-z_sn)-0.5*H.dz;
        xr = fabs(x_pos-x_sn)+0.5*H.dx;
        yr = fabs(y_pos-y_sn)+0.5*H.dy;
        zr = fabs(z_pos-z_sn)+0.5*H.dz;
        rl = sqrt(xl*xl + yl*yl + zl*zl);
        rr = sqrt(xr*xr + yr*yr + zr*zr);
        r = sqrt((x_pos-x_sn)*(x_pos-x_sn) + (y_pos-y_sn)*(y_pos-y_sn) + (z_pos-z_sn)*(z_pos-z_sn));

        // within cluster radius, inject mass and thermal energy
        // entire cell is within sphere
        if (rr < R_c) {
        //if (r < R_c) {
          C.density[id] += rho_dot * H.dt;
          C.Energy[id] += Ed_dot * H.dt;
          #ifdef DE
          C.GasEnergy[id] += Ed_dot * H.dt;
          //Real n = C.density[id]*DENSITY_UNIT/(0.6*MP);
          //Real T = C.GasEnergy[id]*(gama-1.0)*PRESSURE_UNIT/(n*KB);
          //printf("%f %f %f Starburst zone n: %e T:%e.\n", x_pos, y_pos, z_pos, n, T);
          #endif
          //M_dot_tot += rho_dot*H.dx*H.dy*H.dz;
          //E_dot_tot += Ed_dot*H.dx*H.dy*H.dz;
        }
        // on the sphere
        if (rl < R_c && rr > R_c) {
          // quick Monte Carlo to determine weighting
          Ran quickran(50);
          incount = 0;
          for (int mm=0; mm<1000; mm++) {
            // generate a random number between x_pos and dx
            xpoint = xl + H.dx*quickran.doub();
            // generate a random number between y_pos and dy
            ypoint = yl + H.dy*quickran.doub();
            // generate a random number between z_pos and dz
            zpoint = zl + H.dz*quickran.doub();
            // check to see whether the point is within the sphere 
            if (xpoint*xpoint + ypoint*ypoint + zpoint*zpoint < R_c*R_c) incount++;
          }
          weight = incount / 1000.0;
          C.density[id] += rho_dot * H.dt * weight;
          C.Energy[id]  += Ed_dot * H.dt * weight;
          #ifdef DE
          C.GasEnergy[id] += Ed_dot * H.dt * weight;
          //Real n = C.density[id]*DENSITY_UNIT/(0.6*MP);
          //Real T = C.GasEnergy[id]*(gama-1.0)*PRESSURE_UNIT/(n*KB);
          //printf("%f %f %f Starburst zone n: %e T:%e.\n", x_pos, y_pos, z_pos, n, T);
          #endif
          //M_dot_tot += rho_dot*weight*H.dx*H.dy*H.dz;
          //E_dot_tot += Ed_dot*weight*H.dx*H.dy*H.dz;
        }
        // recalculate the timestep for these cells
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
      #ifdef MPI_CHOLLA
      }
      #endif
    }
    }
    }

  }

  /*
  printf("procID: %d M_dot: %e E_dot: %e\n", procID, M_dot_tot, E_dot_tot);
  MPI_Barrier(MPI_COMM_WORLD);
  Real global_M_dot, global_E_dot;
  MPI_Reduce(&M_dot_tot, &global_M_dot, 1, MPI_CHREAL, MPI_SUM, 0, MPI_COMM_WORLD);
  MPI_Reduce(&E_dot_tot, &global_E_dot, 1, MPI_CHREAL, MPI_SUM, 0, MPI_COMM_WORLD);
  chprintf("Total M_dot: %e E_dot: %e \n", global_M_dot, global_E_dot); 
  */

  // compute max inverse of dt
  max_dti = max_vx / H.dx;
  max_dti = fmax(max_dti, max_vy / H.dy);
  max_dti = fmax(max_dti, max_vz / H.dy);

  }

  return max_dti;

}


Real Grid3D::Add_Supernovae_CC85(void)
{
  int i, j, k, id;
  Real x_pos, y_pos, z_pos, r, R_s, z_s, f, t, t1, t2, t3;
  Real M1, M2, M3, E1, E2, E3, M_dot, E_dot, V, rho_dot, Ed_dot;
  Real d_inv, vx, vy, vz, P, cs;
  Real xl, yl, zl, xr, yr, zr, rl, rr;
  int incount, ii;
  Real weight, xpoint, ypoint, zpoint;
  Real max_vx, max_vy, max_vz, max_dti;
  max_dti = max_vx = max_vy = max_vz = 0.0;
  R_s = 0.3; // starburst radius, in kpc
  // High res adiabatic params
  M1 = 1.5e3; 
  E1 = 1.5e42;
  M2 = 12.0e3;
  E2 = 5.4e42;
  M_dot = 0.0;
  E_dot = 0.0;

  // start feedback after 5 Myr, ramp up for 5 Myr, high for 30 Myr, ramp down for 5 Myr
  Real tstart = 5.0;
  Real tramp = 5.0;
  Real thigh = 30.0;
  t = H.t/1000 - tstart;
  t1 = tramp;
  t2 = tramp+thigh;
  t3 = 2*tramp+thigh;
  if (t >= 0) {

  if (t >= 0 && t < t1) {
    M_dot = M1 + (1.0/tramp)*t*(M2-M1); 
    E_dot = E1 + (1.0/tramp)*t*(E2-E1);
  }
  if (t >= t1 && t < t2) {
    M_dot = M2;
    E_dot = E2;
  } 
  if (t >= t2 && t < t3) {
    M_dot = M1 + (1.0/tramp)*(t-t3)*(M1-M2);
    E_dot = E1 + (1.0/tramp)*(t-t3)*(E1-E2);
  }
  if (t >= t3) {
    M_dot = M1;
    E_dot = E1;
  }

  //M_dot = M2;
  //E_dot = E2;

  E_dot = E_dot*TIME_UNIT/(MASS_UNIT*VELOCITY_UNIT*VELOCITY_UNIT); // convert to code units
  V = (4.0/3.0)*PI*R_s*R_s*R_s;
  //V = PI*R_s*R_s*2*z_s;
  f = H.dx*H.dy*H.dz / V;
  rho_dot = f * M_dot / (H.dx*H.dy*H.dz);
  Ed_dot = f * E_dot / (H.dx*H.dy*H.dz);
  //printf("rho_dot: %f  Ed_dot: %f\n", rho_dot, Ed_dot);

  //Real M_dot_tot, E_dot_tot;

  for (k=H.n_ghost; k<H.nz-H.n_ghost; k++) {
    for (j=H.n_ghost; j<H.ny-H.n_ghost; j++) {
      for (i=H.n_ghost; i<H.nx-H.n_ghost; i++) {

        id = i + j*H.nx + k*H.nx*H.ny;

        Get_Position(i, j, k, &x_pos, &y_pos, &z_pos);
        
        // calculate spherical radius
        xl = fabs(x_pos)-0.5*H.dx;
        yl = fabs(y_pos)-0.5*H.dy;
        zl = fabs(z_pos)-0.5*H.dz;
        xr = fabs(x_pos)+0.5*H.dx;
        yr = fabs(y_pos)+0.5*H.dy;
        zr = fabs(z_pos)+0.5*H.dz;
        rl = sqrt(xl*xl + yl*yl + zl*zl);
        rr = sqrt(xr*xr + yr*yr + zr*zr);
        r = sqrt(x_pos*x_pos + y_pos*y_pos + z_pos*z_pos);
        //rl = sqrt(xl*xl + yl*yl);
        //rr = sqrt(xr*xr + yr*yr);
        //r = sqrt(x_pos*x_pos + y_pos*y_pos);

        // within starburst radius, inject mass and thermal energy
        // entire cell is within sphere
        if (rr < R_s) {
        //if (r < R_s) {
          C.density[id] += rho_dot * H.dt;
          C.Energy[id] += Ed_dot * H.dt;
          #ifdef DE
          C.GasEnergy[id] += Ed_dot * H.dt;
          //Real n = C.density[id]*DENSITY_UNIT/(0.6*MP);
          //Real T = C.GasEnergy[id]*(gama-1.0)*PRESSURE_UNIT/(n*KB);
          //printf("%3d %3d %3d Starburst zone n: %e T:%e.\n", i, j, k, n, T);
          #endif
          //M_dot_tot += rho_dot*H.dx*H.dy*H.dz;
          //E_dot_tot += Ed_dot*H.dx*H.dy*H.dz;
          // recalculate the timestep for these cells
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
        // on the sphere
        if (rl < R_s && rr > R_s) {
          // quick Monte Carlo to determine weighting
          Ran quickran(50);
          incount = 0;
          for (ii=0; ii<1000; ii++) {
            // generate a random number between x_pos and dx
            xpoint = xl + H.dx*quickran.doub();
            // generate a random number between y_pos and dy
            ypoint = yl + H.dy*quickran.doub();
            // generate a random number between z_pos and dz
            zpoint = zl + H.dz*quickran.doub();
            // check to see whether the point is within the sphere 
            if (xpoint*xpoint + ypoint*ypoint + zpoint*zpoint < R_s*R_s) incount++;
            //if (xpoint*xpoint + ypoint*ypoint < R_s*R_s) incount++;
          }
          weight = incount / 1000.0;
          C.density[id] += rho_dot * H.dt * weight;
          C.Energy[id]  += Ed_dot * H.dt * weight;
          #ifdef DE
          C.GasEnergy[id] += Ed_dot * H.dt * weight;
          #endif
          // recalculate the timestep for these cells
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
  }

  //printf("%d %e %e\n", procID, M_dot_tot, E_dot_tot);
  //Real global_M_dot, global_E_dot;
  //MPI_Reduce(&M_dot_tot, &global_M_dot, 1, MPI_CHREAL, MPI_SUM, 0, MPI_COMM_WORLD);
  //MPI_Reduce(&E_dot_tot, &global_E_dot, 1, MPI_CHREAL, MPI_SUM, 0, MPI_COMM_WORLD);
  //chprintf("%e %e \n", global_M_dot, global_E_dot*MASS_UNIT*VELOCITY_UNIT*VELOCITY_UNIT/TIME_UNIT); 
  // compute max inverse of dt
  max_dti = max_vx / H.dx;
  max_dti = fmax(max_dti, max_vy / H.dy);
  max_dti = fmax(max_dti, max_vz / H.dy);

  } // t > 0

  return max_dti;

}


void Grid3D::Analysis_Functions(Real *bubble_radius, Real *bubble_mass, Real *bubble_energy, Real *bubble_energy_th) {

  *bubble_radius = 0.0;
  *bubble_mass = 0.0;
  *bubble_energy = 0.0;
  *bubble_energy_th = 0.0;

  Real d, mx, my, mz, P, E, gE;
  Real vx, vy, vz, v, n, T, mu, V;
  Real v_to_kmps = LENGTH_UNIT/TIME_UNIT/1e5;
  Real e_s = MASS_UNIT*LENGTH_UNIT*LENGTH_UNIT / (TIME_UNIT*TIME_UNIT) / 1e51;
  int n_bub = 0;
  mu = 1.27;

  for (int k=H.n_ghost; k<H.nz-H.n_ghost; k++) {
    for (int j=H.n_ghost; j<H.ny-H.n_ghost; j++) {
      for (int i=H.n_ghost; i<H.nx-H.n_ghost; i++) {

        int id = i + j*H.nx + k*H.nx*H.ny;

        d = C.density[id];
        E = C.Energy[id];
        mx = C.momentum_x[id];
        my = C.momentum_y[id];
        mz = C.momentum_z[id];
        E = C.Energy[id];
        P = (E - (0.5/d)*(mx*mx+ my*my+ mz*mz))*(gama-1.0);
        vx = mx/d;
        vy = my/d;
        vz = mz/d;
        v = sqrt(vx*vx + vy*vy + vz*vz);
        n = d*DENSITY_UNIT/(mu*MP);
        #ifdef DE
        gE = C.GasEnergy[id];
        T = C.GasEnergy[id]*(gama-1.0)*PRESSURE_UNIT/(n*KB); 
        #else
        gE = P/(gama-1.0);
        T = P*PRESSURE_UNIT/(n*KB);
        #endif

        // cell is counted as being in the bubble
        if (T > 1e5 || v*v_to_kmps > 1.5) {
          n_bub++;
          *bubble_energy += E;
          *bubble_energy_th += gE;

        }
        if (T > 1e5) {
          *bubble_mass += d*H.dx*H.dy*H.dz;
        }

      }
    }
  }
  V = n_bub*H.dx*H.dy*H.dz;
  *bubble_radius = pow(0.75*V/PI, (1./3.));
  *bubble_energy = *bubble_energy * V * e_s;
  *bubble_energy_th = *bubble_energy_th * V * e_s;

}
