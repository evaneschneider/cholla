/*! \file cluster.cpp
 *  \brief Definitions of functions of the cluster class. */
#ifdef CLUSTERS

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <cmath>
#include <string.h>
#include <time.h>
#include "global.h"
#include "cluster.h"

float S99_table[1000][3];
float cluster_list[N_CL][5];

void Cluster::Initialize(void) {

  flag_on = 0;
  // set the initial cluster mass, radial, azimuthal, and vertical positions from the table
  mass = cluster_list[id][0];
  SF_cl = cluster_list[id][1];
  r_pos = cluster_list[id][2];
  phi_pos = cluster_list[id][3];
  z_pos = cluster_list[id][4];

  // calculate the global x and y positions
  x_pos = r_pos*cos(phi_pos);
  y_pos = r_pos*sin(phi_pos);

  time = 0.0;

  R_cl = 0.03;
  V_cl = (4./3.)*PI*R_cl*R_cl*R_cl;
  
}

void Cluster::Switch(Real t) {

  // For now, just set a constant star formation rate of 20 M_sun/yr
  Real SF_total;
  Real SFR = 20000;
  SF_total = SFR*t;
  
  // if the total star formation is less than the desired total,
  // turn the cluster on
  if (SF_cl < SF_total) {
    flag_on = 1;
  }
  // if the cluster has been on more than 40 Myr, turn it off
  if (time > 40000) {
    flag_on = 0;
  }

}

void Cluster::Rotate(Real dt) {

  Real r_sph, a, v, phi;

  // for halo component, calculate spherical r
  r_sph = sqrt(x_pos * x_pos + y_pos * y_pos + z_pos*z_pos);

  // set properties of halo and disk (these must match gravity_cuda)
  Real a_disk_r, a_halo, a_halo_r;
  Real M_vir, M_d, R_vir, R_d, z_d, R_h, M_h, c_vir, phi_0_h, x;
  // MW model
  //M_vir = 1.0e12; // viral mass of in M_sun
  //M_d = 6.5e10; // viral mass of in M_sun
  //R_d = 3.5; // disk scale length in kpc
  //z_d = 3.5/5.0; // disk scale height in kpc
  //R_vir = 261.; // virial radius in kpc
  //c_vir = 20.0; // halo concentration
  // M82 model
  M_vir = 5.0e10; // viral mass of in M_sun
  M_d = 1.0e10; // mass of disk in M_sun
  R_d = 0.8; // disk scale length in kpc
  z_d = 0.15; // disk scale height in kpc
  R_vir = R_d/0.015; // viral radius in kpc
  c_vir = 10.0; // halo concentration

  M_h = M_vir - M_d; // halo mass in M_sun
  R_h = R_vir / c_vir; // halo scale length in kpc
  phi_0_h = GN * M_h / (log(1.0+c_vir) - c_vir / (1.0+c_vir));
  x = r_sph / R_h;
  
  // calculate acceleration due to NFW halo & Miyamoto-Nagai disk
  a_halo = - phi_0_h * (log(1+x) - x/(1+x)) / (r_sph*r_sph);
  a_halo_r = a_halo*(r_pos/r_sph);
  a_disk_r = - GN * M_d * r_pos * pow(r_pos*r_pos+ pow(R_d + sqrt(z_pos*z_pos + z_d*z_d),2), -1.5);
  // total acceleration is the sum of the halo + disk components
  a = fabs(a_halo_r) + fabs(a_disk_r);
  // radial velocity
  v = sqrt(r_pos*a);
  // how far has the cluster gone?
  phi_pos += v*dt/r_pos;

  // set new cluster center location
  x_pos = r_pos*cos(phi_pos);
  y_pos = r_pos*sin(phi_pos);

  // update the time
  time += dt;

}

// returns M_dot and E_dot (code units)
void Cluster::Get_S99_Fluxes(void) {

  int i;
  Real del_t = 1e5;
  Real tf;
  Real M_slope, E_slope;
  Real f;

  // S99 table is calibrated to a 1e6 M_sun cluster
  f = mass/1e6;
  
  // determine where in the table we are
  tf = (time*1e3)/del_t;
  i = int(tf);
  // interpolate between the table points
  M_slope = S99_table[i+1][1]-S99_table[i][1];
  E_slope = S99_table[i+1][2]-S99_table[i][2];
  
  M_dot = f*(S99_table[i][1] + (tf-i)*M_slope); 
  E_dot = f*pow(10, (S99_table[i][2] + (tf-i)*E_slope) )*TIME_UNIT/(MASS_UNIT*VELOCITY_UNIT*VELOCITY_UNIT);

}


void Load_Cluster_List()
{
  int i;

  FILE *infile;
  char buffer[0x1000];
  char * pch;

  // Read in values from cluster list 
  i=0;
  infile = fopen("./cluster_list.txt", "r");
  if (infile == NULL) {
    printf("Unable to open cluster list.\n");
    exit(1);
  }
  while (fgets(buffer, sizeof(buffer), infile) != NULL)
  {
    if (buffer[0] == '#') {
      continue;
    }
    else {
      pch = strtok(buffer, "\t");
      cluster_list[i][0] = atof(pch);
      while (pch != NULL)
      {
        pch = strtok(NULL, "\t");
        if (pch != NULL)
          cluster_list[i][1] = atof(pch);
        pch = strtok(NULL, "\t");
        if (pch != NULL)
          cluster_list[i][2] = atof(pch);
        pch = strtok(NULL, "\t");
        if (pch != NULL)
          cluster_list[i][3] = atof(pch);
        pch = strtok(NULL, "\t");
        if (pch != NULL)
          cluster_list[i][4] = atof(pch);
      }
      i++;
    }
  }
  fclose(infile);

}

void Load_S99_Table(void) {

  int i;
  int nx = 1000;

  FILE *infile;
  char buffer[0x1000];
  char * pch;

  // Read in data from Starburst 99 calculation (done for a 10^6 M_sun cluster) 
  i=0;
  infile = fopen("./S99_table.txt", "r");
  if (infile == NULL) {
    printf("Unable to open Starburst99 file.\n");
    exit(1);
  }
  while (fgets(buffer, sizeof(buffer), infile) != NULL)
  {
    if (buffer[0] == '#') {
      continue;
    }
    else {
      pch = strtok(buffer, "\t");
      S99_table[i][0] = atof(pch);
      while (pch != NULL)
      {
        pch = strtok(NULL, "\t");
        if (pch != NULL)
          S99_table[i][1] = atof(pch);
        pch = strtok(NULL, "\t");
        if (pch != NULL)
          S99_table[i][2] = atof(pch);
      }
      i++;
    }
  }
  fclose(infile);

}

#endif
