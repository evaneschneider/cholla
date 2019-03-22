/*! \file main.cpp
 *  \brief Program to run the grid code. */

#ifdef MPI_CHOLLA
#include <mpi.h>
#include "mpi_routines.h"
#endif
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "global.h"
#include "grid3D.h"
#include "io.h"
#include "error_handling.h"

#define OUTPUT
//#define CPU_TIME

int main(int argc, char *argv[])
{
  // timing variables
  double start_total, stop_total, start_step, stop_step;
  #ifdef CPU_TIME
  double stop_init, init_min, init_max, init_avg;
  double start_bound, stop_bound, bound_min, bound_max, bound_avg;
  double start_hydro, stop_hydro, hydro_min, hydro_max, hydro_avg;
  double init, bound, hydro;
  init = bound = hydro = 0;
  #endif //CPU_TIME

  // start the total time
  start_total = get_time();

  /* Initialize MPI communication */
  #ifdef MPI_CHOLLA
  InitializeChollaMPI(&argc, &argv);
  #endif /*MPI_CHOLLA*/

  Real dti = 0; // inverse time step, 1.0 / dt

  // input parameter variables
  char *param_file;
  struct parameters P;
  int nfile = 0; // number of output files
  Real outtime = 0; // current output time


  // read in command line arguments
  if (argc != 2)
  {
    chprintf("usage: %s <parameter_file>\n", argv[0]);
    chexit(-1);
  } else {
    param_file = argv[1];
  }

  // create the grid
  Grid3D G;

  // read in the parameters
  parse_params (param_file, &P);
  // and output to screen
  chprintf ("Parameter values:  nx = %d, ny = %d, nz = %d, tout = %f, init = %s, boundaries = %d %d %d %d %d %d\n", 
    P.nx, P.ny, P.nz, P.tout, P.init, P.xl_bcnd, P.xu_bcnd, P.yl_bcnd, P.yu_bcnd, P.zl_bcnd, P.zu_bcnd);
  chprintf ("Output directory:  %s\n", P.outdir);


  // initialize the grid
  G.Initialize(&P);
  chprintf("Local number of grid cells: %d %d %d %d\n", G.H.nx_real, G.H.ny_real, G.H.nz_real, G.H.n_cells);


  // Set initial conditions and calculate first dt
  chprintf("Setting initial conditions...\n");
  G.Set_Initial_Conditions(P);
  chprintf("Initial conditions set.\n");
  // set main variables for Read_Grid inital conditions
  if (strcmp(P.init, "Read_Grid") == 0) {
    dti = C_cfl / G.H.dt;
    outtime += G.H.t;
    nfile = P.nfile*P.nfull;
  }
  
  #ifdef CPU_TIME
  G.Timer.Initialize();
  #endif
  
  #ifdef GRAVITY
  G.Initialize_Gravity(&P);
  #endif
  
  #ifdef PARTICLES
  G.Initialize_Particles(&P);
  #endif
  
  #ifdef COSMOLOGY
  G.Initialize_Cosmology(&P);
  #endif
  
  #ifdef COOLING_GRACKLE
  G.Initialize_Grackle(&P);
  #endif

  #ifdef GRAVITY
  G.Compute_Gravitational_Potential( &P);
  #endif

  // set boundary conditions (assign appropriate values to ghost cells)
  chprintf("Setting boundary conditions...\n");
  G.Set_Boundary_Conditions_All(P);
  chprintf("Boundary conditions set.\n");  
  
  #ifdef PARTICLES
  G.Get_Particles_Accelration();
  #endif

  chprintf("Dimensions of each cell: dx = %f dy = %f dz = %f\n", G.H.dx, G.H.dy, G.H.dz);
  chprintf("Ratio of specific heats gamma = %f\n",gama);
  chprintf("Nstep = %d  Timestep = %f  Simulation time = %f\n", G.H.n_step, G.H.dt, G.H.t);


  #ifdef OUTPUT
  if (strcmp(P.init, "Read_Grid") != 0 || G.H.Output_Now ) {
  // write the initial conditions to file
  chprintf("Writing initial conditions to file...\n");
  WriteData(G, P, nfile);
  }
  // add one to the output file count
  nfile++;
  #endif //OUTPUT
  // increment the next output time
  outtime += P.outstep;

  #ifdef CPU_TIME
  stop_init = get_time();
  init = stop_init - start_total;
  #ifdef MPI_CHOLLA
  init_min = ReduceRealMin(init);
  init_max = ReduceRealMax(init);
  init_avg = ReduceRealAvg(init);
  chprintf("Init  min: %9.4f  max: %9.4f  avg: %9.4f\n", init_min, init_max, init_avg);
  #else
  printf("Init %9.4f\n", init);
  #endif //MPI_CHOLLA
  #endif //CPU_TIME

  // Evolve the grid, one timestep at a time
  chprintf("Starting calculations.\n");
  // while (G.H.t < P.tout)
  while (G.H.n_step < 20)
  {
    chprintf("n_step: %d \n", G.H.n_step + 1 );
    // get the start time
    start_step = get_time();
    
    // calculate the timestep
    G.set_dt(dti);

    if (G.H.t + G.H.dt > outtime) 
    {
      G.H.dt = outtime - G.H.t;
    }

    // This is not necessary, G.set_dt finds the global max_dti, dt is the same in all processes by now
    // #ifdef MPI_CHOLLA
    // G.H.dt = ReduceRealMin(G.H.dt);
    // #endif
    
    
    #ifdef PARTICLES
    //Advance the particles KDK( first step )
    G.Advance_Particles( 1 );
    #endif
    
    // Advance the grid by one timestep
    dti = G.Update_Hydro_Grid();
    
    // update the simulation time ( t += dt )
    G.Update_Time();
    
        
    #ifdef GRAVITY
    //Compute Gravitational potential for next step
    G.Compute_Gravitational_Potential( &P);
    #endif

    // add one to the timestep count
    G.H.n_step++;

    // set boundary conditions for next time step 
    G.Set_Boundary_Conditions_All(P);
    
    #ifdef PARTICLES
    //Advance the particles KDK( second step )
    G.Advance_Particles( 2 );
    #endif

    
    #ifdef CPU_TIME
    G.Timer.Print_Times();
    #endif


    // get the time to compute the total timestep
    stop_step = get_time();
    stop_total = get_time();
    G.H.t_wall = stop_total-start_total;
    #ifdef MPI_CHOLLA
    G.H.t_wall = ReduceRealMax(G.H.t_wall);
    #endif 
    chprintf("n_step: %d   sim time: %10.7f   sim timestep: %7.4e  timestep time = %9.3f ms   total time = %9.4f s\n\n", 
      G.H.n_step, G.H.t, G.H.dt, (stop_step-start_step)*1000, G.H.t_wall);

    if (G.H.t == outtime || G.H.Output_Now )
    {
      #ifdef OUTPUT
      /*output the grid data*/
      WriteData(G, P, nfile);
      // add one to the output file count
      nfile++;
      #endif //OUTPUT
      // update to the next output time
      outtime += P.outstep;      
    }
    
    #ifdef CPU_TIME
    G.Timer.n_steps += 1;
    #endif

    #ifdef COSMOLOGY
    if ( G.Cosmo.current_a >= G.Cosmo.scale_outputs[G.Cosmo.n_outputs-1] ) {
      chprintf( "\nReached Last Cosmological Output: Ending Simulation\n");
      break;
    }
    #endif
/*
    // check for failures
    for (int i=G.H.n_ghost; i<G.H.nx-G.H.n_ghost; i++) {
      for (int j=G.H.n_ghost; j<G.H.ny-G.H.n_ghost; j++) {
        for (int k=G.H.n_ghost; k<G.H.nz-G.H.n_ghost; k++) {
          int id = i + j*G.H.nx + k*G.H.nx*G.H.ny;
          if (G.C.density[id] < 0.0 || G.C.density[id] != G.C.density[id]) {
            printf("Failure in cell %d %d %d. Density %e\n", i, j, k, G.C.density[id]);
            #ifdef MPI_CHOLLA
            MPI_Finalize();
            chexit(-1);
            #endif
            exit(0);
          }
        }
      }
    }
*/   

  } /*end loop over timesteps*/
  
  
  #ifdef CPU_TIME
  G.Timer.Get_Average_Times();
  G.Timer.Print_Average_Times( P );
  #endif

  // free the grid
  G.Reset();

  #ifdef MPI_CHOLLA
  MPI_Finalize();
  #endif /*MPI_CHOLLA*/

  return 0;

}
