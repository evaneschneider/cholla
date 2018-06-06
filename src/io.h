#ifndef IO_CHOLLA_H
#define IO_CHOLLA_H

#include"global.h"
#include"grid3D.h"

/* Output the history variables to file. */
void WriteHistory(Grid3D G, struct parameters P);

/* Write the data */
void WriteData(Grid3D G, parameters P);

/* Output the grid data to file. */
void OutputData(Grid3D G, struct parameters P);

/* Output a projection of the grid data to file. */
void OutputProjectedData(Grid3D G, parameters P);

/* Output a rotated projection of the grid data to file. */
void OutputRotatedProjectedData(Grid3D G, parameters P);

/* Output xy, xz, and yz slices of the grid data to file. */
void OutputSlices(Grid3D G, parameters P);

/* MPI-safe printf routine */
int chprintf(const char * __restrict sdata, ...);

#endif /*IO_CHOLLA_H*/
