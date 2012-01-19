//#include <mpi.h>

//#include <stdio.h>
#include <iostream>
//#include <stdlib.h>

#include "simulation_defines.hpp"
#include "memory/domain.hpp"

using namespace std;
using namespace MDSIM;


int main( int argc, char *argv[] )
{
  //int rank, size;

  //MPI_Init( &argc, &argv ); /* starts MPI */
  //MPI_Comm_rank( MPI_COMM_WORLD, &rank ); /* get current process id */
  //MPI_Comm_size( MPI_COMM_WORLD, &size ); /* get number of processes */

  //printf( "Hello world from process %d of %d\n", rank, size );
  
  // set up simulation area with cells
  memory::Domain<double> myDomain( simParams::cellsX, simParams::cellsY, 0, 0);
  
  // initialize particles
  
  
  
  //MPI_Finalize();
  return 0;
}
