#include <mpi.h>

//#include <stdio.h>
#include <iostream>
//#include <stdlib.h>

#include "simulation_defines.hpp"
#include "memory/domain.hpp"
#include "memory/particle.hpp"
#include "memory/vector3D.hpp"
#include "simulation_defines _constexpr.hpp"

using namespace std;
using namespace MDSIM;

int
main( int argc, char *argv[] )
{
  int rank, size;

  MPI_Init( &argc, &argv ); /* starts MPI */
  MPI_Comm_rank( MPI_COMM_WORLD, &rank ); /* get current process id */
  MPI_Comm_size( MPI_COMM_WORLD, &size ); /* get number of processes */

  printf( "Hello world from process %d of %d\n", rank, size );
  
  MPI_Barrier( MPI_COMM_WORLD );
  
  /// \todo Init MPI 1-D Communication Matrix ( Option: periodic)
  ///       in "y" direction above each other
  /// \todo Write Communicator for particle positions: copy border to ghost
  /// \todo (l8r) Write load balancing routine for domain size in y:
  ///       - weighting function (left, right, min-relative-difference of n-%)
  ///       - swap cell-lines to left or right

  // set up simulation area with cells
  /// \todo set option: x-direction periodic (y-direction option by MPI)
  /// \todo capsule physics as template parameter
  /// \todo (optional) capsule particle attributes / sheme, too
  memory::Domain<double> myDomain( simParams::cellsX, simParams::cellsY, 0, 0 );

  memory::vector3D<double> rE( 0.0, simParams::distance_Sun_Earth, 0.0 );
  memory::vector3D<double> vE( simParams::EarthSpeed, 0.0, 0.0 );
  memory::Particle<double> earth( rE, vE, simParams::mass * simParams::partialEarthSun );

  memory::vector3D<double> rS( 0.0, 0.0, 0.0 );
  memory::vector3D<double> vS( 0.0, 0.0, 0.0 );
  memory::Particle<double> sun( rS, vS, simParams::mass );

  // initialize particles
  myDomain.addParticle( earth );
  myDomain.addParticle( sun );

  for( double t = 0.0; t < simParams::simTime; t += simParams::dt )
  {
    myDomain.resetForces();
    myDomain.calculateForces();
    myDomain.moveParticles();
    myDomain.coutParticlePos();
  }

  
  MPI_Finalize();
  return 0;
}
