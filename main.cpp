#include <mpi.h>

//#include <stdio.h>
#include <iostream>
//#include <stdlib.h>

#include "simulation_defines.hpp"
#include "memory/domain.hpp"
#include "memory/particle.hpp"
#include "memory/vector3D.hpp"
#include "communicator/mpi_communicator.hpp"

using namespace std;
using namespace MDSIM;

int
main( int argc, char *argv[] )
{ 
  communicator::MPI_Communicator& comm =
    communicator::MPI_Communicator::getInstance();
  
  /// \todo Write Communicator for particle positions: copy border to ghost
  /// \todo (l8r) Write load balancing routine for domain size in y:
  ///       - weighting function (left, right, min-relative-difference of n-%)
  ///       - swap cell-lines to left or right

  /// \todo capsule physics as template parameter
  /// \todo (optional) capsule particle attributes / sheme, too

  // set up simulation area with cells
  const int myXCells  = simParams::cellsX;
  const int myXOffset = 0;
        int myYCells  = floor( double(simParams::cellsY) / comm.getSize() );
  const int myYOffset = myYCells * comm.getRank();
  
  if( simParams::cellsY < comm.getSize() )
  {
    std::cout << "Error: CellsY( " << simParams::cellsY
              << " ) < Size( " << comm.getSize() << " )"
              << std::endl;
    return 1;
  }
  
  // expand at last cpu if not Mod 0
  if( comm.getRank() == comm.getSize() -1 && comm.getSize()*myYCells != simParams::cellsY )
  {
    const int expandYatLastRank = simParams::cellsY - comm.getSize()*myYCells;
    myYCells += expandYatLastRank;
    std::cout << "Info: Added " << expandYatLastRank
              << " YCellLines to last Rank because " << simParams::cellsY
              << "%" << comm.getSize() << " != 0 (= cellsY global % size )"
              << std::endl;
  }
  
  std::cout << "Info: MyDomain (" << myXOffset << " - "
            << myXOffset + myXCells -1 << "; "
            << myYOffset << " - " << myYOffset + myYCells -1 << ")"
            << " at Rank " << comm.getRank()
            << std::endl;
  
  memory::Domain<double> myDomain( myXCells, myYCells, myXOffset, myYOffset );
  
  MPI_Barrier( MPI_COMM_WORLD );


  // Initialize Particles
  /// \todo from density and velocity function
  memory::vector3D<double> rE( 0.0, simParams::distance_Sun_Earth, 0.0 );
  memory::vector3D<double> vE( simParams::EarthSpeed, 0.0, 0.0 );
  memory::Particle<double> earth( rE, vE, simParams::mass * simParams::partialEarthSun );

  memory::vector3D<double> rS( 0.0, 0.0, 0.0 );
  memory::vector3D<double> vS( 0.0, 0.0, 0.0 );
  memory::Particle<double> sun( rS, vS, simParams::mass );

  // initialize particles
  myDomain.addParticle( earth );
  //myDomain.addParticle( sun );

  myDomain.coutParticlePos();
  for( double t = 0.0; t < simParams::simTime; t += simParams::dt )
  {
    myDomain.resetForces();
    
    /// \todo also for NxN with next neighbors!
    myDomain.calculateForces();
    myDomain.clearGhostCells();
    
    myDomain.moveParticles();
    myDomain.mapParticlesToCells();
    
    // in-Domain: periodic movement
    std::list<memory::Particle<double>* > pIn =
      myDomain.getArea( myDomain.Left | myDomain.Right,
                        myDomain.AreaGhost );
    myDomain.moveInnerDomainPeriodic( pIn,
                                      myDomain.XPeriodic );
    
    // out-of-Domain: MPI Swap to next real Domain and Cloning of
    //                border particles to neighbor ghosts
    std::list<memory::Particle<double>* > pOutTop =
      myDomain.getArea( myDomain.Top,
                        myDomain.AreaBorder | myDomain.AreaGhost );
    /// \todo send particle to top neighbor (pos, vel, mass)
    std::list<memory::Particle<double>* > pOutBottom =
      myDomain.getArea( myDomain.Bottom,
                        myDomain.AreaBorder | myDomain.AreaGhost );
    /// \todo send particle to bottom neighbor (pos, vel, mass)
    
    /// \todo receive particles
    
    myDomain.coutParticlePos();
  }

  
  return 0;
}
