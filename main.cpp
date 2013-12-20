#include <mpi.h>

#include <iostream>

#include "simulation_defines.hpp"
#include "memory/domain.hpp"
#include "memory/particle.hpp"
#include "memory/vector3D.hpp"
#include "memory/init_domain.hpp"
#include "communicator/mpi_communicator.hpp"
#include "physics/force_gravity.hpp"
#include "physics/initial_particles.hpp"
//#include "simulation_defines_constexpr.hpp"

using namespace std;
using namespace MDSIM;

int
main( int argc, char *argv[] )
{
  communicator::MPI_Communicator& comm =
    communicator::MPI_Communicator::getInstance();

  /// \todo (l8r) Write load balancing routine for domain size in y:
  ///       - weighting function (left, right, min-relative-difference of n-%)
  ///       - swap cell-lines to left or right

  /// \todo capsule physics as template parameter
  /// \todo (optional) capsule particle attributes / sheme, too

  // set up simulation area with cells
  memory::Domain<double> myDomain = init_domain();

  MPI_Barrier( MPI_COMM_WORLD );


  // Initialize Particles
  //physics::init_SunEarth( myDomain );
  //physics::init_Benchmark( myDomain, 10 );

  memory::vector3D<double> origin1( 50.0 * simParams::LJ,  55.0 * simParams::LJ, 0.0 );
  memory::vector3D<double> origin2( 50.0 * simParams::LJ,  95.0 * simParams::LJ, 0.0 );
  memory::vector3D<double> v0( 1.0 * simParams::HyadenSpeed / 50.0,
                               1.0 * simParams::HyadenSpeed / 8.0,
                               0.0 );

  physics::init_StarCluster( myDomain,
                             origin1,
                             v0,
                             15.0 * simParams::LJ,
                             500 );
  physics::init_StarCluster( myDomain,
                             origin2,
                             v0 * (-1.0),
                             15.0 * simParams::LJ,
                             500 );

  typename communicator::MPI_Communicator::handle hSendToTop = comm.getNullHandle();
  typename communicator::MPI_Communicator::handle hSendToBot = comm.getNullHandle();

  int lastPercent = 0;

  // init Force calc
  myDomain.resetForces();
  myDomain.calculateForces<physics::force_gravity<double> >(
    myDomain.AreaBorder | myDomain.AreaCore | myDomain.AreaHardCore );


  // Main Loop
  int i = 0;
  for( double t = 0.0; t < simParams::simTime; t += simParams::dt )
  {
    i++;
    if( i % simParams::output == 0 )
    {
      myDomain.coutParticlePos( i / simParams::output, false );
      //comm.cleanBuffers();
    }
    //myDomain.coutParticleNum();

    myDomain.clearGhostCells();
    myDomain.moveParticles();
    myDomain.mapParticlesToCells();

    // in-Domain: periodic movement
    myDomain.moveInnerDomainPeriodic( myDomain.XPeriodic );

    // out-of-Domain: MPI Swap to next real Domain and Cloning of
    //                border particles to neighbor ghosts
    std::list<memory::Particle<double> > pOutTop =
      myDomain.getArea( myDomain.Top,
                        myDomain.AreaGhost | myDomain.AreaBorder );
    hSendToTop = comm.sendParticles( pOutTop, comm.Top, false );
    pOutTop.clear();

    std::list<memory::Particle<double> > pOutBottom =
      myDomain.getArea( myDomain.Bottom,
                        myDomain.AreaGhost | myDomain.AreaBorder );
    hSendToBot = comm.sendParticles( pOutBottom, comm.Bottom, false );
    pOutBottom.clear();

    // calc forces for HardCore
    myDomain.resetForces();
    myDomain.calculateForces<physics::force_gravity<double> >(
      myDomain.AreaHardCore );

    // receive particles
    std::vector<memory::Particle<double> > pRecv;
    comm.receiveParticles( pRecv, false );

    // add received particles
    myDomain.addParticle( pRecv, true );
    pRecv.clear();

    // in-Domain: build x-ghosts
    myDomain.createInnerDomainGhosts( myDomain.XPeriodic );

    // calc forces for Border and Core
    myDomain.calculateForces<physics::force_gravity<double> >(
      myDomain.AreaBorder | myDomain.AreaCore );

    // clear send buffers and wait for send finish
    comm.finishSend( hSendToTop );
    comm.finishSend( hSendToBot );

    // calc forces (not overlapped with communication)
    //myDomain.resetForces();
    //myDomain.calculateForces<physics::force_gravity<double> >(
    //  myDomain.AreaHardCore | myDomain.AreaBorder | myDomain.AreaCore );


    // Output process each n percent
    if( int( t / simParams::simTime * (100 / simParams::outPercent ) ) > lastPercent )
    {
      lastPercent = int( t / simParams::simTime * ( 100 / simParams::outPercent ) );
      if( comm.getRank() == 0 )
        std::cout << lastPercent * simParams::outPercent << "% ("
                  << t << "/" << simParams::simTime << ")" << std::endl;
      myDomain.coutParticleNum();
    }
  }


  return 0;
}
