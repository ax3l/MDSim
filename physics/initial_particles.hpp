#ifndef INITIAL_PARTICLES_HPP
#define INITIAL_PARTICLES_HPP

#include <ctime>
#include <cstdlib>

#include "../simulation_defines.hpp"
#include "../memory/vector3D.hpp"
#include "../memory/particle.hpp"
#include "../memory/domain.hpp"

namespace MDSIM
{
  namespace physics
  {

    /// Init Sun and Earth
    ///
    template <typename floatType>
    void
    init_SunEarth( memory::Domain<floatType>& myDomain )
    {
      const memory::vector3D<double> rE( 0.0, simParams::distance_Sun_Earth, 0.0 );
      const memory::vector3D<double> vE( simParams::EarthSpeed, 0.0, 0.0 );
      //memory::vector3D<double> vE( 0.0, simParams::EarthSpeed, 0.0 );
      memory::Particle<double> earth( rE, vE, simParams::mass * simParams::partialEarthSun );

      const memory::vector3D<double> rS( 0.0, 0.0, 0.0 );
      const memory::vector3D<double> vS( 0.0, 0.0, 0.0 );
      memory::Particle<double> sun( rS, vS, simParams::mass );

      // initialize particles
      myDomain.addParticle( earth );
      myDomain.addParticle( sun );
    }
    
    /// Init Benchmark
    ///
    template <typename floatType>
    void
    init_Benchmark( memory::Domain<floatType>& myDomain,
                    int particlesPerCell = 10 )
    {
      // start of my Domain
      double x0;
      double y0 = myDomain.getFirstCellPos().y;
      double step = simParams::cutoff;
      // end of my Domain
      double xE = myDomain.getLastCellPos().x;
      double yE = myDomain.getLastCellPos().y;
      
      // init random generator
      time_t seed;
      //seed = time (NULL);
      seed = 958838454; // deterministic init for benchmark
      srand ( seed );
      
      while( y0 < yE )
      {
        x0 = myDomain.getFirstCellPos().x;
        while( x0 < xE )
        {
          for( int i = 0; i < particlesPerCell; i++ )
          {
            double rx = double( rand() ) / double( RAND_MAX );
            double ry = double( rand() ) / double( RAND_MAX );
            
            const memory::vector3D<double> r( x0 + rx, y0 + ry, 0.0 );
            const memory::vector3D<double> v( 0.0, 0.0, 0.0 );
            memory::Particle<double> particle( r, v, simParams::mass );

            myDomain.addParticle( particle );
          }
          x0 += step;
        }
        y0 += step;
      }
    }

  } // namespace physics
} // namespace MDSIM

#endif // INITIAL_PARTICLES_HPP