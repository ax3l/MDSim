#ifndef FORCE_GRAVITY_HPP
#define FORCE_GRAVITY_HPP

//#include "<list>"
#include "../simulation_defines.hpp"
#include "../memory/vector3D.hpp"
#include "../memory/particle.hpp"

namespace MDSIM
{
  namespace physics
  {

    /// Force Model for Gravity
    ///
    /// Functor
    ///
    /// \tparam floatType for float or double precision
    ///
    template <typename floatType>
    class force_gravity
    {
    private:
      
    public:
      /// Operator () which returns the force on Particle1
      ///
      /// \param[in] std::list<memory::Particle<floatType> >::iterator,
      ///                 Iterator to Particle 1
      /// \param[in] std::list<memory::Particle<floatType> >::iterator,
      ///                 Iterator to Particle 2
      /// \param[in] memory::vector3D<floatType> cellOffset
      ///                 CellOffset for Particle 2 in Respect to Particle 1
      /// \return memory::vector3D<floatType> with Force on Particle 1
      ///
      typename memory::vector3D<floatType>
      operator( )( typename std::list<memory::Particle<floatType> >::iterator& p1,
                   typename std::list<memory::Particle<floatType> >::iterator& p2,
                   typename memory::vector3D<floatType>& cellOff ) const
      {
        
        if( p2 == p1 )
          return memory::vector3D<floatType>( 0.0, 0.0, 0.0 );

        // F = G*m1*m2/r^2
        memory::vector3D<floatType> r( p2->getPosition( ) + cellOff - p1->getPosition( ) );
        memory::vector3D<floatType> er( r / sqrt( r.abs2( ) ) );

        const floatType fAbsTmp = simParams::G * p1->getMass( ) * p2->getMass( );
        const floatType fAbs = fAbsTmp / r.abs2( );
        //std::cout << "Abs: " << r.abs2() << std::endl;

        memory::vector3D<floatType> force( fAbs * er.x, fAbs * er.y, 0.0 );
        
        return force;
      }
      
    };

  } // namespace physics
} // namespace MDSIM

#endif // FORCE_GRAVITY_HPP