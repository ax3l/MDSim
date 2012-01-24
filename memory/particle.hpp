#ifndef PARTICLE_HPP
#define PARTICLE_HPP

#include "vector3D.hpp"

namespace MDSIM
{
  namespace memory
  {

    /// Particle
    ///
    /// A standard particle with its Attributes.
    ///
    /// \tparam floatType Double of Float: for Particle Attributes

    template <typename floatType>
    class Particle
    {
    private:
      // position in cell
      // [0., 1.)
      vector3D<floatType> _position;
      vector3D<floatType> _velocity;
      /// \todo calculate with impulses using SRT
      vector3D<floatType> _force;
      floatType _mass;

    public:

      Particle( const vector3D<floatType> pos,
                const vector3D<floatType> vel,
                const floatType mass ) :
      _position( pos ),
      _velocity( vel ),
      _force( 0. ),
      _mass( mass )
      {
      }

      Particle( const Particle<floatType>& p ) :
      _position( p.getPosition( ) ),
      _velocity( p.getVelocity( ) ),
      _force( p.getForce( ) ),
      _mass( p.getMass( ) )
      {
      }

      ~Particle( )
      {
      }

      inline vector3D<floatType> getPosition( ) const
      {
        return _position;
      }

      inline vector3D<floatType> getVelocity( ) const
      {
        return _velocity;
      }

      inline vector3D<floatType> getForce( ) const
      {
        return _force;
      }

      inline floatType getMass( ) const
      {
        return _mass;
      }
      
      inline void resetForce( )
      {
        _force.x = 0.0;
        _force.y = 0.0;
        _force.z = 0.0;
      }

      inline void addPosition( const vector3D<floatType>& off )
      {
        _position += off;
      }

      inline void addVelocity( const vector3D<floatType>& off )
      {
        _velocity += off;
      }

      inline void addForce( const vector3D<floatType>& off )
      {
        _force += off;
      }
    };

  } // namespace memory
} // namespace MDSIM

#endif // PARTICLE_HPP