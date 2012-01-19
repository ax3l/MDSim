#ifndef CELL_HPP
#define CELL_HPP

#include <list>
#include "particle.hpp"

namespace MDSIM {
  namespace memory {
  
    /// Cell
    ///
    /// A standard cell which stores the information about the
    /// particles within it.
    ///
    /// \tparam floatType Double of Float: for Particle Attributes
    template <typename floatType>
    class Cell
    {
    private:
      // Cell Position (Cell Number) on the Node
      int _localX;
      int _localY;
      
      std::list<Particle<floatType> > _particles;
      
    public:
      Cell( const int localX,
            const int localY ):
        _localX( localX ),
        _localY( localY )
      {}
      
      Cell( const Cell<floatType>& c ):
        _localX( c.getX() ),
        _localY( c.getY() )
      {}
      
      ~Cell()
      {
        _particles.clear();
      }
      
      inline int getX() const { return _localX; }
      inline int getY() const { return _localY; }
      
      inline void addParticle( const Particle<floatType>& p ) { _particles.push_back( p ); }
      inline void clearParticles( )
      {
        _particles.clear();
      }
      inline void getParticleList( std::list<Particle<floatType> >* _particleList )
      {
        _particleList = &_particles;
      }
      
    };

  } // namespace memory
} // namespace MDSIM

#endif