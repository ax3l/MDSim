#ifndef DOMAIN_HPP
#define DOMAIN_HPP

#include <vector>
#include <list>
#include "cell.hpp"
#include "particle.hpp"

namespace MDSIM {
  namespace memory {
  
    /// Domain
    ///
    /// The local 2D grid of cells
    ///
    /// \tparam floatType Double of Float: for Particle Attributes
    template <typename floatType>
    class Domain
    {
    private:
      std::vector<Cell<floatType> > _cellMatrix;
      
      int _totalSizeX, _totalSizeY, _x0, _y0;
      
    public:
      /// Constructor for Domain
      ///
      /// Creates sizeX+2 times sizeY+2 cells (with ghosts)
      ///
      /// \param[in] sizeX number of cells in x direction
      /// \param[in] sizeY number of cells in y direction
      /// \param[in] x0    x cell number for first non-ghost cell
      /// \param[in] y0    y cell number for first non-ghost cell
      Domain( const int sizeX,
              const int sizeY,
              const int x0,
              const int y0 ) :
        _totalSizeX( sizeX +2 ),
        _totalSizeY( sizeY +2 ),
        _x0( x0 ),
        _y0( y0 )
      {
        _cellMatrix.reserve( sizeX * sizeY );
        
        for( int x = x0-1; x < x0+sizeX+1; x++ )
          for( int y = y0-1; y < y0+sizeY+1; y++ )
          {
            Cell<floatType> c( x, y );
            _cellMatrix.push_back( c );
          }
      }
      
      ~Domain()
      {
        _cellMatrix.clear();
      }
      
      inline void clearGhostCells()
      {
        for( int x = _x0-1; x < _x0 + _totalSizeX; x++ )
          for( int y = _y0-1; y < _y0 + _totalSizeY; y++ )
            if( x == _x0-1 || x == _totalSizeX-1 ||
                y == _y0-1 || y == _totalSizeY-1 )
              _cellMatrix.at( y*_totalSizeX + x ).clearParticles();
      }
      
      inline void mapParticlesToCells()
      {
        std::list<Particle<floatType> >* curParticleList;
        typename std::list<Particle<floatType> >::iterator p;
        
        for( int x = _x0-1; x < _x0 + _totalSizeX; x++ )
          for( int y = _y0-1; y < _y0 + _totalSizeY; y++ )
          {
            _cellMatrix.at( y*_totalSizeX + x ).getParticleList( curParticleList );
            
            // check for <0. and >=1. in local pos of particle
            for( p = *curParticleList.begin(); p!=*curParticleList.end(); p++ )
            {
                std::cout << p->getMass() << std::endl;
                /// \todo check for <0. and >=1. in local pos of particle
            }
            
            /// \todo add to other cell (left, right, top, bottom)
            /// \todo remove particle from list
          }
      }
      
    };

  } // namespace memory
} // namespace MDSIM

#endif // DOMAIN_HPP