#ifndef DOMAIN_HPP
#define DOMAIN_HPP

#include <vector>
#include <list>
#include <cmath>
#include <iostream>
#include "cell.hpp"
#include "particle.hpp"
#include "vector3D.hpp"

namespace MDSIM
{
  namespace memory
  {

    /// Domain
    ///
    /// The local 2D grid of cells
    ///
    /// \tparam floatType Double of Float: for Particle Attributes
    ///
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
      ///
      Domain( const int sizeX,
              const int sizeY,
              const int x0,
              const int y0 ) :
      _totalSizeX( sizeX + 2 ),
      _totalSizeY( sizeY + 2 ),
      _x0( x0 ),
      _y0( y0 )
      {
        _cellMatrix.reserve( sizeX * sizeY );

        for( int x = x0 - 1; x < x0 + sizeX + 1; x++ )
          for( int y = y0 - 1; y < y0 + sizeY + 1; y++ )
          {
            Cell<floatType> c( x, y );
            _cellMatrix.push_back( c );
          }
      }

      ~Domain( )
      {
        _cellMatrix.clear( );
      }
      
      enum {
        NotInDomain = -1
      };
      
      /// Get the Real Position of the Beginning of this Domain
      ///
      /// \return vector3D<floatType> with offset for the first Cell in the Domain
      ///
      inline vector3D<floatType> getFirstCellPos( const bool includeGhosts = false ) const
      {
        int ghosts = 0;
        if( includeGhosts == true ) ghosts = -1;
        
        vector3D<floatType> r( ( _x0 + ghosts ) * simParams::cutoff,
                               ( _y0 + ghosts ) * simParams::cutoff,
                               0.0  );
        return r;
      }
      
      /// Get the Real Position of the End of this Domain
      ///
      /// \return vector3D<floatType> with the position of the last cells "right" border
      ///
      inline vector3D<floatType> getLastCellPos( const bool includeGhosts = false ) const
      {
        int ghosts = 2;
        if( includeGhosts == true ) ghosts = 1;
        
        vector3D<floatType> r( ( _x0 + _totalSizeX - ghosts ) * simParams::cutoff,
                               ( _y0 + _totalSizeY - ghosts ) * simParams::cutoff,
                               0.0  );
        return r;
      }
      
      /// Get the Cell number for a position
      ///
      /// \param[in] vector3D<floatType> position
      /// \param[out] bool isGhost the position is a ghost cell
      /// \return int between [0; NrOfNonGhostCellsInDomain-1] or NotInDomain
      inline int getCellNr( const vector3D<floatType>& pos, bool& isGhost ) const
      {
        const int xCell = floor( ( pos.x - getFirstCellPos( true ).x )
                                 / simParams::cutoff );
        const int yCell = floor( ( pos.y - getFirstCellPos( true ).y )
                                 / simParams::cutoff );
        
        if( xCell < 0 || xCell >= _totalSizeX ||
            yCell < 0 || yCell >= _totalSizeY    )
          return NotInDomain;
        
        if( xCell == 0 || xCell == _totalSizeX -1 ||
            yCell == 0 || yCell == _totalSizeY -1   )
          isGhost = true;
        
        return yCell * _totalSizeY + xCell;
      }

      /// Delete the particles within ghost cells
      ///
      inline void clearGhostCells( )
      {
        for( int x = 0; x < _totalSizeX; x++ )
          for( int y = 0; y < _totalSizeY; y++ )
            if( x == 0 || x == _totalSizeX - 1 ||
                y == 0 || y == _totalSizeY - 1 )
              _cellMatrix.at( y * _totalSizeX + x ).clearParticles( );
      }

      inline void mapParticlesToCells( )
      {
        std::list<Particle<floatType> >* curParticleList;
        typename std::list<Particle<floatType> >::iterator p;

        // walk trough all cells in this domain which (with ghosts)
        for( int x = 0; x < _totalSizeX; x++ )
          for( int y = 0; y < _totalSizeY; y++ )
          {
            _cellMatrix.at( y * _totalSizeX + x ).getParticleList( curParticleList );

            // check for <0. and >=1. in local pos of particle
            for( p = curParticleList->begin( ); p != curParticleList->end( ); p++ )
            {
              std::cout << p->getMass( ) << std::endl;
              /// \todo check for <0. and >=1. in local pos of particle
            }

            /// \todo add to other cell (left, right, top, bottom)
            /// \todo remove particle from list
          }
      }

      /// Add a particle to this Domain
      ///
      /// Adds the particle to the corresponding cell within this domain.
      ///
      /// \param[in] Particle<floatType> particle with attributes
      /// \param[in] bool allowGhost, allows to add a particle to the ghost
      ///                             cell regions
      ///
      inline void addParticle( const Particle<floatType>& p,
                               const bool allowGhost = false )
      {
        bool isGhost = false;
        const int cellNr = getCellNr( p.getPosition(), isGhost );
        
        if( cellNr == NotInDomain || 
            ( isGhost == true && allowGhost == false ) )
        {
          std::cout << "ERROR: Position( " << p.getPosition().x << " "
                    << p.getPosition().y << " " << p.getPosition().z
                    << " ) not in Domain! (isGhost: " << isGhost << ")"
                    << std::endl;
        } else
        {
          _cellMatrix.at( cellNr ).addParticle( p );
          //std::cout << "Added particle to cell nr. " << cellNr << std::endl;
        }
      }

      /// Resets the Forces for each Particle
      ///
      /// (without ghosts)
      ///
      inline void resetForces( )
      {
        std::list<Particle<floatType> >* curParticleList;
        typename std::list<Particle<floatType> >::iterator p;

        // walk trough all cells in this domain which are not ghosts
        for( int x = 1; x < _totalSizeX - 1; x++ )
          for( int y = 1; y < _totalSizeY - 1; y++ )
          {
            curParticleList = _cellMatrix.at( y * _totalSizeX + x ).getParticleList( );

            // walk through all particles within this cell
            for( p = curParticleList->begin( ); p != curParticleList->end( ); p++ )
            {
              p->resetForce();
            }
          }
      }
      
      inline void calculateForces( )
      {
        std::list<Particle<floatType> >* curParticleList;
        typename std::list<Particle<floatType> >::iterator p1, p2;

        // walk trough all cells in this domain which are not ghosts
        for( int x = 1; x < _totalSizeX - 1; x++ )
          for( int y = 1; y < _totalSizeY - 1; y++ )
          {
            curParticleList = _cellMatrix.at( y * _totalSizeX + x ).getParticleList( );

            // calculcate in-cell-forces - NxN
            for( p1 = curParticleList->begin( ); p1 != curParticleList->end( ); p1++ )
            {
              p1->resetForce();
              // calculcate in-cell-forces - NxN
              for( p2 = p1; p2 != curParticleList->end( ); p2++ )
              {
                if( p2 == p1 ) continue;
                
                // F = G*m1*m2/r^2
                vector3D<floatType> r( p2->getPosition() - p1->getPosition() );
                vector3D<floatType> er( r / sqrt( r.abs2( ) ) );
                
                const floatType fAbsTmp = simParams::G * p1->getMass( ) * p2->getMass( );
                const floatType fAbs = fAbsTmp / r.abs2();
                
                vector3D<floatType> force( fAbs * er.x, fAbs * er.y, 0.0 );
                //std::cout << "Er: " << force.x << " " << force.y << " " << force.z << std::endl;
                //std::cout << "G: " << simParams::G << std::endl;
                
                p1->addForce( force );
                p2->addForce( force*(-1.0) );
              }
              
              /// \todo forces with particles in neighbor cells
              
            }
          }
      }
      
      inline void moveParticles( )
      {
        std::list<Particle<floatType> >* curParticleList;
        typename std::list<Particle<floatType> >::iterator p;

        // walk trough all cells in this domain which are not ghosts
        for( int x = 1; x < _totalSizeX - 1; x++ )
          for( int y = 1; y < _totalSizeY - 1; y++ )
          {
            curParticleList = _cellMatrix.at( y * _totalSizeX + x ).getParticleList( );

            // edit every velocity and move to next position
            for( p = curParticleList->begin( ); p != curParticleList->end( ); p++ )
            {
              // Newton: F = dp/dt = m * dv/dt
              // dv = F / m * dt
              p->addVelocity( p->getForce() / p->getMass() * simParams::dt );
              
              // Move Particle: ds = v * dt
              p->addPosition( p->getVelocity() * simParams::dt );
              //std::cout << "Move: " << p->getVelocity().x << " " << p->getVelocity().y << " " << p->getVelocity().z << std::endl;
              
            }
          }
      }
      
      inline void coutParticlePos( )
      {
        std::list<Particle<floatType> >* curParticleList;
        typename std::list<Particle<floatType> >::iterator p;

        // walk trough all cells in this domain which are not ghosts
        for( int x = 1; x < _totalSizeX -1; x++ )
          for( int y = 1; y < _totalSizeY -1; y++ )
          {
            curParticleList = _cellMatrix.at( y * _totalSizeX + x ).getParticleList( );

            // edit every velocity and move to next position
            for( p = curParticleList->begin( ); p != curParticleList->end( ); p++ )
            {
              const vector3D<floatType> r( p->getPosition() );
              std::cout << r.x << " " << r.y << " " << r.z << std::endl;
            }
          }
      }

    };

  } // namespace memory
} // namespace MDSIM

#endif // DOMAIN_HPP