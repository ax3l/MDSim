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

      inline void clearGhostCells( )
      {
        for( int x = _x0 - 1; x < _x0 + _totalSizeX; x++ )
          for( int y = _y0 - 1; y < _y0 + _totalSizeY; y++ )
            if( x == _x0 - 1 || x == _totalSizeX - 1 ||
                y == _y0 - 1 || y == _totalSizeY - 1 )
              _cellMatrix.at( y * _totalSizeX + x ).clearParticles( );
      }

      inline void mapParticlesToCells( )
      {
        std::list<Particle<floatType> >* curParticleList;
        typename std::list<Particle<floatType> >::iterator p;

        // walk trough all cells in this domain which (with ghosts)
        for( int x = _x0 - 1; x < _x0 + _totalSizeX; x++ )
          for( int y = _y0 - 1; y < _y0 + _totalSizeY; y++ )
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

      inline void addParticle( const Particle<floatType>& p )
      {
        /// \todo this!
        _cellMatrix.at( 4 ).addParticle( p );
      }

      inline void resetForces( )
      {
        std::list<Particle<floatType> >* curParticleList;
        typename std::list<Particle<floatType> >::iterator p;

        // walk trough all cells in this domain which are not ghosts
        for( int x = 1; x < _totalSizeX - 1; x++ )
          for( int y = 1; y < _totalSizeY - 1; y++ )
          {
            curParticleList = _cellMatrix.at( y * _totalSizeX + x ).getParticleList( );

            // calculcate in-cell-forces - NxN
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