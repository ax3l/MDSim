#ifndef DOMAIN_HPP
#define DOMAIN_HPP

#include <vector>
#include <list>
#include <cmath>
#include <iostream>

#include "cell.hpp"
#include "particle.hpp"
#include "vector3D.hpp"
#include "../simulation_defines.hpp"

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
              const int y0 );

      ~Domain( );
      
      enum {
        NotInDomain = -1
      };
      
      enum {
        Top    = 1u, // Y: Line 0
        Bottom = 2u, // Y: Line ( _totalSizeY -1 )
        Left   = 5u, // X: Column 0
        Right  = 7u  // X: Column ( _totalSizeX -1 )
      };
      
      enum {
        XPeriodic = 1u,
        YPeriodic = 2u
      };
      
      enum {
        AreaGhost  = 1u,
        AreaBorder = 2u//,
        //AreaCore   = 5u
      };
      
      /// Get the Real Position of the Beginning of this Domain
      ///
      /// \return vector3D<floatType> with offset for the first Cell in the Domain
      ///
      vector3D<floatType> getFirstCellPos( const bool includeGhosts = false ) const;
      
      /// Get the Real Position of the End of this Domain
      ///
      /// \return vector3D<floatType> with the position of the last cells "right" border
      ///
      vector3D<floatType> getLastCellPos( const bool includeGhosts = false ) const;
      
      /// Get the Cell number for a position
      ///
      /// \param[in] vector3D<floatType> global position
      /// \param[out] vector3D<floatType> local position in cell if
      ///                                 positionToLocal is set to true
      /// \param[out] bool isGhost the position is a ghost cell
      /// \param[in] bool positionToLocal set the position to the local position
      ///                 within this cell
      /// \return int between [0; NrOfNonGhostCellsInDomain-1] or NotInDomain
      int getCellNr( vector3D<floatType>& pos,
                            bool& isGhost,
                            const bool positionToLocal = false ) const;

      /// Delete the particles within ghost cells
      ///
      void clearGhostCells( );
      
      /// Get Particles within an area (f.e. Domain-leavers, ghosts, border)
      ///
      /// \param[out] std::list<Particle<floatType>* > list of pointers to
      ///                                              particles
      /// \param[in] unsigned int direction; Top, Bottom, Left or Right
      //  \param[in] unsigned int area; AreaGhost or AreaBorder
      ///
      void getArea( std::list<Particle<floatType>* >& p,
                    const unsigned int direction,
                    const unsigned int area );
      
      /// Mark Particles for Removal from this Domain
      ///
      /// Marks particles as inactive.
      /// Removal \see removeParticles
      ///
      /// \param[in] std::list<Particle<floatType>* > list of pointers to
      ///                                             particles
      ///
      //void markParticlesForRemoval( std::list<Particle<floatType>* >& p );
      
      /// Move Particles from Ghost to opposite border area
      ///
      /// Move Particles from Ghost to opposite border area with a specified
      /// direction
      ///
      /// \param[in] std::list<Particle<floatType>* > list of pointers to
      ///                                             particles
      /// \param[in const unsigned int periodic, XPeriodic or YPeriodic
      ///
      void moveInnerDomainPeriodic( std::list<Particle<floatType>* >& p,
                                    const unsigned int periodic );
      
      /// Remove Particles from this Domain
      ///
      /// Removes inaktive marked particles.
      ///
      //void removeParticles( );

      /// \todo this... and impl., too?
      ///
      void mapParticlesToCells( );

      /// Add a particle to this Domain
      ///
      /// Adds the particle to the corresponding cell within this domain.
      ///
      /// \param[in] Particle<floatType> particle with attributes
      /// \param[in] bool allowGhost, allows to add a particle to the ghost
      ///                             cell regions
      ///
      void addParticle( const Particle<floatType>& p,
                        const bool allowGhost = false );

      /// Resets the Forces for each Particle
      ///
      /// (without ghosts)
      ///
      void resetForces( );
      
      /// Calculate the NxN Forces
      ///
      /// \todo extract the Physics to physics/ and into classes
      /// \todo calculate forces with particles in neighbor cells
      ///
      void calculateForces( );
      
      /// Move Particles
      ///
      /// Read the particles force vector and move it (non-SRT, Newtonian)
      ///
      void moveParticles( );
      
      /// Write out the global position for each particle
      ///
      void coutParticlePos( );

    };
    
#include "domain.tpp"

  } // namespace memory
} // namespace MDSIM

#endif // DOMAIN_HPP