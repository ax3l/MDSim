#ifndef MPI_COMMUNICATOR_HPP
#define MPI_COMMUNICATOR_HPP

#include <mpi.h>

#include <iostream>
#include <map>
#include <utility> // for pair
#include <vector>

#include "../memory/particle.hpp"
#include "../memory/vector3D.hpp"

namespace MDSIM
{
  namespace communicator
  {
    /// Singleton Class for MPI Communicator
    ///
    class MPI_Communicator
    {
    private:
      int _rank, _size;
      MPI_Comm _comm_cart;
      std::map<MPI_Request, std::vector<double>* > _dataOut;
      
      // Hide Constructors
      MPI_Communicator()
      {
        MPI_Init( NULL, NULL ); /* starts MPI */
        MPI_Comm_rank( MPI_COMM_WORLD, &_rank ); /* get current process id */
        MPI_Comm_size( MPI_COMM_WORLD, &_size ); /* get number of processes */
        
        const int ndims = 1;
        int nrRanks[ ndims ] = { _size };
        int periods[ ndims ] = { true };
        const int reorder = false; // could be done with appropriate hardware supp.

        MPI_Cart_create( MPI_COMM_WORLD,
                         ndims,
                         nrRanks,
                         periods,
                         reorder,
                         &_comm_cart );
        
        std::cout << "Info: MPI initialized..." << std::endl;
      }
      
      MPI_Communicator(const MPI_Communicator&) {}
      MPI_Communicator& operator=(const MPI_Communicator&) {}
      
      // Hide Destructor
      ~MPI_Communicator()
      {
        MPI_Finalize();
        std::cout << "Info: MPI finalized..." << std::endl;
      }
      
    public:
      typedef MPI_Request handle;
      
      enum Direction {
        Top    = 1u, // Y: Line 0
        Bottom = 2u  //, // Y: Line ( _totalSizeY -1 )
        //Left   = 4u, // X: Column 0
        //Right  = 8u  // X: Column ( _totalSizeX -1 )
      };
      
      /// Create or get this Communicator
      ///
      /// \return instance as static MPI_Communicator
      ///
      static inline MPI_Communicator& getInstance()
      {
        static MPI_Communicator instance; 
        return instance;
      }
      
      int getRank( ) const
      {
        return _rank;
      }
      
      int getSize( ) const
      {
        return _size;
      }
      
      /// Send Particles to an other node
      ///
      /// \param[in] std::list<memory::Particle<double> > byRef with global
      ///                                                 positions
      /// \param[in] const unsigned int direction, Top or Bottom
      /// \param[in] bool posOnly send only positions or complete particle
      ///
      handle
      sendParticles( std::list<memory::Particle<double> >& p,
                     const unsigned int direction,
                     const bool posOnly)
      {
        std::cout << "Error: Not implemented!" << std::endl;
        
        int elemPerParticle;
        if( posOnly )
          elemPerParticle = 3; // 3xPos
        else
          elemPerParticle = 7; // 3xPos, 3xVel, 1xMass
        
        std::vector<double>* sendData = new std::vector<double>;
        sendData->reserve( p.size()*elemPerParticle );
        
        for( std::list<memory::Particle<double> >::iterator it = p.begin();
             it != p.end(); it++ )
        {
          if( posOnly )
          {
            const memory::vector3D<double> r( it->getPosition() );
            sendData->push_back( r.x );
            sendData->push_back( r.y );
            sendData->push_back( r.z );
          }
          else
          {
            const memory::vector3D<double> r( it->getPosition() );
            const memory::vector3D<double> v( it->getVelocity() );
            double m = it->getMass();
            sendData->push_back( r.x );
            sendData->push_back( r.y );
            sendData->push_back( r.z );
            sendData->push_back( v.x );
            sendData->push_back( v.y );
            sendData->push_back( v.z );
            sendData->push_back( m );
          }
        }
        
        //if( ( direction & this->Top ) == this->Top )
        /// \todo get recv rank
        
        int dest = 0; // todo!
        int tag  = int(posOnly);
        MPI_Request mySend;
        
        MPI_Isend( &(*sendData->begin()),
                   sendData->size(),
                   MPI_DOUBLE,
                   dest,
                   tag,
                   MPI_COMM_WORLD,
                   &mySend );

        _dataOut.insert(
          std::pair<MPI_Request, std::vector<double>* >( mySend, sendData ) );
        
        return mySend;
      }
      
      /// Finish Send and free buffers
      ///
      void
      finishSend( handle h )
      {
        std::cout << "Error: Not implemented!" << std::endl;
        /// MPI_Test oder _Probe
        ///
        /// map call destructor for value to handle
        /// then erase map element with handle
      }
      
      /// Receive Particles
      ///
      void
      receiveParticles( std::vector<memory::Particle<double> >& p )
      {
        std::cout << "Error: Not implemented!" << std::endl;
        
        /// New Method: Receive
        ///   MPI_Probe
        ///   MPI_Get_Count auf MPI_DOUBLE
        ///   MPI_Recv
        /// 
      }
      
    };

  } // namespace communicator
} // namespace MDSIM

#endif // MPI_COMMUNICATOR_HPP