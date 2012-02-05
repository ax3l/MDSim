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
      //MPI_Comm _comm_cart;
      std::map<MPI_Request, std::vector<double>* > _dataOut;
      
      // Hide Constructors
      MPI_Communicator()
      {
        MPI_Init( NULL, NULL ); /* starts MPI */
        MPI_Comm_rank( MPI_COMM_WORLD, &_rank ); /* get current process id */
        MPI_Comm_size( MPI_COMM_WORLD, &_size ); /* get number of processes */
        
        //const int ndims = 1;
        //int nrRanks[ ndims ] = { _size };
        //int periods[ ndims ] = { true };
        //const int reorder = false; // could be done with appropriate hardware supp.

        //MPI_Cart_create( MPI_COMM_WORLD,
        //                 ndims,
        //                nrRanks,
        //                 periods,
        //                 reorder,
        //                 &_comm_cart );
        
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
      
      /// Parse Double Buffer from recv to particle list
      /// 
      /// \param[in] recvBuf as std::vector<double>
      /// \param[in] posOnly as bool - 3 attr or 7 attr particle?
      /// \param[out] 
      void
      parseDoubleToParticle( std::vector<double>& recvBuf,
                             bool posOnly,
                             std::vector<memory::Particle<double> >& p )
      {
        for( std::vector<double>::iterator it = recvBuf.begin( );
             it != recvBuf.end( );
             it++ )
        {
          double x = *it;
          ++it;
          double y = *it;
          ++it;
          double z = *it;
          ++it;
          const memory::vector3D<double> r( x, y, z );

          x = *it;
          ++it;
          y = *it;
          ++it;
          z = *it;
          ++it;
          const memory::vector3D<double> v( x, y, z );

          const memory::Particle<double> tmpP( r, v, *it );

          p.push_back( tmpP );
        }
      }
      
    public:
      typedef MPI_Request handle;
      
      handle
      getNullHandle( )
      {
        return MPI_REQUEST_NULL;
      }
      
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
      /// \param[in] const unsigned int direction, Top xor Bottom
      /// \param[in] bool posOnly send only positions or complete particle
      ///
      handle
      sendParticles( std::list<memory::Particle<double> >& p,
                     const unsigned int direction,
                     const bool posOnly)
      { 
        int elemPerParticle;
        if( posOnly )
          elemPerParticle = 4; // 3xPos + 1xMass
        else
          elemPerParticle = 7; // 3xPos, 3xVel, 1xMass
        
        std::vector<double>* sendData = new std::vector<double>;
        
        //if( p.size() > 0 )
        //      std::cout << "Info: Send Particles " << p.size() << std::endl;
        
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
            sendData->push_back( it->getMass() );
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
        
        int dest;
        if( ( direction & this->Top ) == this->Top )
        {
          dest = _rank - 1;
          if( dest == -1 )
            dest += _size;
        }

        if( ( direction & this->Bottom ) == this->Bottom )
        {
          dest = _rank + 1;
          if( dest == _size )
            dest = 0;
        }
        
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
      finishSend( handle& h )
      { 
        if( h = MPI_REQUEST_NULL )
          return;
        
        MPI_Status status;
        MPI_Wait( &h,
                  &status );
        
        // clean destruct vector
        _dataOut[h]->clear();
        delete _dataOut[h];
        // delete map element
        _dataOut.erase( h );
      }
      
      /// Receive Particles
      ///
      void
      receiveParticles( std::vector<memory::Particle<double> >& p,
                        bool includeOnlyParticles )
      {
        if( includeOnlyParticles )
          std::cout << "Error: Not implemented!" << std::endl;
        
        int flagTop = int(false);
        int flagBot = int(false);
        
        bool finishTop = false;
        bool finishBot = false;
        
        int countTop;
        int countBot;
        
        MPI_Status statusTop;
        MPI_Status statusBot;
        
        int tag        = int(false);
        //int tagPosOnly = int(true);

        int srcBot = _rank + 1;
        if( srcBot == _size )
          srcBot = 0;
        
        int srcTop = _rank - 1;
        if( srcTop == -1 )
          srcTop += _size;
        
        while( !flagTop || !flagBot )
        {
          // TOP
          if( !flagTop )
            MPI_Iprobe( srcTop,
                        tag,
                        MPI_COMM_WORLD,
                        &flagTop,
                        &statusTop );
          if( flagTop && !finishTop )
          { 
            MPI_Get_count( &statusTop,
                           MPI_DOUBLE,
                           &countTop );
            
            //if( countTop > 0 )
            //  std::cout << "Info: Receive(" << _rank
            //            << ") from Top(" << srcTop << ") "
            //            << countTop/7 << " particles" << std::endl;
            
            p.reserve( p.size() + countTop/7 );
            
            std::vector<double> newP;
            newP.resize( countTop );
            
            MPI_Recv( &(*newP.begin()),
                      countTop, 
                      MPI_DOUBLE,
                      srcTop,
                      tag,
                      MPI_COMM_WORLD,
                      &statusTop );
            
            // put into particle list
            parseDoubleToParticle( newP,
                                   false,
                                   p );
            newP.clear();
            finishTop = true;
          }
          
          // BOTTOM
          if( !flagBot )
            MPI_Iprobe( srcBot,
                        tag,
                        MPI_COMM_WORLD,
                        &flagBot,
                        &statusBot );
          if( flagBot && !finishBot )
          { 
            MPI_Get_count( &statusBot,
                           MPI_DOUBLE,
                           &countBot );
            
            //if( countBot > 0 )
            //  std::cout << "Info: Receive(" << _rank
            //            << ") from Bot(" << srcBot << ") "
            //            << countBot/7 << " particles" << std::endl;
            
            p.reserve( p.size() + countBot/7 );
            
            std::vector<double> newP;
            newP.resize( countBot );
            
            MPI_Recv( &(*newP.begin()),
                      countBot, 
                      MPI_DOUBLE,
                      srcBot,
                      tag,
                      MPI_COMM_WORLD,
                      &statusBot );
            
            // put into particle list
            parseDoubleToParticle( newP,
                                   false,
                                   p );
            newP.clear();
            finishBot = true;
          }
        }
      }
      
      /// MPI_Reduce with Sum
      ///
      /// \param[in] int i, local input
      /// \return int, global sum
      ///
      int
      globalSumReduce( int& i )
      {
        int o;
        MPI_Reduce( &i,
                    &o,
                    1,
                    MPI_INTEGER,
                    MPI_SUM,
                    0,
                    MPI_COMM_WORLD );
        
        return o;
      }
      
    };

  } // namespace communicator
} // namespace MDSIM

#endif // MPI_COMMUNICATOR_HPP