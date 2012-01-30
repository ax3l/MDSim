#ifndef MPI_COMMUNICATOR_HPP
#define MPI_COMMUNICATOR_HPP

#include <mpi.h>
#include <iostream>

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
      void
      sendParticles( std::list<memory::Particle<double> >& p,
                     const unsigned int direction,
                     const bool posOnly)
      {
        std::cout << "Error: Not implemented!" << std::endl;
        
        //if( ( direction & this->Top ) == this->Top )
        
        /// \todo ....
        /// use different tag for posOnly true or false
        /// create long array (or std::vector) with doubles
        ///
        /// Send:
        ///   MPI_CART_SHIFT
        ///   MPI_ISEND
        ///   -> return the ISend Request as a handle and allow
        ///      the user to call a barrier later with it
        /// 
        /// New Method: Receive
        ///   MPI_Probe
        ///   MPI_Get_Count auf MPI_DOUBLE
        ///   MPI_Recv
      }
      
    };

  } // namespace communicator
} // namespace MDSIM

#endif // MPI_COMMUNICATOR_HPP