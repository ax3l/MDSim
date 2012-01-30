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
    };

  } // namespace communicator
} // namespace MDSIM

#endif // MPI_COMMUNICATOR_HPP