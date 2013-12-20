#pragma once

#include <mpi.h>
#include <iostream>
#include <cstdlib>

#include "../simulation_defines.hpp"
#include "domain.hpp"
#include "../communicator/mpi_communicator.hpp"

using namespace std;
using namespace MDSIM;


memory::Domain<double> init_domain()
{

  communicator::MPI_Communicator& comm =
    communicator::MPI_Communicator::getInstance();

  const int myXCells  = simParams::cellsX;
  const int myXOffset = 0;
        int myYCells  = floor( double(simParams::cellsY) / comm.getSize() );
  const int myYOffset = myYCells * comm.getRank();

  if( simParams::cellsY < comm.getSize() )
  {
    std::cout << "Error: CellsY( " << simParams::cellsY
              << " ) < Size( " << comm.getSize() << " )"
              << std::endl;
    exit( EXIT_FAILURE );
  }

  // expand at last cpu if not Mod 0
  if( comm.getRank() == comm.getSize() -1 && comm.getSize()*myYCells != simParams::cellsY )
  {
    const int expandYatLastRank = simParams::cellsY - comm.getSize()*myYCells;
    myYCells += expandYatLastRank;
    std::cout << "Info: Added " << expandYatLastRank
              << " YCellLines to last Rank because " << simParams::cellsY
              << "%" << comm.getSize() << " != 0 (= cellsY global % size )"
              << std::endl;
  }

  std::cout << "Info: MyDomain  cellNr[" << myXOffset << " - "
            << myXOffset + myXCells -1 << "; "
            << myYOffset << " - " << myYOffset + myYCells -1 << "]"
            << " at Rank " << comm.getRank()
            << std::endl;

  memory::Domain<double> myDomain( myXCells, myYCells, myXOffset, myYOffset );

  std::cout << "Info: MyDomain posReal[" << myDomain.getFirstCellPos().x << " - "
            << myDomain.getLastCellPos().x << "; "
            << myDomain.getFirstCellPos().y << " - " << myDomain.getLastCellPos().y << ")"
            << " at Rank " << comm.getRank()
            << std::endl;

  return myDomain;
}
