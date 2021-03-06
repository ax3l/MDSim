#ifndef DOMAIN_TPP
#define DOMAIN_TPP

template <typename floatType>
Domain<floatType>::Domain( const int sizeX,
                           const int sizeY,
                           const int x0,
                           const int y0 ) :
  _totalSizeX( sizeX + 2 ),
  _totalSizeY( sizeY + 2 ),
  _x0( x0 ),
  _y0( y0 )
{
  if( sizeX < 4 )
    std::cout << "ERROR: SizeX must be at least 4 (is: " << sizeX << ")"
              << std::endl;
  if( sizeY < 4 )
    std::cout << "ERROR: SizeY must be at least 4 (is: " << sizeY << ")"
              << std::endl;
  
  _cellMatrix.reserve( sizeX * sizeY );

  for( int x = x0 - 1; x < x0 + sizeX + 1; x++ )
    for( int y = y0 - 1; y < y0 + sizeY + 1; y++ )
    {
      Cell<floatType> c( x, y );
      _cellMatrix.push_back( c );
    }
}

template <typename floatType>
Domain<floatType>::~Domain( )
{
  _cellMatrix.clear( );
}

template <typename floatType>
vector3D<floatType>
Domain<floatType>::getFirstCellPos( const bool includeGhosts ) const
{
  int ghosts = 0;
  if( includeGhosts == true ) ghosts = -1;

  vector3D<floatType> r( ( _x0 + ghosts ) * simParams::cutoff,
                         ( _y0 + ghosts ) * simParams::cutoff,
                         0.0 );
  return r;
}

template <typename floatType>
vector3D<floatType>
Domain<floatType>::getLastCellPos( const bool includeGhosts ) const
{
  int ghosts = 2;
  if( includeGhosts == true ) ghosts = 1;

  vector3D<floatType> r( ( _x0 + _totalSizeX - ghosts ) * simParams::cutoff,
                         ( _y0 + _totalSizeY - ghosts ) * simParams::cutoff,
                         0.0 );
  return r;
}

template <typename floatType>
int
Domain<floatType>::getCellNr( vector3D<floatType>& pos,
                              bool& isGhost,
                              const bool positionToLocal ) const
{
  const int xCell = floor( ( pos.x - getFirstCellPos( true ).x )
                           / simParams::cutoff );
  const int yCell = floor( ( pos.y - getFirstCellPos( true ).y )
                           / simParams::cutoff );

  if( xCell < 0 || xCell >= _totalSizeX ||
      yCell < 0 || yCell >= _totalSizeY )
    return NotInDomain;

  if( xCell == 0 || xCell == _totalSizeX - 1 ||
      yCell == 0 || yCell == _totalSizeY - 1 )
    isGhost = true;

  if( positionToLocal )
  {
    const vector3D<floatType> offset( double( xCell -1 ) * simParams::cutoff
                                      + getFirstCellPos().x,
                                      double( yCell -1 ) * simParams::cutoff
                                      + getFirstCellPos().y,
                                      0.0 );
    pos -= offset;
    //std::cout << "Info: local pos ("
    //          << pos.x << ", " << pos.y << ")" << std::endl;
  }

  return yCell * _totalSizeX + xCell;
}

template <typename floatType>
void
Domain<floatType>::clearGhostCells( )
{
  for( int x = 0; x < _totalSizeX; x++ )
    for( int y = 0; y < _totalSizeY; y++ )
      if( x == 0 || x == _totalSizeX - 1 ||
          y == 0 || y == _totalSizeY - 1 )
        _cellMatrix.at( y * _totalSizeX + x ).clearParticles( );
}

template <typename floatType>
std::list<Particle<floatType> >
Domain<floatType>::getArea( const unsigned int direction,
                            const unsigned int area )
{
  std::list<Particle<floatType> > p;
  
  std::list<Particle<floatType> >* curParticleList;
  typename std::list<Particle<floatType> >::iterator pInCell;
  
  for( int y = 0; y < _totalSizeY; y++ )
  {
    for( int x = 0; x < _totalSizeX; x++ )
    {
      bool ok = false;
      
      // TOP
      if( ( direction & this->Top ) == this->Top 
          && y <= 1 )
      {
        if( ( area & this->AreaGhost ) == this->AreaGhost
          && y == 0 )
          ok = true;
        if( ( area & this->AreaBorder ) == this->AreaBorder
            && y == 1 )
          ok = true;
      }
      
      // BOTTOM
      if( ( direction & this->Bottom ) == this->Bottom
          && y >= _totalSizeY -2 )
      {
        if( ( area & this->AreaGhost ) == this->AreaGhost
            && y == _totalSizeY -1 )
          ok = true;
        if( ( area & this->AreaBorder ) == this->AreaBorder
            && y == _totalSizeY -2 )
          ok = true;
      }
      
      // LEFT
      if( ( direction & this->Left ) == this->Left
          && x <= 1 )
      {
        if( ( area & this->AreaGhost ) == this->AreaGhost
            && x == 0 )
          ok = true;
        if( ( area & this->AreaBorder ) == this->AreaBorder
            && x == 1 )
          ok = true;
      }
      
      // RIGHT
      if( ( direction & this->Right ) == this->Right
          && x >= _totalSizeX -2 )
      {
        if( ( area & this->AreaGhost ) == this->AreaGhost
            && x == _totalSizeX - 1 )
          ok = true;
        if( ( area & this->AreaBorder ) == this->AreaBorder
            && x == _totalSizeX - 2 )
          ok = true;
      }
      
      if( ok == false )
        continue;
      
      //std::cout << "ok! (" << x  << ", "  << x  << ")" << std::endl;
      
      // add particles of this cell
      curParticleList = _cellMatrix.at( y * _totalSizeX + x ).getParticleList( );
      const vector3D<floatType> cellOffset( double( x -1 ) * simParams::cutoff,
                                            double( y -1 ) * simParams::cutoff,
                                            0.0 );
      
      for( pInCell = curParticleList->begin( ); pInCell != curParticleList->end( ); pInCell++ )
      {
        //std::cout << "Info: Found particle!" << std::endl;
        p.push_back( (*pInCell) );
        
        const vector3D<floatType> r( getFirstCellPos( )   // Domain Offset
                                     + cellOffset );      // Cell in Domain
        
        p.back().addPosition( r );
        
        vector3D<floatType> yCheck( p.back().getPosition() );
        
        // y-periodic:
        if( yCheck.y > ( simParams::cellsY -1 )*simParams::cutoff )
          yCheck.y -= simParams::cellsY * simParams::cutoff;
        else if( yCheck.y < simParams::cutoff )
          yCheck.y += simParams::cellsY * simParams::cutoff;
        
        p.back().setPosition( yCheck );
      }
    }
  }
  
  return p;
}

template <typename floatType>
void
Domain<floatType>::moveInnerDomainPeriodic( const unsigned int periodic )
{
  if( ( periodic & this->XPeriodic ) != this->XPeriodic )
    std::cout << "Error: Not implemented!" << std::endl;
  
  std::list<Particle<floatType> >* curParticleList;
  std::list<Particle<floatType> >* moveToParticleList;
  typename std::list<Particle<floatType> >::iterator pInCell;

  // go through LEFT and RIGHT ghosts
  for( int y = 0; y < _totalSizeY; y++ )
    for( int x = 0; x < _totalSizeX; x+= _totalSizeX -1 )
    {
      curParticleList = _cellMatrix.at( y * _totalSizeX + x ).getParticleList( );

      for( pInCell = curParticleList->begin( ); pInCell != curParticleList->end( ); pInCell++ )
      {
        
        if( x == 0 )
        {
          moveToParticleList = _cellMatrix.at( y * _totalSizeX + _totalSizeX -2 ).getParticleList( );
          moveToParticleList->push_back( (*pInCell) );
          //std::cout << "Info: x-periodic to cell ("
          //          << _totalSizeX -2 << ", " << y << ")" << std::endl;
        }
        if( x == _totalSizeX -1 )
        {
          moveToParticleList = _cellMatrix.at( y * _totalSizeX + 1 ).getParticleList( );
          moveToParticleList->push_back( (*pInCell) );
          //std::cout << "Info: x-periodic to cell ("
          //          << _totalSizeX +1 << ", " << y << ")" << std::endl;
        }
      }
      curParticleList->clear();
    }
}

template <typename floatType>
void
Domain<floatType>::createInnerDomainGhosts( const unsigned int periodic )
{
  if( ( periodic & this->XPeriodic ) != this->XPeriodic )
    std::cout << "Error: Not implemented!" << std::endl;
  
  std::list<Particle<floatType> >* curParticleList;
  std::list<Particle<floatType> >* copyToParticleList;
  typename std::list<Particle<floatType> >::iterator pInCell;

  // go through LEFT and RIGHT border
  for( int y = 0; y < _totalSizeY; y++ )
    for( int x = 1; x < _totalSizeX; x+= _totalSizeX -3 )
    {
      curParticleList = _cellMatrix.at( y * _totalSizeX + x ).getParticleList( );

      for( pInCell = curParticleList->begin( ); pInCell != curParticleList->end( ); pInCell++ )
      {
        
        if( x == 1 )
        {
          // select ghost (right)
          copyToParticleList = _cellMatrix.at( y * _totalSizeX + _totalSizeX -1 ).getParticleList( );
          copyToParticleList->push_back( (*pInCell) );
        }
        if( x == _totalSizeX -2 )
        {
          // select ghost (left)
          copyToParticleList = _cellMatrix.at( y * _totalSizeX + 0 ).getParticleList( );
          copyToParticleList->push_back( (*pInCell) );
        }
      }
    }
}

template <typename floatType>
void
Domain<floatType>::mapParticlesToCells( )
{
  std::list<Particle<floatType> >* curParticleList;
  std::list<Particle<floatType> >* moveToParticleList;
  typename std::list<Particle<floatType> >::iterator p;
  
  int xOff;
  int yOff;

  // walk trough all cells in this domain which (with ghosts)
  for( int y = 0; y < _totalSizeY; y++ )
    for( int x = 0; x < _totalSizeX; x++ )
    {
      curParticleList = _cellMatrix.at( y * _totalSizeX + x ).getParticleList( );

      // check for <0. and >=1. in local pos of particle
      for( p = curParticleList->begin( ); p != curParticleList->end( ); p++ )
      {
        vector3D<floatType> r_local( p->getPosition( ) );
        
        xOff = 0;
        yOff = 0;
        
        if( r_local.x >= 1.0 )
          xOff=1;
        if( r_local.x <  0.0 )
          xOff=-1;
        if( r_local.y >= 1.0 )
          yOff=1;
        if( r_local.y <  0.0 )
          yOff=-1;
        
        if( x + xOff >= _totalSizeX || x + xOff < 0 )
          std::cout << "ERROR: x even out of Ghost" << std::endl;
        if( y + yOff >= _totalSizeY || y + yOff < 0 )
          std::cout << "ERROR: y even out of Ghost" << std::endl;
        
        if( xOff != 0 || yOff != 0 )
        {
          //std::cout << "Info: map particle to next cell ("
          //          << ( x + xOff ) << ", " << ( y + yOff ) << ")" << std::endl;
          
          vector3D<floatType> r_Off( -1.0*double(xOff), -1.0*double(yOff), 0.0 );
          p->addPosition( r_Off );
          
          moveToParticleList = _cellMatrix.at( ( y + yOff ) * _totalSizeX + ( x + xOff ) ).getParticleList( );
          moveToParticleList->push_back( ( *p ) );

          p = curParticleList->erase( p );
          if( p != curParticleList->begin( ) )
            p--;
        }
      }
    }
}

template <typename floatType>
void
Domain<floatType>::addParticle( const Particle<floatType>& p,
                                const bool allowGhost )
{
  Particle<floatType> p_local( p );
  vector3D<floatType> r_local( p_local.getPosition( ) );

  bool isGhost = false;
  const int cellNr = getCellNr( r_local, isGhost, true );

  if( cellNr == NotInDomain ||
      ( isGhost == true && allowGhost == false ) )
  {
    communicator::MPI_Communicator& comm =
      communicator::MPI_Communicator::getInstance( );
    int rank = comm.getRank( );
    
    std::cout << "ERROR: Position( " << p.getPosition( ).x << " "
      << p.getPosition( ).y << " " << p.getPosition( ).z
      << " ) not in Domain! (isGhost: " << isGhost << ")"
      << " Rank(" << rank << ")"
      << std::endl;
  }
  else
  {
    p_local.setPosition( r_local );
    _cellMatrix.at( cellNr ).addParticle( p_local );
    
    //communicator::MPI_Communicator& comm =
    //  communicator::MPI_Communicator::getInstance();
    //int rank = comm.getRank();
        
    //std::cout << "Added particle to cell nr. " << cellNr
    //          << " rank(" << rank << ")" << std::endl;
  }
}

template <typename floatType>
void
Domain<floatType>::addParticle( std::vector<memory::Particle<double> >& p,
                                 const bool allowGhost )
{
  for( std::vector<memory::Particle<double> >::iterator itNewP = p.begin( );
       itNewP != p.end( );
       itNewP++ )
  {
    addParticle( ( *itNewP ), allowGhost );
  }
}

template <typename floatType>
void
Domain<floatType>::resetForces( )
{
  std::list<Particle<floatType> >* curParticleList;
  typename std::list<Particle<floatType> >::iterator p;

  // walk trough all cells in this domain which are not ghosts
  for( int y = 1; y < _totalSizeY - 1; y++ )
    for( int x = 1; x < _totalSizeX - 1; x++ )
    {
      curParticleList = _cellMatrix.at( y * _totalSizeX + x ).getParticleList( );

      // walk through all particles within this cell
      for( p = curParticleList->begin( ); p != curParticleList->end( ); p++ )
      {
        p->resetForce( );
      }
    }
}

template <typename floatType>
template <typename ForceModel>
void
Domain<floatType>::calculateForces( const unsigned int area )
{
  int x0, xE;
  int y0, yE;
  
  if( area == this->AreaHardCore )
  {
    x0 = 3;
    y0 = 3;
    xE = _totalSizeX - 3;
    yE = _totalSizeY - 3;
  }
  else if( area == ( this->AreaBorder | this->AreaCore ) )
  {
    x0 = 1;
    y0 = 1;
    xE = _totalSizeX - 1;
    yE = _totalSizeY - 1;
  }
  else if( ( area & this->AreaBorder )   == this->AreaBorder &&
           ( area & this->AreaCore )     == this->AreaCore &&
           ( area & this->AreaHardCore ) == this->AreaHardCore &&
           ( area & this->AreaGhost )    != this->AreaGhost     )
  {
    x0 = 1;
    y0 = 1;
    xE = _totalSizeX - 1;
    yE = _totalSizeY - 1;
  }
  else
    std::cout << "Error: Not implemented!" << std::endl;
  
  std::list<Particle<floatType> >* curParticleList;
  std::list<Particle<floatType> >* neighborParticleList;
  typename std::list<Particle<floatType> >::iterator p1, p2;
  
  ForceModel forceModel;

  // walk trough all cells in this domain which are not ghosts
  for( int y = y0; y < yE; y++ )
    for( int x = x0; x < xE; x++ )
    {
      // Not HardCore?
      if( ( area & this->AreaHardCore ) != this->AreaHardCore &&
          ( x >= 3 && y >= 3 && x <_totalSizeX - 3 && y < _totalSizeY - 3 ) )
        continue;
      
      curParticleList = _cellMatrix.at( y * _totalSizeX + x ).getParticleList( );

      // calculcate in-cell-forces - NxN
      for( p1 = curParticleList->begin( ); p1 != curParticleList->end( ); p1++ )
      {
        // walk through next neighbor cells
        for( int yN = y -1; yN <= y +1; yN++ )
          for( int xN = x -1; xN <= x +1; xN++ )
          {
            // my Cell
            if( yN == y && xN == x )
            {
              // calculcate in-cell-forces - NxN
              for( p2 = p1; p2 != curParticleList->end( ); p2++ )
              {
                if( p2 == p1 ) continue;
                
                vector3D<floatType> cellOff( 0.0, 0.0, 0.0 );
                vector3D<floatType> force( forceModel( p1, p2, cellOff ) );

                p1->addForce( force );
                p2->addForce( force * ( -1.0 ) );
                //std::cout << "Force2N: " << force.abs2() << " " << p1->getMass() << std::endl;
              }
            }
            else // neighbor cells
            {
              neighborParticleList = _cellMatrix.at( yN * _totalSizeX + xN ).getParticleList( );
              for( p2 = neighborParticleList->begin( ); p2 != neighborParticleList->end( ); p2++ )
              {
                vector3D<floatType> cellOff( xN - x, yN - y, 0.0 );
                vector3D<floatType> force( forceModel( p1, p2, cellOff ) );

                p1->addForce( force );
                //std::cout << "Force2N: " << force.abs2() << " " << p1->getMass() << std::endl;
                /// \todo speed improvement:
                ///       implement newtons third law, which requires an
                ///       intelligent mapper/walk through domain
              }
            }
            
          }

      }
    }
}

template <typename floatType>
void
Domain<floatType>::moveParticles( )
{
  std::list<Particle<floatType> >* curParticleList;
  typename std::list<Particle<floatType> >::iterator p;
  
  const memory::vector3D<floatType> nullVec( 0.0, 0.0, 0.0 );

  // walk trough all cells in this domain which are not ghosts
  for( int y = 1; y < _totalSizeY - 1; y++ )
    for( int x = 1; x < _totalSizeX - 1; x++ )
    {
      curParticleList = _cellMatrix.at( y * _totalSizeX + x ).getParticleList( );

      // edit every velocity and move to next position
      for( p = curParticleList->begin( ); p != curParticleList->end( ); p++ )
      {
        // Newton: F = dp/dt = m * dv/dt
        // dv = F / m * dt
        p->addVelocity( p->getForce( ) / p->getMass( ) * simParams::dt );

        // Move Particle: ds = v * dt
        memory::vector3D<floatType> dr( p->getVelocity( ) * simParams::dt );
        
        // Ok, SRT / too small cell size fix:
        if( dr.abs2() >= simParams::cutoff*simParams::cutoff )
        {
          dr.x = 0.1 * simParams::cutoff;
          dr.y = 0.1 * simParams::cutoff;
          p->setVelocity( nullVec );
        }
        
        p->addPosition( dr );
        //std::cout << "Move: " << p->getVelocity().x << " " << p->getVelocity().y << " " << p->getVelocity().z << std::endl;

      }
    }
}

template <typename floatType>
void
Domain<floatType>::coutParticlePos( const int timestep, const bool ghostsToo )
{
  std::list<Particle<floatType> >* curParticleList;
  typename std::list<Particle<floatType> >::iterator p;

  communicator::MPI_Communicator& comm =
    communicator::MPI_Communicator::getInstance( );
  int rank = comm.getRank( );
  
  int ghostOffset = 0;
  if( !ghostsToo )
    ghostOffset = 1;

  // walk trough all cells in this domain which are not ghosts
  for( int y = ghostOffset; y < _totalSizeY - ghostOffset; y++ )
    for( int x = ghostOffset; x < _totalSizeX - ghostOffset; x++ )
    {
      curParticleList = _cellMatrix.at( y * _totalSizeX + x ).getParticleList( );
      const vector3D<floatType> cellOffset( double( x -1 ) * simParams::cutoff,
                                            double( y -1 ) * simParams::cutoff,
                                            0.0 );

      // edit every velocity and move to next position
      for( p = curParticleList->begin( ); p != curParticleList->end( ); p++ )
      {
        const vector3D<floatType> r( p->getPosition( )
                                     + getFirstCellPos( ) // Domain Offset
                                     + cellOffset );      // Cell in Domain
        
        int isGhost = int( y==0 || y == _totalSizeY -1 ||
                           x==0 || x == _totalSizeX -1   );
        std::cout << r.x << " " << r.y << " " << r.z
                  << " ghost(" << isGhost << ")"
                  << " cell(" << x << ", " << y << ")"
                  << " rank(" << rank << ")" << " time(" << timestep << ")"
                  << std::endl;
      }
    }
}

template <typename floatType>
void
Domain<floatType>::coutParticleNum( )
{
  std::list<Particle<floatType> >* curParticleList;
  typename std::list<Particle<floatType> >::iterator p;
  
  communicator::MPI_Communicator& comm =
    communicator::MPI_Communicator::getInstance( );
  int rank = comm.getRank( );
  
  int anzPart = 0;

  // walk trough all cells in this domain which are not ghosts
  for( int y = 1; y < _totalSizeY - 1; y++ )
    for( int x = 1; x < _totalSizeX - 1; x++ )
    {
      curParticleList = _cellMatrix.at( y * _totalSizeX + x ).getParticleList( );

      // edit every velocity and move to next position
      for( p = curParticleList->begin( ); p != curParticleList->end( ); p++ )
      {
        anzPart++;
      }
    }
  
  int anzGlobal = comm.globalSumReduce( anzPart );
  
  if( rank == 0 )
    std::cout << "Particles global " << anzGlobal << std::endl;
}


#endif // DOMAIN_TPP