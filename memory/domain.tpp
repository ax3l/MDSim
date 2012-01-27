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
    const vector3D<floatType> offset( double( xCell ) * simParams::cutoff,
                                      double( yCell ) * simParams::cutoff,
                                      0.0 );
    pos -= offset;
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
void
Domain<floatType>::getArea( std::list<Particle<floatType>* >& p,
                            const unsigned int direction,
                            const unsigned int area )
{
  std::cout << "ERROR: not implemented!" << std::endl;
}

template <typename floatType>
void
Domain<floatType>::moveInnerDomainPeriodic( std::list<Particle<floatType>* >& p,
                                            const unsigned int periodic )
{
  std::cout << "ERROR: not implemented!" << std::endl;
}

template <typename floatType>
void
Domain<floatType>::mapParticlesToCells( )
{
  std::list<Particle<floatType> >* curParticleList;
  typename std::list<Particle<floatType> >::iterator p;

  // walk trough all cells in this domain which (with ghosts)
  for( int y = 0; y < _totalSizeY; y++ )
    for( int x = 0; x < _totalSizeX; x++ )
    {
      curParticleList = _cellMatrix.at( y * _totalSizeX + x ).getParticleList( );

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
    std::cout << "ERROR: Position( " << p.getPosition( ).x << " "
      << p.getPosition( ).y << " " << p.getPosition( ).z
      << " ) not in Domain! (isGhost: " << isGhost << ")"
      << std::endl;
  }
  else
  {
    p_local.setPosition( r_local );
    _cellMatrix.at( cellNr ).addParticle( p_local );
    //std::cout << "Added particle to cell nr. " << cellNr << std::endl;
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
void
Domain<floatType>::calculateForces( )
{
  std::list<Particle<floatType> >* curParticleList;
  typename std::list<Particle<floatType> >::iterator p1, p2;

  // walk trough all cells in this domain which are not ghosts
  for( int y = 1; y < _totalSizeY - 1; y++ )
    for( int x = 1; x < _totalSizeX - 1; x++ )
    {
      curParticleList = _cellMatrix.at( y * _totalSizeX + x ).getParticleList( );

      // calculcate in-cell-forces - NxN
      for( p1 = curParticleList->begin( ); p1 != curParticleList->end( ); p1++ )
      {
        p1->resetForce( );
        // calculcate in-cell-forces - NxN
        for( p2 = p1; p2 != curParticleList->end( ); p2++ )
        {
          if( p2 == p1 ) continue;

          // F = G*m1*m2/r^2
          vector3D<floatType> r( p2->getPosition( ) - p1->getPosition( ) );
          vector3D<floatType> er( r / sqrt( r.abs2( ) ) );

          const floatType fAbsTmp = simParams::G * p1->getMass( ) * p2->getMass( );
          const floatType fAbs = fAbsTmp / r.abs2( );

          vector3D<floatType> force( fAbs * er.x, fAbs * er.y, 0.0 );
          //std::cout << "Er: " << force.x << " " << force.y << " " << force.z << std::endl;
          //std::cout << "G: " << simParams::G << std::endl;

          p1->addForce( force );
          p2->addForce( force * ( -1.0 ) );
        }

        /// \todo forces with particles in neighbor cells

      }
    }
}

template <typename floatType>
void
Domain<floatType>::moveParticles( )
{
  std::list<Particle<floatType> >* curParticleList;
  typename std::list<Particle<floatType> >::iterator p;

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
        p->addPosition( p->getVelocity( ) * simParams::dt );
        //std::cout << "Move: " << p->getVelocity().x << " " << p->getVelocity().y << " " << p->getVelocity().z << std::endl;

      }
    }
}

template <typename floatType>
void
Domain<floatType>::coutParticlePos( )
{
  std::list<Particle<floatType> >* curParticleList;
  typename std::list<Particle<floatType> >::iterator p;

  // walk trough all cells in this domain which are not ghosts
  for( int y = 1; y < _totalSizeY - 1; y++ )
    for( int x = 1; x < _totalSizeX - 1; x++ )
    {
      curParticleList = _cellMatrix.at( y * _totalSizeX + x ).getParticleList( );
      const vector3D<floatType> cellOffset( double( x ) * simParams::cutoff,
                                            double( y ) * simParams::cutoff,
                                            0.0 );

      // edit every velocity and move to next position
      for( p = curParticleList->begin( ); p != curParticleList->end( ); p++ )
      {
        const vector3D<floatType> r( p->getPosition( )
                                     + getFirstCellPos( ) // Domain Offset
                                     + cellOffset ); // Cell in Domain)
        std::cout << r.x << " " << r.y << " " << r.z << std::endl;
      }
    }
}


#endif // DOMAIN_TPP