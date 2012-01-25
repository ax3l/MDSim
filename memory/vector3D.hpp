#ifndef VECTOR3D_HPP
#define VECTOR3D_HPP

namespace MDSIM
{
  namespace memory
  {

    /// Vector 3D
    ///
    /// A simple 3D floating type with x,y,z
    ///
    /// \tparam floatType float or double for the type of x,y,z

    template <typename floatType>
    class vector3D
    {
    public:
      floatType x;
      floatType y;
      floatType z;

      //       vector3D( ) :
      //         x(0.), y(0.), z(0.) {};

      vector3D( const floatType a ) :
      x( a ), y( a ), z( a )
      {
      };

      vector3D( const floatType a, const floatType b, const floatType c ) :
      x( a ), y( b ), z( c )
      {
      };

      vector3D( const vector3D<floatType>& v ) :
      x( v.x ), y( v.y ), z( v.z )
      {
      };

      ~vector3D( )
      {
      }

      floatType abs2( )
      {
        return this->x * this->x
          + this->y * this->y
          + this->z * this->z;
      }

      static floatType abs2( const vector3D& rhs )
      {
        return rhs.x * rhs.x
          + rhs.y * rhs.y
          + rhs.z * rhs.z;
      }

      vector3D operator+(const vector3D& rhs )
      {
        vector3D result( this->x + rhs.x,
                         this->y + rhs.y,
                         this->z + rhs.z );
        return result;
      }
      
      vector3D operator-(const vector3D& rhs )
      {
        vector3D result( this->x - rhs.x,
                         this->y - rhs.y,
                         this->z - rhs.z );
        return result;
      }
      
      vector3D operator*(const floatType& rhs )
      {
        vector3D result( this->x*rhs,
                         this->y*rhs,
                         this->z*rhs );
        return result;
      }
      
      vector3D operator/(const floatType& rhs )
      {
        vector3D result( this->x/rhs,
                         this->y/rhs,
                         this->z/rhs );
        return result;
      }

      vector3D& operator+=(const vector3D& rhs )
      {
        this->x += rhs.x;
        this->y += rhs.y;
        this->z += rhs.z;
        return *this;
      }
      
      vector3D& operator-=(const vector3D& rhs )
      {
        this->x -= rhs.x;
        this->y -= rhs.y;
        this->z -= rhs.z;
        return *this;
      }

    };

  } // namespace memory
} // namespace MDSIM

#endif // VECTOR3D_HPP