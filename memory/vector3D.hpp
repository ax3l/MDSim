#ifndef VECTOR3D_HPP
#define VECTOR3D_HPP

namespace MDSIM {
  namespace memory {
    
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
        x(a), y(a), z(a) {};

      vector3D ( const vector3D<floatType>& v ) :
        x(v.x), y(v.y), z(v.z) {};
        
      ~vector3D() {}
    };

  } // namespace memory
} // namespace MDSIM

#endif // VECTOR3D_HPP