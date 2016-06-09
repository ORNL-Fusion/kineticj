#ifndef ROTATION_HPP
#define ROTATION_HPP

#include "c3vec.hpp"

C3Vec rot_XYZ_to_abp ( const C3Vec A_XYZ, const C3Vec bUnit_XYZ, const int direction );
void transpose ( float A[][3] );

#include "rotation.tpp"

#endif
