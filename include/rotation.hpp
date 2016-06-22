#ifndef ROTATION_HPP
#define ROTATION_HPP

#include "c3vec.hpp"

C3<float> rot_XYZ_to_abp ( const C3<float> A_XYZ, const C3<float> bUnit_XYZ, const int direction );
void transpose ( float A[][3] );

#include "rotation.tpp"

#endif
