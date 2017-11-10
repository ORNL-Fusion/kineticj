#ifndef ROTATION_HPP
#define ROTATION_HPP

#if defined(__CUDACC__) 
#define HOST __host__ 
#define DEVICE __device__
#else
#define HOST 
#define DEVICE
#endif

#if defined(__CUDACC__) 
#define PRAGMA #pragma hd_warning_disable 
#else
#define PRAGMA
#endif

#include <cmath>
#include "c3vec.hpp"

C3<float> rot_XYZ_to_abp ( const C3<float> A_XYZ, const C3<float> bUnit_XYZ, const int direction );
C3<float> rot_axis_angle ( const C3<float> v, const C3<float> u, const float th_deg );

void transpose ( float A[][3] );

#include "rotation.tpp"

#endif
