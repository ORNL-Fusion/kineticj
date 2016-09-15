#ifndef ROTATION_HPP
#define ROTATION_HPP

#ifdef __CUDACC__
#define HOST __host__ 
#define DEVICE __device__
#else
#define HOST 
#define DEVICE
#endif

#ifdef __CUDACC__
#define PRAGMA #pragma hd_warning_disable 
#else
#define PRAGMA
#endif

#include "c3vec.hpp"

C3<float> rot_XYZ_to_abp ( const C3<float> A_XYZ, const C3<float> bUnit_XYZ, const int direction );
void transpose ( float A[][3] );

#include "rotation.tpp"

#endif
