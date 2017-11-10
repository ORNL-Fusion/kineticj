#ifndef INTERP_HPP
#define INTERP_HPP

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

#if defined(__CUDACC__) || defined(__THRUST)
#include <thrust/complex.h>
#endif

#include <cmath>

#include "interp.tpp"

#endif




