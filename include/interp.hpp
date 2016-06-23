#ifndef INTERP_HPP
#define INTERP_HPP

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

#ifdef __CUDACC__
#include <thrust/complex.h>
#endif

#include "interp.tpp"

#endif




