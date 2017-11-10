#ifndef GETFIELDS_HPP
#define GETFIELDS_HPP

#include "cparticle.hpp"
#include "c3vec.hpp"
#include "interp.hpp"
#include "rotation.hpp"

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

using namespace std;

#if defined(__CUDACC__) || defined(__THRUST)
HOST DEVICE
C3<thrust::complex<float> > getE1orB1_XYZ_fromCYL(CParticle& p_XYZ, float *rVec, C3<thrust::complex<float> > *E1Vec_CYL, int nR, int nPhi);
#endif

HOST
C3<std::complex<float> > getE1orB1_XYZ_fromCYL(CParticle& p_XYZ, float *rVec, C3<std::complex<float> > *E1Vec_CYL, int nR, int nPhi);

#if defined(__CUDACC__) || defined(__THRUST)
HOST DEVICE
C3<thrust::complex<float> > getE1orB1_XYZ_fromXYZ(CParticle& p_XYZ, float *rVec, C3<thrust::complex<float> > *E1Vec_CYL, int nR, float ky, float kz);
#endif

HOST
C3<std::complex<float> > getE1orB1_XYZ_fromXYZ(CParticle& p_XYZ, float *rVec, C3<std::complex<float> > *E1Vec_CYL, int nR, float ky, float kz);


#include "getFields.tpp"

// Functor to wrap these 

#if defined(__CUDACC__) || defined(__THRUST)
struct getPerturbedField_device 
{

    float *r;
    C3<thrust::complex<float> > *field_CYL;
    int nR;
    int nPhi;
    float ky;
    float kz;
    float weight;
    float wrf;
    float t;

    getPerturbedField_device( float *_r, C3<thrust::complex<float> > *_field_CYL, int _nR, int _nPhi, float _ky, float _kz, float _weight, float _wrf, float _t) : 
            r(_r), field_CYL(_field_CYL), nR(_nR), nPhi(_nPhi), ky(_ky), kz(_kz), weight(_weight), wrf(_wrf), t(_t) {}

    HOST DEVICE
    C3<thrust::complex<float> > operator() (CParticle &p) {
        thrust::complex<float> _i(0, 1);
#if CYLINDRICAL_INPUT_FIELDS >=1
        C3<thrust::complex<float> > E1_XYZ = weight * thrust::exp(-_i * wrf * t) * getE1orB1_XYZ_fromCYL(p, r, field_CYL, nR, nPhi);
#else
        C3<thrust::complex<float> > E1_XYZ = weight * thrust::exp(-_i * wrf * t) * getE1orB1_XYZ_fromXYZ(p, r, field_CYL, nR, ky, kz);
#endif
        C3<thrust::complex<float> > field_XYZ = E1_XYZ * (1 - p.status);
        return field_XYZ;
    }
};
#endif

struct getPerturbedField 
{

    float *r;
    C3<std::complex<float> > *field_CYL;
    int nR;
    int nPhi;
    float ky;
    float kz;
    float weight;
    float wrf;
    float t;

    getPerturbedField( float *_r, C3<std::complex<float> > *_field_CYL, int _nR, int _nPhi, float _ky, float _kz, float _weight, float _wrf, float _t) : 
            r(_r), field_CYL(_field_CYL), nR(_nR), nPhi(_nPhi), ky(_ky), kz(_kz), weight(_weight), wrf(_wrf), t(_t) {}

    HOST 
    C3<std::complex<float> > operator() (CParticle &p) {
        std::complex<float> _i(0, 1);
#if CYLINDRICAL_INPUT_FIELDS >=1
        C3<std::complex<float> > E1_XYZ = weight * std::exp(-_i * wrf * t) * getE1orB1_XYZ_fromCYL(p, r, field_CYL, nR, nPhi);
#else
        C3<std::complex<float> > E1_XYZ = weight * std::exp(-_i * wrf * t) * getE1orB1_XYZ_fromXYZ(p, r, field_CYL, nR, ky, kz);
#endif
        C3<std::complex<float> > field_XYZ = E1_XYZ * (1 - p.status);
        return field_XYZ;
    }
};


#endif
