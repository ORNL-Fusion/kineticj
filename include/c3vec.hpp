#ifndef C3VEC_HPP
#define C3VEC_HPP

#include <complex> 
#include <vector>
#include <iostream>
#include <assert.h>
#include <algorithm> // for std::max_element
#include "cparticle.hpp"

#ifdef __CUDACC__
#include <thrust/complex.h>
#endif

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

template <typename T>
class C3 {
		public:
				T c1, c2, c3;

                PRAGMA
                HOST DEVICE
                C3 (int _C) {c1=_C;c2=_C;c3=_C;};

                PRAGMA
                HOST DEVICE
				C3 () {c1=T(0);c2=T(0);c3=T(0);};

                PRAGMA
                HOST DEVICE
				C3 ( T _c1, T _c2, T _c3 ) {c1=_c1;c2=_c2;c3=_c3;};

                // These assignment operators are implemented here because
                // I was unable to get them to work when the implementation
                // was in the c3vec.tpp file.

                PRAGMA
                HOST DEVICE
                C3<T>& operator=(const C3<T>& rhs)
                {
                    if (this != &rhs) {
                        c1 = rhs.c1;
                        c2 = rhs.c2;
                        c3 = rhs.c3;
                    }
                    return *this;
                }

                // Removing the self-assignment check, this may be dangerours.
                
                PRAGMA
                template <typename T2>
                HOST DEVICE
                C3<T>& operator=(const C3<T2>& rhs)
                {
                    //if (this != &rhs) {
                        c1 = rhs.c1;
                        c2 = rhs.c2;
                        c3 = rhs.c3;
                    //}
                    return *this;
                }


                HOST DEVICE
				C3& operator += (const C3 &rhs);
                HOST DEVICE
				C3& operator += (const float &rhs);
                HOST DEVICE
				C3& operator -= (const C3 &rhs);
                HOST DEVICE
				C3& operator -= (const float &rhs);
                HOST DEVICE
				C3& operator *= (const C3 &rhs);
                HOST DEVICE
				C3& operator *= (const float &rhs);
                HOST DEVICE
				C3& operator /= (const C3 &rhs);
                HOST DEVICE
				C3& operator /= (const float &rhs);


                HOST DEVICE
				C3 operator + (const C3 &other);
                HOST DEVICE
				C3 operator + (const float &other);
                HOST DEVICE
				C3 operator - (const C3 &other);
                HOST DEVICE
				C3 operator - (const float &other);
                HOST DEVICE
				C3 operator * (const C3 &other);
                HOST DEVICE
				C3 operator * (const float &other);

                template <typename T2>
                HOST DEVICE
				friend C3<T2> operator * (const float &other, const C3<T2> &rhs);

                HOST DEVICE
				C3 operator / (const C3 &other);
                HOST DEVICE
				C3 operator / (const float &other);
};

#include "c3vec.tpp" // since templated functions need to be in headers

// Functor to wrap cross and dot 

template <typename T>
struct vCross 
{
    C3<T> operator() (CParticle &p, C3<T> &field) {
        C3<float> thisVel_XYZ(p.v_c1, p.v_c2, p.v_c3);
        C3<T> result = cross(thisVel_XYZ,field);
        return result;
    }
};

template <typename T, typename T2>
struct doDotProduct 
{
    std::complex<float> operator() (C3<T> &a, C3<T2> &b) {
        std::complex<float> result = dot(a,b);
        return result;
    }
};

struct runningIntegral
{
    float intFactor;
    runningIntegral( float _intFactor ) : intFactor(_intFactor) {}

    std::complex<float> operator() (std::complex<float> &integral, std::complex<float> &y) {
        std::complex<float> result = integral + intFactor * y;
        return result;
    }
};

struct multiplyByChargeOverMass
{
    std::complex<float> operator() (std::complex<float> &integral, CParticle &p) {
        std::complex<float> result = -(float)(p.q / p.m) * integral;
        return result;
    }
};

struct multiplyByCharge
{
    std::complex<float> operator() (std::complex<float> &x, CParticle &p) {
        std::complex<float> result = -(float)(p.q) * x;
        return result;
    }
};

#endif
