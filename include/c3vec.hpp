#ifndef C3VEC_HPP
#define C3VEC_HPP

#include <complex> 
#include <cmath>
#include <vector>
#include <iostream>
#include <assert.h>
#include <algorithm> // for std::max_element
#include "cparticle.hpp"

#if defined(__CUDACC__) || defined(__THRUST)
#include <thrust/complex.h>
#endif

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

                PRAGMA
                HOST DEVICE
				C3 ( const C3 &C3_2 ) {c1=C3_2.c1;c2=C3_2.c2;c3=C3_2.c3;};

                PRAGMA
                template <typename T2>
                HOST DEVICE
				C3 ( const C3<T2> &C3_2 ) {c1=C3_2.c1;c2=C3_2.c2;c3=C3_2.c3;};

                // These assignment operators are implemented here because
                // I was unable to get them to work when the implementation
                // was in the c3vec.tpp file.

                PRAGMA
                HOST DEVICE
                C3& operator=(const C3& rhs)
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
                C3& operator=(const C3<T2>& rhs)
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
    HOST DEVICE
    C3<T> operator() (CParticle &p, C3<T> &field) {
        C3<float> thisVel_XYZ(p.v_c1, p.v_c2, p.v_c3);
        C3<T> result = cross(thisVel_XYZ,field);
        return result;
    }
};

struct doDotProduct 
{
    HOST 
    std::complex<float> operator() (C3<std::complex<float> > &a, C3<float> &b) {
        std::complex<float> result = dot(a,b);
        return result;
    }
};

#if defined(__CUDACC__) || defined(__THRUST)
struct doDotProduct_device
{
    HOST DEVICE
    thrust::complex<float> operator() (C3<thrust::complex<float> > &a, C3<float> &b) {
        thrust::complex<float> result = dot(a,b);
        return result;
    }
};
#endif

template <typename T>
struct runningIntegral
{
    float intFactor;
    runningIntegral( float _intFactor ) : intFactor(_intFactor) {}

    PRAGMA
    HOST DEVICE
    T operator() (T &integral, T &y) {
        T result = integral + intFactor * y;
        return result;
    }
};

template <typename T>
struct multiplyByChargeOverMass
{
    PRAGMA
    HOST DEVICE
    T operator() (T &integral, CParticle &p) {
        T result = -(float)(p.q / p.m) * integral;
        return result;
    }
};

template <typename T>
struct multiplyByCharge
{
    PRAGMA
    HOST DEVICE
    T operator() (T &x, CParticle &p) {
        T result = -(float)(p.q) * x; // Why is this minus sign here?
        return result;
    }
};

#endif
