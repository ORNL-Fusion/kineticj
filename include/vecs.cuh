#ifndef VECS_H
#define VECS_H

#include "cusp/complex.h"

// used for complex
using namespace cusp;

class CParticle_PODS {
        public:
                float c1, c2, c3, v_c1, v_c2, v_c3;
                int number;
                float weight;
                int status;
                double m, q;
                double amu;
                int Z;
};

class C3VecI {
		public:
				complex<float> c1, c2, c3;

	        __host__ __device__ C3VecI ();
		__host__ __device__	C3VecI ( complex<float> _c1, complex<float> _c2, complex<float> _c3 );

		__host__ __device__      C3VecI& operator = (const C3VecI &rhs);
		__host__ __device__ 	C3VecI& operator += (const C3VecI &rhs);
		__host__ __device__ 	C3VecI& operator += (const float &rhs);
		__host__ __device__ 	C3VecI& operator -= (const C3VecI &rhs);
		__host__ __device__ 	C3VecI& operator -= (const float &rhs);
		__host__ __device__ 	C3VecI& operator *= (const C3VecI &rhs);
		__host__ __device__ 	C3VecI& operator *= (const float &rhs);
		__host__ __device__ 	C3VecI& operator /= (const C3VecI &rhs);
		__host__ __device__ 	C3VecI& operator /= (const float &rhs);

		__host__ __device__ 	C3VecI operator + (const C3VecI &other);
		__host__ __device__ 	C3VecI operator + (const float &other);
		__host__ __device__ 	C3VecI operator - (const C3VecI &other);
		__host__ __device__ 	C3VecI operator - (const float &other);
		__host__ __device__ 	C3VecI operator * (const C3VecI &other);
		__host__ __device__ 	C3VecI operator * (const float &other);
		__host__ __device__ 	friend C3VecI operator * (const float &other, const C3VecI &rhs);
		__host__ __device__ 	C3VecI operator / (const C3VecI &other);
		__host__ __device__ 	C3VecI operator / (const float &other);
};

class C3Vec {
		public:
				float c1, c2, c3;

            	     __host__ __device__ C3Vec ();
                     __host__ __device__ C3Vec ( float _c1, float _c2, float _c3 );
                     __host__ __device__ C3Vec ( int _arg );

		     __host__ __device__ C3Vec& operator = (const C3Vec &rhs);
		     __host__ __device__ C3Vec& operator = (const float &rhs);
	 	     __host__ __device__ C3Vec& operator += (const C3Vec &rhs);
                     __host__ __device__ C3Vec& operator += (const float &rhs);
                     __host__ __device__ C3Vec& operator -= (const C3Vec &rhs);
                     __host__ __device__ C3Vec& operator -= (const float &rhs);
                     __host__ __device__ C3Vec& operator *= (const C3Vec &rhs);
                     __host__ __device__ C3Vec& operator *= (const float &rhs);
                     __host__ __device__ C3Vec& operator /= (const C3Vec &rhs);
                     __host__ __device__ C3Vec& operator /= (const float &rhs);
                     __host__ __device__ C3Vec operator + (const C3Vec &other);
                     __host__ __device__ C3Vec operator + (const float &other);
                     __host__ __device__ C3Vec operator - (const C3Vec &other);
		     __host__ __device__ C3Vec operator - (const float &other);
                     __host__ __device__ C3Vec operator * (const C3Vec &other);
                     __host__ __device__ C3Vec operator * (const float &other);
                     __host__ __device__ friend C3Vec operator * (const float &other, const C3Vec &rhs);
                     __host__ __device__ C3Vec operator / (const C3Vec &other);
                     __host__ __device__ C3Vec operator / (const float &other);
};

#endif
