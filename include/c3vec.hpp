#ifndef C3VEC_HPP
#define C3VEC_HPP

#include <complex> 
#include <vector>
#include <iostream>
#include <assert.h>
#include <algorithm> // for std::max_element

class C3VecI {
		public:
				std::complex<float> c1, c2, c3;

                C3VecI (int _const) {c1=_const;c2=_const;c3=_const;};
				C3VecI () {c1=std::complex<float>(0.0f,0.0f);c2=std::complex<float>(0.0f,0.0f);c3=std::complex<float>(0.0f,0.0f);};
				C3VecI ( std::complex<float> _c1, std::complex<float> _c2, std::complex<float> _c3 ) {c1=_c1;c2=_c2;c3=_c3;};

				C3VecI& operator = (const C3VecI &rhs);
				C3VecI& operator += (const C3VecI &rhs);
				C3VecI& operator += (const float &rhs);
				C3VecI& operator -= (const C3VecI &rhs);
				C3VecI& operator -= (const float &rhs);
				C3VecI& operator *= (const C3VecI &rhs);
				C3VecI& operator *= (const float &rhs);
				C3VecI& operator /= (const C3VecI &rhs);
				C3VecI& operator /= (const float &rhs);


				C3VecI operator + (const C3VecI &other);
				C3VecI operator + (const float &other);
				C3VecI operator - (const C3VecI &other);
				C3VecI operator - (const float &other);
				C3VecI operator * (const C3VecI &other);
				C3VecI operator * (const float &other);
				friend C3VecI operator * (const float &other, const C3VecI &rhs);
				C3VecI operator / (const C3VecI &other);
				C3VecI operator / (const float &other);
};

class C3Vec {
		public:
				float c1, c2, c3;

                C3Vec (int _const) {c1=_const;c2=_const;c3=_const;};
				C3Vec () {c1=0;c2=0;c3=0;};
				C3Vec ( float _c1, float _c2, float _c3 ) {c1=_c1;c2=_c2;c3=_c3;};

				C3Vec& operator = (const C3Vec &rhs);
				C3Vec& operator += (const C3Vec &rhs);
				C3Vec& operator += (const float &rhs);
				C3Vec& operator -= (const C3Vec &rhs);
				C3Vec& operator -= (const float &rhs);
				C3Vec& operator *= (const C3Vec &rhs);
				C3Vec& operator *= (const float &rhs);
				C3Vec& operator /= (const C3Vec &rhs);
				C3Vec& operator /= (const float &rhs);

				C3Vec operator + (const C3Vec &other);
				C3Vec operator + (const float &other);
				C3Vec operator - (const C3Vec &other);
				C3Vec operator - (const float &other);
				C3Vec operator * (const C3Vec &other);
				C3Vec operator * (const float &other);
				friend C3Vec operator * (const float &other, const C3Vec &rhs);
				C3Vec operator / (const C3Vec &other);
				C3Vec operator / (const float &other);
};

C3Vec operator* ( const float &other, const C3Vec &rhs );
C3VecI operator* ( const float &other, const C3VecI &rhs );
C3VecI operator* ( const std::complex<float> &other, const C3VecI &rhs );
C3Vec operator+ ( const C3Vec &other, const C3Vec &rhs);
C3VecI operator+ ( const C3VecI &other, const C3VecI &rhs);

std::ostream& operator<< ( std::ostream &os, const C3Vec &v );
std::ostream& operator<< ( std::ostream &os, const C3VecI &v );
std::vector<C3Vec> operator- ( const std::vector<C3Vec> &other, const C3Vec &rhs);
std::vector<C3Vec> operator+ ( const std::vector<C3Vec> &other, const C3Vec &rhs);
std::vector<C3VecI> operator+ ( const std::vector<C3VecI> &other, const C3VecI &rhs);
std::vector<C3Vec> operator- ( const std::vector<C3Vec> &other, const std::vector<C3Vec> &rhs);
std::vector<C3Vec> operator+ ( const std::vector<C3Vec> &other, const std::vector<C3Vec> &rhs);
std::vector<C3Vec> operator* ( const std::vector<C3Vec> &other, const std::vector<float> &rhs);
float mag ( const C3Vec &in );
float mag ( const C3VecI &in );
C3Vec pow ( const C3Vec &in, const int arg );
C3Vec sqrt ( const C3Vec &in );
float dot ( const C3Vec &Y, const C3Vec &X );
std::complex<float> dot ( const C3VecI &Y, const C3Vec &X );
C3Vec atan2 ( const C3Vec &Y, const C3Vec &X );
C3Vec cross ( const C3Vec A, const C3Vec B );
int isnan ( const C3Vec arg );
int isnan ( const C3VecI arg );
int isinf ( const C3Vec arg );
int isinf ( const C3VecI arg );

#include "c3vec.tpp" // since templated functions need to be in headers

#endif
