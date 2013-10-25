#include "constants.hpp"
#include "grid_sizes.hpp"
#include "cusp/complex.h"
#include "vecs.cuh"
#include <vector>

// used for complex
using cusp::complex;
using std::vector;

__host__ __device__ C3Vec::C3Vec () {
    c1=0;
    c2=0;
    c3=0;
}

__host__ __device__ C3Vec::C3Vec ( float _c1, float _c2, float _c3 )
{
    c1=_c1;
    c2=_c2;
    c3=_c3;
}

__host__ __device__ C3Vec::C3Vec ( int _arg ) 
{
    c1=_arg;
    c2=_arg;
    c3=_arg;
}

__host__ __device__ C3VecI::C3VecI ( complex<float> _c1, complex<float> _c2, complex<float> _c3 ) 
{
    c1=_c1;
    c2=_c2;
    c3=_c3;
}

__host__ __device__ C3VecI::C3VecI ()
{
    c1=complex<float>(0.0f,0.0f);
    c2=complex<float>(0.0f,0.0f);
    c3=complex<float>(0.0f,0.0f);
}

__host__ __device__
C3Vec& C3Vec::operator= (const C3Vec &rhs ) {
		if (this != &rhs) {
				c1 = rhs.c1;
				c2 = rhs.c2;
				c3 = rhs.c3;
		}
		return *this;
}

__host__ __device__ 
C3Vec& C3Vec::operator= (const float &rhs){
	c1 = rhs;
	c2 = rhs;
	c3 = rhs;
	return *this;
}

__host__ __device__
C3VecI& C3VecI::operator= (const C3VecI &rhs ) {
		if (this != &rhs) {
				c1 = rhs.c1;
				c2 = rhs.c2;
				c3 = rhs.c3;
		}
		return *this;
}
__host__ __device__
C3Vec& C3Vec::operator+= (const C3Vec &rhs ) {
		c1 = c1 + rhs.c1;
		c2 = c2 + rhs.c2;
		c3 = c3 + rhs.c3;
		return *this;
}

__host__ __device__
C3Vec& C3Vec::operator+= (const float &rhs ) {
		c1 = c1 + rhs;
		c2 = c2 + rhs;
		c3 = c3 + rhs;
		return *this;
}
__host__ __device__
C3Vec& C3Vec::operator-= (const C3Vec &rhs ) {
		c1 = c1 - rhs.c1;
		c2 = c2 - rhs.c2;
		c3 = c3 - rhs.c3;
		return *this;
}
__host__ __device__
C3Vec& C3Vec::operator-= (const float &rhs ) {
		c1 = c1 - rhs;
		c2 = c2 - rhs;
		c3 = c3 - rhs;
		return *this;
}
__host__ __device__
C3VecI& C3VecI::operator-= (const C3VecI &rhs ) {
		c1 = c1 - rhs.c1;
		c2 = c2 - rhs.c2;
		c3 = c3 - rhs.c3;
		return *this;
}
__host__ __device__
C3VecI& C3VecI::operator-= (const float &rhs ) {
		c1 = c1 - rhs;
		c2 = c2 - rhs;
		c3 = c3 - rhs;
		return *this;
}
__host__ __device__
C3Vec& C3Vec::operator*= (const C3Vec &rhs ) {
		c1 *= rhs.c1;
		c2 *= rhs.c2;
		c3 *= rhs.c3;
		return *this;
}
__host__ __device__
C3Vec& C3Vec::operator*= (const float &rhs ) {
		c1 *= rhs;
		c2 *= rhs;
		c3 *= rhs;
		return *this;
}
__host__ __device__
C3Vec& C3Vec::operator/= (const C3Vec &rhs ) {
		c1 /= rhs.c1;
		c2 /= rhs.c2;
		c3 /= rhs.c3;
		return *this;
}
__host__ __device__
C3Vec& C3Vec::operator/= (const float &rhs ) {
		c1 /= rhs;
		c2 /= rhs;
		c3 /= rhs;
		return *this;
}
__host__ __device__
C3VecI& C3VecI::operator/= (const C3VecI &rhs ) {
		c1 /= rhs.c1;
		c2 /= rhs.c2;
		c3 /= rhs.c3;
		return *this;
}
__host__ __device__
C3VecI& C3VecI::operator/= (const float &rhs ) {
		c1 /= rhs;
		c2 /= rhs;
		c3 /= rhs;
		return *this;
}
__host__ __device__
C3Vec C3Vec::operator+ (const C3Vec &other) {
		return C3Vec(this->c1+other.c1,this->c2+other.c2,this->c3+other.c3);
}
__host__ __device__
C3Vec C3Vec::operator+ (const float &other) {
		return C3Vec(*this)+=other;
}
__host__ __device__
C3Vec C3Vec::operator- (const C3Vec &other) {
		return C3Vec(*this)-=other;
}
__host__ __device__
C3Vec C3Vec::operator- (const float &other) {
		return C3Vec(*this)-=other;
}
__host__ __device__
C3VecI C3VecI::operator- (const C3VecI &other) {
		return C3VecI(*this)-=other;
}
__host__ __device__
C3VecI C3VecI::operator- (const float &other) {
		return C3VecI(*this)-=other;
}
__host__ __device__
C3Vec C3Vec::operator* (const C3Vec &other) {
		return C3Vec(*this)*=other;
}
__host__ __device__
C3Vec C3Vec::operator* (const float &other) {
		return C3Vec(*this)*=other;
}
__host__ __device__
C3Vec C3Vec::operator/ (const C3Vec &other) {
		return C3Vec(*this)/=other;
}
__host__ __device__
C3Vec C3Vec::operator/ (const float &other) {
		return C3Vec(*this)/=other;
}
__host__ __device__
C3VecI C3VecI::operator/ (const C3VecI &other) {
		return C3VecI(*this)/=other;
}
__host__ __device__
C3VecI C3VecI::operator/ (const float &other) {
		return C3VecI(*this)/=other;
}
// C3VecI 
__host__ __device__
C3VecI& C3VecI::operator+= (const C3VecI &rhs ) {
		c1 += rhs.c1;
		c2 += rhs.c2;
		c3 += rhs.c3;
		return *this;
}
__host__ __device__
C3VecI& C3VecI::operator+= (const float &rhs ) {
		c1 += rhs;
		c2 += rhs;
		c3 += rhs;
		return *this;
}
__host__ __device__
C3VecI& C3VecI::operator*= (const C3VecI &rhs ) {
		c1 *= rhs.c1;
		c2 *= rhs.c2;
		c3 *= rhs.c3;
		return *this;
}
__host__ __device__
C3VecI& C3VecI::operator*= (const float &rhs ) {
		c1 *= rhs;
		c2 *= rhs;
		c3 *= rhs;
		return *this;
}
__host__ __device__
C3VecI C3VecI::operator+ (const C3VecI &other) {
		return C3VecI(this->c1+other.c1,this->c2+other.c2,this->c3+other.c3);
}
__host__ __device__
C3VecI C3VecI::operator+ (const float &other) {
		return C3VecI(*this)+=other;
}
__host__ __device__
C3VecI C3VecI::operator* (const C3VecI &other) {
		return C3VecI(*this)*=other;
}
__host__ __device__
C3VecI C3VecI::operator* (const float &other) {
		return C3VecI(*this)*=other;
}

// Global (not member) functions for lhs operators
__host__ __device__
C3Vec operator* ( const float &other, const C3Vec &rhs ) {
		return C3Vec(rhs)*=other;
}
__host__ __device__
C3VecI operator* ( const float &other, const C3VecI &rhs ) {
		return C3VecI(rhs)*=other;
}
__host__ __device__
C3Vec operator+ ( const C3Vec &other, const C3Vec &rhs) {
		return C3Vec(other.c1+rhs.c1,other.c2+rhs.c2,other.c3+rhs.c3);
}
__host__ __device__
C3VecI operator+ ( const C3VecI &other, const C3VecI &rhs) {
		return C3VecI(other.c1+rhs.c1,other.c2+rhs.c2,other.c3+rhs.c3);
}
__host__ __device__
C3Vec pow ( const C3Vec &in, const int arg ) {
		C3Vec out;
		out.c1 = pow(in.c1,arg);
		out.c2 = pow(in.c2,arg);
		out.c3 = pow(in.c3,arg);
		return out;
}
__host__ __device__
C3Vec sqrt ( const C3Vec &in ) {
		C3Vec out;
		out.c1 = sqrt(in.c1);
		out.c2 = sqrt(in.c2);
		out.c3 = sqrt(in.c3);
		return out;
}
__host__ __device__
C3Vec atan2 ( const C3Vec &Y, const C3Vec &X ) {
		C3Vec out;
		out.c1 = atan2(Y.c1,X.c1);
		out.c2 = atan2(Y.c2,X.c2);
		out.c3 = atan2(Y.c3,X.c3);
		return out;
}

/*
__host__
float maxC3VecAbs ( const vector<C3Vec> &input ) {

        vector<float> inputAbs(input.size());
        for(int i=0;i<input.size();i++) {
                inputAbs[i] = sqrt(pow(input[i].c1,2)+pow(input[i].c2,2)+pow(input[i].c3,2));
        }
        return *max_element(inputAbs.begin(),inputAbs.end());
}
__host__
C3Vec intC3VecArray ( const vector<float> &x, const vector<C3Vec> &f ) {

        C3Vec result;
        float h = x[1]-x[0];
        for(int i=1;i<f.size();i++) {
                result += h/2.0*(f[i-1]+f[i]);
        }

        return result;
}
__host__
C3VecI intC3VecArray ( const vector<float> &x, const vector<C3VecI> &f ) {

        C3VecI result;
        float h = x[1]-x[0];
        for(int i=1;i<f.size();i++) {
                result += h/2.0*(f[i-1]+f[i]);
        }

        return result;
}
__host__
C3VecI intC3VecArray ( const float x[], const vector<C3VecI> &f ) {

        C3VecI result;
        float h = x[1]-x[0];
        for(int i=1;i<f.size();i++) {
                result += h/2.0*(f[i-1]+f[i]);
        }

        return result;
}
__host__
C3Vec intC3VecArray ( const float x[], const vector<C3Vec> &f ) {

        C3Vec result;
        float h = x[1]-x[0];
        for(int i=1;i<f.size();i++) {
                result += h/2.0*(f[i-1]+f[i]);
        }

        return result;
}
*/
