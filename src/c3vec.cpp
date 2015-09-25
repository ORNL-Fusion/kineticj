#include "c3vec.hpp"

C3Vec& C3Vec::operator= (const C3Vec &rhs ) {
		if (this != &rhs) {
				c1 = rhs.c1;
				c2 = rhs.c2;
				c3 = rhs.c3;
		}
		return *this;
}
C3VecI& C3VecI::operator= (const C3VecI &rhs ) {
		if (this != &rhs) {
				c1 = rhs.c1;
				c2 = rhs.c2;
				c3 = rhs.c3;
		}
		return *this;
}
C3Vec& C3Vec::operator+= (const C3Vec &rhs ) {
		c1 += rhs.c1;
		c2 += rhs.c2;
		c3 += rhs.c3;
		return *this;
}

C3Vec& C3Vec::operator+= (const float &rhs ) {
		c1 += rhs;
		c2 += rhs;
		c3 += rhs;
		return *this;
}

C3Vec& C3Vec::operator-= (const C3Vec &rhs ) {
		c1 -= rhs.c1;
		c2 -= rhs.c2;
		c3 -= rhs.c3;
		return *this;
}

C3Vec& C3Vec::operator-= (const float &rhs ) {
		c1 -= rhs;
		c2 -= rhs;
		c3 -= rhs;
		return *this;
}
C3VecI& C3VecI::operator-= (const C3VecI &rhs ) {
		c1 -= rhs.c1;
		c2 -= rhs.c2;
		c3 -= rhs.c3;
		return *this;
}

C3VecI& C3VecI::operator-= (const float &rhs ) {
		c1 -= rhs;
		c2 -= rhs;
		c3 -= rhs;
		return *this;
}

C3Vec& C3Vec::operator*= (const C3Vec &rhs ) {
		c1 *= rhs.c1;
		c2 *= rhs.c2;
		c3 *= rhs.c3;
		return *this;
}

C3Vec& C3Vec::operator*= (const float &rhs ) {
		c1 *= rhs;
		c2 *= rhs;
		c3 *= rhs;
		return *this;
}

C3Vec& C3Vec::operator/= (const C3Vec &rhs ) {
		c1 /= rhs.c1;
		c2 /= rhs.c2;
		c3 /= rhs.c3;
		return *this;
}

C3Vec& C3Vec::operator/= (const float &rhs ) {
		c1 /= rhs;
		c2 /= rhs;
		c3 /= rhs;
		return *this;
}

C3VecI& C3VecI::operator/= (const C3VecI &rhs ) {
		c1 /= rhs.c1;
		c2 /= rhs.c2;
		c3 /= rhs.c3;
		return *this;
}

C3VecI& C3VecI::operator/= (const float &rhs ) {
		c1 /= rhs;
		c2 /= rhs;
		c3 /= rhs;
		return *this;
}
C3Vec C3Vec::operator+ (const C3Vec &other) {
		return C3Vec(this->c1+other.c1,this->c2+other.c2,this->c3+other.c3);
}

C3Vec C3Vec::operator+ (const float &other) {
		return C3Vec(*this)+=other;
}

C3Vec C3Vec::operator- (const C3Vec &other) {
		return C3Vec(*this)-=other;
}

C3Vec C3Vec::operator- (const float &other) {
		return C3Vec(*this)-=other;
}

C3VecI C3VecI::operator- (const C3VecI &other) {
		return C3VecI(*this)-=other;
}

C3VecI C3VecI::operator- (const float &other) {
		return C3VecI(*this)-=other;
}

C3Vec C3Vec::operator* (const C3Vec &other) {
		return C3Vec(*this)*=other;
}

C3Vec C3Vec::operator* (const float &other) {
		return C3Vec(*this)*=other;
}

C3Vec C3Vec::operator/ (const C3Vec &other) {
		return C3Vec(*this)/=other;
}

C3Vec C3Vec::operator/ (const float &other) {
		return C3Vec(*this)/=other;
}

C3VecI C3VecI::operator/ (const C3VecI &other) {
		return C3VecI(*this)/=other;
}

C3VecI C3VecI::operator/ (const float &other) {
		return C3VecI(*this)/=other;
}
// C3VecI 

C3VecI& C3VecI::operator+= (const C3VecI &rhs ) {
		c1 += rhs.c1;
		c2 += rhs.c2;
		c3 += rhs.c3;
		return *this;
}

C3VecI& C3VecI::operator+= (const float &rhs ) {
		c1 += rhs;
		c2 += rhs;
		c3 += rhs;
		return *this;
}

C3VecI& C3VecI::operator*= (const C3VecI &rhs ) {
		c1 *= rhs.c1;
		c2 *= rhs.c2;
		c3 *= rhs.c3;
		return *this;
}

C3VecI& C3VecI::operator*= (const float &rhs ) {
		c1 *= rhs;
		c2 *= rhs;
		c3 *= rhs;
		return *this;
}
C3VecI C3VecI::operator+ (const C3VecI &other) {
		return C3VecI(this->c1+other.c1,this->c2+other.c2,this->c3+other.c3);
}

C3VecI C3VecI::operator+ (const float &other) {
		return C3VecI(*this)+=other;
}
C3VecI C3VecI::operator* (const C3VecI &other) {
		return C3VecI(*this)*=other;
}

C3VecI C3VecI::operator* (const float &other) {
		return C3VecI(*this)*=other;
}

// Global (not member) functions for lhs operators

C3Vec operator* ( const float &other, const C3Vec &rhs ) {
		return C3Vec(rhs)*=other;
}

C3VecI operator* ( const float &other, const C3VecI &rhs ) {
		return C3VecI(rhs)*=other;
}

C3VecI operator* ( const std::complex<float> &other, const C3VecI &rhs ) {
        C3VecI tmp;
        tmp.c1 = other * rhs.c1;
        tmp.c2 = other * rhs.c2;
        tmp.c3 = other * rhs.c3;
		return tmp;
}

C3Vec operator+ ( const C3Vec &other, const C3Vec &rhs) {
		return C3Vec(other.c1+rhs.c1,other.c2+rhs.c2,other.c3+rhs.c3);
}
C3VecI operator+ ( const C3VecI &other, const C3VecI &rhs) {
		return C3VecI(other.c1+rhs.c1,other.c2+rhs.c2,other.c3+rhs.c3);
}


