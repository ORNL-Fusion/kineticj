// Current hacks:
// the E-field is being multiplied by some factor (e.g. 1e3) in rk4_evalf

#include <string>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <cstdlib>
#include <netcdf>
#include <vector>
#include <algorithm>
#include <complex>
#include "constants.hpp"
#include <libconfig.h++>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <assert.h>
#include <omp.h>

#if CLOCK >= 1 
#include <ctime>
#endif

#if USEPAPI >= 1
#include <papi.h>
#endif

#if LOWMEM_USEPAPI >= 1
#include <papi.h>
#endif


//#include <google/profiler.h>

#ifdef __CUDA_ARCH__
#define PRINT cuPrintf 
#else
#define PRINT printf
#endif

using namespace std;
using namespace netCDF;
using namespace exceptions;

class CSpecies {
		public:
				double m, q;
				double amu;
				int Z;
				string name;

				CSpecies () {};
				CSpecies ( const double amu, const int Z);
				CSpecies ( const double amu, const int Z, const char *_s);
};

CSpecies::CSpecies ( const double _amu, const int _Z ) {

		amu = _amu;
		Z = _Z;
		m = amu * _mi;
		q = Z * _e;
		name = " ";
}

CSpecies::CSpecies ( const double _amu, const int _Z, const char *_s ) {
		amu = _amu;
		Z = _Z;
		m = amu * _mi;
		q = Z * _e;
		name = string(_s);
}

class CParticle: public CSpecies {
		public:
				float c1, c2, c3, v_c1, v_c2, v_c3;
				int number;
				float weight;
				int status;
                //float df0_dvx, df0_dvy, df0_dvz;
                float dvx, dvy, dvz, d3v;
                float vPar, vPer, gyroPhase, u, vTh;
                float vAlp, vBet, phs, rfPhase;

				CParticle ();
				CParticle ( double _amu, int _Z);
				CParticle (float c1, float c2, float c3, 
								float v_c1, float v_c2, float v_c3, 
								double _amu, int _Z, float _weight, float _rfPhase);
    
				CParticle (CSpecies _species);
};

CParticle::CParticle () {
		status = 0;
}

CParticle::CParticle ( double _amu, int _Z ): CSpecies(_amu,_Z) {
		status = 0;
}

CParticle::CParticle 
	(float _c1, float _c2, float _c3, 
	 float _v_c1, float _v_c2, float _v_c3, 
	 double _amu, int _Z, float _weight, float _rfPhase) :
	CSpecies (_amu, _Z) {

	c1 = _c1;
	c2 = _c2;
	c3 = _c3;

	v_c1 = _v_c1;
	v_c2 = _v_c2;
	v_c3 = _v_c3;

	weight = _weight;
    
    rfPhase = _rfPhase;
	
	status = 0;
}

CParticle::CParticle ( CSpecies _species ) : CSpecies(_species) {
}

class C3VecI {
		public:
				complex<float> c1, c2, c3;


                C3VecI (int _const) {c1=_const;c2=_const;c3=_const;};
				C3VecI () {c1=complex<float>(0.0f,0.0f);c2=complex<float>(0.0f,0.0f);c3=complex<float>(0.0f,0.0f);};
				C3VecI ( complex<float> _c1, complex<float> _c2, complex<float> _c3 ) {c1=_c1;c2=_c2;c3=_c3;};

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

C3VecI operator* ( const complex<float> &other, const C3VecI &rhs ) {
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

ostream& operator<< ( ostream &os, const C3Vec &v ) {
        os << v.c1 << ", "<< v.c2 << ", " << v.c3;
        return os;
}

ostream& operator<< ( ostream &os, const C3VecI &v ) {
        os << v.c1 << ", "<< v.c2 << ", " << v.c3;
        return os;
}


vector<C3Vec> operator- ( const vector<C3Vec> &other, const C3Vec &rhs) {
		vector<C3Vec> out(other.size());
		for(int i=0;i<other.size();i++) {
				out[i].c1 = other[i].c1 - rhs.c1;
				out[i].c2 = other[i].c2 - rhs.c2;
				out[i].c3 = other[i].c3 - rhs.c3;
		}
		return out;
}

vector<C3Vec> operator+ ( const vector<C3Vec> &other, const C3Vec &rhs) {
		vector<C3Vec> out(other.size());
		for(int i=0;i<other.size();i++) {
				out[i].c1 = other[i].c1 + rhs.c1;
				out[i].c2 = other[i].c2 + rhs.c2;
				out[i].c3 = other[i].c3 + rhs.c3;
		}
		return out;
}
vector<C3VecI> operator+ ( const vector<C3VecI> &other, const C3VecI &rhs) {
		vector<C3VecI> out(other.size());
		for(int i=0;i<other.size();i++) {
				out[i].c1 = other[i].c1 + rhs.c1;
				out[i].c2 = other[i].c2 + rhs.c2;
				out[i].c3 = other[i].c3 + rhs.c3;
		}
		return out;
}
vector<C3Vec> operator- ( const vector<C3Vec> &other, const vector<C3Vec> &rhs) {
		assert(other.size()==rhs.size());
		vector<C3Vec> out(other.size());
		for(int i=0;i<other.size();i++) {
				out[i].c1 = other[i].c1 - rhs[i].c1;
				out[i].c2 = other[i].c2 - rhs[i].c2;
				out[i].c3 = other[i].c3 - rhs[i].c3;
		}
		return out;
}

vector<C3Vec> operator+ ( const vector<C3Vec> &other, const vector<C3Vec> &rhs) {
		assert(other.size()==rhs.size());
		vector<C3Vec> out(other.size());
		for(int i=0;i<other.size();i++) {
				out[i].c1 = other[i].c1 + rhs[i].c1;
				out[i].c2 = other[i].c2 + rhs[i].c2;
				out[i].c3 = other[i].c3 + rhs[i].c3;
		}
		return out;
}

vector<C3Vec> operator* ( const vector<C3Vec> &other, const vector<float> &rhs) {
		assert(other.size()==rhs.size());
		vector<C3Vec> out(other.size());
		for(int i=0;i<other.size();i++) {
				out[i].c1 = other[i].c1 * rhs[i];
				out[i].c2 = other[i].c2 * rhs[i];
				out[i].c3 = other[i].c3 * rhs[i];
		}
		return out;
}

float mag ( const C3Vec &in ) {
		return sqrt(pow(in.c1,2)+pow(in.c2,2)+pow(in.c3,2));
}

C3Vec pow ( const C3Vec &in, const int arg ) {
		C3Vec out;
		out.c1 = pow(in.c1,arg);
		out.c2 = pow(in.c2,arg);
		out.c3 = pow(in.c3,arg);
		return out;
}

C3Vec sqrt ( const C3Vec &in ) {
		C3Vec out;
		out.c1 = sqrt(in.c1);
		out.c2 = sqrt(in.c2);
		out.c3 = sqrt(in.c3);
		return out;
}

float dot ( const C3Vec &Y, const C3Vec &X ) {
		return Y.c1*X.c1 + Y.c2*X.c2 + Y.c3*X.c3;
}

complex<float> dot ( const C3VecI &Y, const C3Vec &X ) {
		return Y.c1*X.c1 + Y.c2*X.c2 + Y.c3*X.c3;
}

C3Vec atan2 ( const C3Vec &Y, const C3Vec &X ) {
		C3Vec out;
		out.c1 = atan2(Y.c1,X.c1);
		out.c2 = atan2(Y.c2,X.c2);
		out.c3 = atan2(Y.c3,X.c3);
		return out;
}

C3VecI cross ( const C3VecI A, const C3Vec B ) {

        C3VecI answer;
        answer.c1 =  (A.c2*B.c3 - A.c3*B.c2);
        answer.c2 = -(A.c1*B.c3 - A.c3*B.c1);
        answer.c3 =  (A.c1*B.c2 - A.c2*B.c1);
        return answer;
}

C3VecI cross ( const C3VecI A, const C3VecI B ) {

        C3VecI answer;
        answer.c1 =  (A.c2*B.c3 - A.c3*B.c2);
        answer.c2 = -(A.c1*B.c3 - A.c3*B.c1);
        answer.c3 =  (A.c1*B.c2 - A.c2*B.c1);
        return answer;
}

C3VecI cross ( const C3Vec A, const C3VecI B ) {

        C3VecI answer;
        answer.c1 =  (A.c2*B.c3 - A.c3*B.c2);
        answer.c2 = -(A.c1*B.c3 - A.c3*B.c1);
        answer.c3 =  (A.c1*B.c2 - A.c2*B.c1);
        return answer;
}

C3Vec cross ( const C3Vec A, const C3Vec B ) {

        C3Vec answer;
        answer.c1 =  (A.c2*B.c3 - A.c3*B.c2);
        answer.c2 = -(A.c1*B.c3 - A.c3*B.c1);
        answer.c3 =  (A.c1*B.c2 - A.c2*B.c1);
        return answer;
}

float maxC3VecAbs ( const vector<C3Vec> &input ) {

	vector<float> inputAbs(input.size());
	for(int i=0;i<input.size();i++) {
		inputAbs[i] = sqrt(pow(input[i].c1,2)+pow(input[i].c2,2)+pow(input[i].c3,2));
	}
	return *max_element(inputAbs.begin(),inputAbs.end());
}

complex<float> intVecArray ( const vector<float> &x, const vector<complex<float> > &f ) {

	complex<float> result;
	float h = x[1]-x[0];
	for(int i=1;i<f.size();i++) {
		result += h/2.0f*(f[i-1]+f[i]);
	}

	return result;
}

C3Vec intVecArray ( const vector<float> &x, const vector<C3Vec> &f ) {

	C3Vec result;
	float h = x[1]-x[0];
	for(int i=1;i<f.size();i++) {
		result += h/2.0*(f[i-1]+f[i]);
	}

	return result;
}

C3VecI intVecArray ( const vector<float> &x, const vector<C3VecI> &f ) {

	C3VecI result;
	float h = x[1]-x[0];
	for(int i=1;i<f.size();i++) {
		result += h/2.0*(f[i-1]+f[i]);
	}

	return result;
}

void kj_print ( const C3Vec arg, string name ) {
    cout << name<<".c1: " << arg.c1 << "  "<< name<<".c2: " << arg.c2 << "  "<<name<<".c3: " << arg.c3 << endl;
    return;
}

void kj_print ( const float arg, string name ) {
    cout << name <<": "<< arg << endl;
    return;
}

int isnan ( const C3Vec arg ) {
    int answer = 0;
    if(isnan(arg.c1)) answer = 1;
    if(isnan(arg.c2)) answer = 1;
    if(isnan(arg.c3)) answer = 1;
    return answer;
}

int isinf ( const C3Vec arg ) {
    int answer = 0;
    if(isinf(arg.c1)) answer = 1;
    if(isinf(arg.c2)) answer = 1;
    if(isinf(arg.c3)) answer = 1;
    return answer;
}

C3Vec XYZ_to_CYL ( const C3Vec xyz ) {
        C3Vec cyl;
        cyl.c1 = sqrt(pow(xyz.c1,2)+pow(xyz.c2,2));
        cyl.c2 = atan2(xyz.c2,xyz.c1);
        cyl.c3 = xyz.c3;
        return cyl;
}

C3Vec CYL_to_XYZ ( const C3Vec cyl ) {
        C3Vec xyz;
        xyz.c1 = cyl.c1*cos(cyl.c2);
        xyz.c2 = cyl.c1*sin(cyl.c2);
        xyz.c3 = cyl.c3;
        return xyz;
}

class fieldMeshClass {
    public:
        vector<float> r;
        vector<float> z;
    
        fieldMeshClass ( vector<float> r, vector<float> z);
        fieldMeshClass ( vector<float> r);
};

fieldMeshClass::fieldMeshClass ( vector<float> _r, vector<float> _z) {
		r = _r;
		z = _z;
}

fieldMeshClass::fieldMeshClass ( vector<float> _r) {
		r = _r;
}

////1D w/o a Cparticle input
template<class TYPE>
TYPE kj_interp ( const C3Vec &Loc, const fieldMeshClass &fieldMesh, const vector<TYPE> &fVec, int &status ) {

    status = 0;
    
    float x = Loc.c1;
    
	float _x, x0, x1;
	float xTmp;
    xTmp = x;

#if _PARTICLE_BOUNDARY == 1
	if(x<fieldMesh.r.front()||x>fieldMesh.r.back()) {
			// Particle absorbing walls
#if DEBUG_INTERP >= 2
			cout<<"Particle absorbed at "<<x<<endl;
            cout<<"x:"<<x<<endl;
            cout<<"xVec.front():"<<fieldMesh.r.front()<<endl;
            cout<<"xVec.back():"<<fieldMesh.r.back()<<endl;
#endif
			status = 1;
			return TYPE(0);
	}
#elif _PARTICLE_BOUNDARY == 2
			// Periodic
            float xRange = fieldMesh.r.back() - fieldMesh.r.front();
            int N;
			if(xTmp<fieldMesh.r.front()){
                N = int(floor( (fieldMesh.r.front() - xTmp)/xRange )) + 1;
                xTmp = xTmp + N*xRange;
               // xTmp = xVec.back()-(xVec.front()-xTmp);
             }
			if(xTmp>fieldMesh.r.back()){
                N = int(floor( (xTmp - fieldMesh.r.back())/xRange )) + 1;
                xTmp = xTmp - N*xRange;
            }
    
#elif _PARTICLE_BOUNDARY == 3
			// Particle reflecting walls
			if(xTmp<fieldMesh.r.front()) xTmp = fieldMesh.r.front()+(fieldMesh.r.front()-xTmp);
			if(xTmp>fieldMesh.r.back()) xTmp = fieldMesh.r.back()-(xTmp-fieldMesh.r.back());
#endif

#if DEBUG_INTERP >= 1    
	if(status>0){
#if DEBUG_INTERP >= 2
			cout<<"ERROR: Should never get here with _PARTICLE_BOUNDARY ==2|3"<<endl;
#endif
			return TYPE(0);
	}
#endif
	//else
	//{
		_x = (xTmp-fieldMesh.r.front())/(fieldMesh.r.back()-fieldMesh.r.front())*(fieldMesh.r.size()-1);
	//}

	x0 = floor(_x);
	x1 = ceil(_x);
	
	// Catch for particle at point
	if(x0==x1) {
#if DEBUG_INTERP >= 2
        cout << "status version of kj_interp" << endl;
		cout << "x0: " << x0 << " x1: " <<x1<< " _x: "<<_x << endl;
		cout << "Particle at point catch: " << x0/x1 << "  "  << abs(1.0-x0/x1) << endl;
#endif
		return fVec[x0];
	}
	else {
		TYPE f0 = fVec[x0];
		TYPE f1 = fVec[x1];

		return f0+(_x-x0)*(f1-f0)/(x1-x0);
	}
}

////1D w/ a Cparticle input
template<class TYPE>
TYPE kj_interp ( const C3Vec &Loc, const fieldMeshClass &fieldMesh, const vector<TYPE> &fVec, CParticle &p, int &status ) {

    status = 0;
    float x = Loc.c1;
	float _x, x0, x1;
	float xTmp;
	xTmp = x;

#if _PARTICLE_BOUNDARY == 1
	if(x<fieldMesh.r.front()||x>fieldMesh.r.back()) {
			// Particle absorbing walls
#if DEBUG_INTERP >= 1
            if(fieldMesh.r.size()!=fVec.size()) {
#if DEBUG_INTERP >= 2
                cout<<"ERROR: xVec and fVec are not the same size for interpolation!"<<endl;
#endif
                p.status = 1;
                status = 1;
                return TYPE(0);
            }
#if DEBUG_INTERP >= 2
			cout<<"Particle absorbed at "<<x<<endl;
            cout<<"Particle number: "<<p.number<<endl;
            cout<<"x:"<<x<<endl;
            cout<<"xVec.front():"<<xVec.front()<<endl;
            cout<<"xVec.back():"<<xVec.back()<<endl;
#endif
#endif
            status = 1;
			p.status = 1;
			return TYPE(0);
	}
#elif _PARTICLE_BOUNDARY == 2
			// Periodic 
            float xRange = fieldMesh.r.back() - fieldMesh.r.front();
            int N;
			if(xTmp<fieldMesh.r.front()){
                N = int(floor( (fieldMesh.r.front() - xTmp)/xRange )) + 1;
                xTmp = xTmp + N*xRange;
               // xTmp = xVec.back()-(xVec.front()-xTmp);
             }
			if(xTmp>fieldMesh.r.back()){
                N = int(floor( (xTmp - fieldMesh.r.back())/xRange )) + 1;
                xTmp = xTmp - N*xRange;
            }
#elif _PARTICLE_BOUNDARY == 3
			// Particle reflecting walls
			if(xTmp<fieldMesh.r.front()) xTmp = fieldMesh.r.front()+(fieldMesh.r.front()-xTmp);
			if(xTmp>fieldMesh.r.back()) xTmp = fieldMesh.r.back()-(xTmp-fieldMesh.r.back());
#endif

#if DEBUG_INTERP >= 1    
	if(p.status>0){
#if DEBUG_INTERP >= 2
			cout<<"ERROR: Should never get here with _PARTICLE_BOUNDARY ==2|3"<<endl;
#endif
            status = 1;
			return TYPE(0);
	}
#endif
	//else
	//{
		_x = (xTmp-fieldMesh.r.front())/(fieldMesh.r.back()-fieldMesh.r.front())*(fieldMesh.r.size()-1);
	//}

	x0 = floor(_x);
	x1 = ceil(_x);
	
	// Catch for particle at point
	if(x0==x1) {
#if DEBUG_INTERP >= 2
        cout << "Particle version of kj_interp" << endl;
		cout << "x0: " << x0 << " x1: " <<x1<< " _x: "<<_x << endl;
		cout << "Particle at point catch: " << x0/x1 << "  "  << abs(1.0-x0/x1) << endl;
#endif
		return fVec[x0];
	}
	else {
		TYPE f0 = fVec[x0];
		TYPE f1 = fVec[x1];
#if DEBUG_INTERP >=2
        //cout << "kj_interp: " << endl;
        //if(typeid(TYPE)==typeid(C3Vec)) cout << "Type is C3Vec" << endl;
        //if(typeid(TYPE)==typeid(float)) cout << "Type is float" << endl;
        //kj_print(x0,"x0");
        //kj_print(x1,"x1");
        //kj_print(y0,"y0");
        //kj_print(y1,"y1");
        //cout << endl;
        if(x0>fVec.size()-1||x0<0||x1>fVec.size()-1||x1<1) {
                cout<<"ERROR: interpolation point off the end of array"<<endl;
                cout<<"x.front: "<<fieldMesh.r.front()<<endl;
                cout<<"x.back: "<<fieldMesh.r.back()<<endl;
                cout<<"x: "<<x<<endl;
                cout<<"xVec.size(): "<<fieldMesh.r.size()<<endl;
                cout<<"fVec.size(): "<<fVec.size()<<endl;
                ++p.status;
                status = 1;
                return TYPE(0);
        }
#endif
        TYPE result = f0+(_x-x0)*(f1-f0)/(x1-x0);

#if DEBUG_INTERP >=1
        if(isnan(result)) {
#if DEBUG_INTERP >= 2
                cout<<"ERROR: interpolation produced a NaN"<<endl;
#endif
                ++p.status;
                status = 1;
                return TYPE(0);
        } 
        if(isinf(result)) {
#if DEBUG_INTERP >= 2
                cout<<"ERROR: interpolation produced a INF"<<endl;
#endif
                ++p.status;
                status = 1;
                return TYPE(0);
        }
#endif
		return result;
	}

}

/////////////////////////////////////////////////////////// CHANGE x TO r FOR 2D INTERPOLATION, OR USE MORE GENERAL c1, c2, c3 ???
//// 2D templated interp w/ Cparticle input
template<class TYPE2>
TYPE2 kj_interp ( const C3Vec &Loc, const fieldMeshClass &fieldMesh, const vector< vector< TYPE2 > > &fVec, CParticle &p, int &status ) {

    status = 0;
    float x = Loc.c1;
    float z = Loc.c3;
    
  //  cout << "x   " << x << endl;
  //  cout << "z   " << z << endl;
    
	float _x, x0, x1, _z, z0, z1;
    
    int ix, iz;
    
	float xTmp = x;
	float zTmp = z;

#if _PARTICLE_BOUNDARY == 1
    
	if(x < fieldMesh.r.front()||x>fieldMesh.r.back() || z < fieldMesh.z.front() || z > fieldMesh.z.back() ) {
			// Particle absorbing walls
#if DEBUG_INTERP >= 1
            if(fieldMesh.r.size()!=fVec.size()) {
#if DEBUG_INTERP >= 2
                cout<<"ERROR: xVec and fVec are not the same size for interpolation!"<<endl;
#endif
                p.status = 1;
                status = 1;
                return TYPE2(0);
            }
#if DEBUG_INTERP >= 2
			cout<<"Particle absorbed at "<<x<<endl;
            cout<<"Particle number: "<<p.number<<endl;
            cout<<"x:"<<x<<endl;
            cout<<"xVec.front():"<<fieldMesh.r.front()<<endl;
            cout<<"xVec.back():"<<fieldMesh.r.back()<<endl;
#endif
#endif
            status = 1;
			p.status = 1;
			return TYPE2(0);
	}
    
#elif _PARTICLE_BOUNDARY == 2
			// Periodic 
            float xRange = fieldMesh.r.back() - fieldMesh.r.front();
            float zRange = fieldMesh.z.back() - fieldMesh.z.front();

            int N;
			if(xTmp<fieldMesh.r.front()){
                N = int(floor( (fieldMesh.r.front() - xTmp)/xRange )) + 1;
                xTmp = xTmp + N*xRange;
               // xTmp = xVec.back()-(xVec.front()-xTmp);
             }
			if(xTmp>fieldMesh.r.back()){
                N = int(floor( (xTmp - fieldMesh.r.back())/xRange )) + 1;
                xTmp = xTmp - N*xRange;
            }
    
			if(zTmp<fieldMesh.z.front()){
                N = int(floor( (fieldMesh.z.front() - zTmp)/zRange )) + 1;
                zTmp = zTmp + N*zRange;
               // xTmp = xVec.back()-(xVec.front()-xTmp);
             }
			if(zTmp>fieldMesh.z.back()){
                N = int(floor( (zTmp - fieldMesh.z.back())/zRange )) + 1;
                zTmp = zTmp - N*zRange;
            }
    
#elif _PARTICLE_BOUNDARY == 3
			// Particle reflecting walls
			if(xTmp<fieldMesh.r.front()) xTmp = fieldMesh.r.front()+(fieldMesh.r.front()-xTmp);
			if(xTmp>fieldMesh.r.back()) xTmp = fieldMesh.r.back()-(xTmp-fieldMesh.r.back());
    
			if(zTmp<fieldMesh.z.front()) zTmp = fieldMesh.z.front()+(fieldMesh.z.front()-zTmp);
			if(zTmp>fieldMesh.z.back()) zTmp = fieldMesh.z.back()-(zTmp-fieldMesh.z.back());
#endif

#if DEBUG_INTERP >= 1    
	if(p.status>0){
#if DEBUG_INTERP >= 2
			cout<<"ERROR: Should never get here with _PARTICLE_BOUNDARY ==2|3"<<endl;
#endif
            status = 1;
			return TYPE2(0);
	}
#endif
	//else
	//{
		_x = (xTmp-fieldMesh.r.front())/(fieldMesh.r.back()-fieldMesh.r.front())*(fieldMesh.r.size()-1);
		_z = (zTmp-fieldMesh.z.front())/(fieldMesh.z.back()-fieldMesh.z.front())*(fieldMesh.z.size()-1);
	//}

	//x0 = floor(_x);
	ix = (int) floor(_x);
    //x1 = ceil(_x);

	//z0 = floor(_z);
	iz = (int) floor(_z);
	//z1 = ceil(_z);
	// Catch for particle at point
//	if(x0==x1) {
//#if DEBUG_INTERP >= 2
//        cout << "Particle version of kj_interp" << endl;
//		cout << "x0: " << x0 << " x1: " <<x1<< " _x: "<<_x << endl;
//		cout << "Particle at point catch: " << x0/x1 << "  "  << abs(1.0-x0/x1) << endl;
//#endif
//		return fVec[x0][z0];
//	}
//	else {

		TYPE2 f00 = fVec[ix][iz];
		TYPE2 f01 = fVec[ix][iz + 1];
		TYPE2 f10 = fVec[ix + 1][iz];
		TYPE2 f11 = fVec[ix + 1][iz + 1];
        
#if DEBUG_INTERP >=2
        //cout << "kj_interp: " << endl;
        //if(typeid(TYPE)==typeid(C3Vec)) cout << "Type is C3Vec" << endl;
        //if(typeid(TYPE)==typeid(float)) cout << "Type is float" << endl;
        //kj_print(x0,"x0");
        //kj_print(x1,"x1");
        //kj_print(y0,"y0");
        //kj_print(y1,"y1");
        //cout << endl;
        if(x0>fVec.size()-1||x0<0||x1>fVec.size()-1||x1<1) {
                cout<<"ERROR: interpolation point off the end of array"<<endl;
                cout<<"x.front: "<<fieldMesh.r.front()<<endl;
                cout<<"x.back: "<<fieldMesh.r.back()<<endl;
                cout<<"x: "<<x<<endl;
                cout<<"xVec.size(): "<<fieldMesh.r.size()<<endl;
                cout<<"fVec.size(): "<<fVec.size()<<endl;
                ++p.status;
                status = 1;
                return TYPE2(0);
        }
#endif
        x0 = fieldMesh.r[ix];
        x1 = fieldMesh.r[ix + 1];

        z0 = fieldMesh.z[iz];
        z1 = fieldMesh.z[iz + 1];
        
       // TYPE2 result = (1.0/( (x1 - x0)*(z1 - z0) ))*(f00*(x1 - _x)*(z1 - _z) + f01*(_x - x0)*(z1 - _z) + f10*(x1 - _x)*(_z - z0) + f11*(_x - x0)*(z - z0) );
        TYPE2 result = (1.0/( (x1 - x0)*(z1 - z0) ))*(f00*(x1 - x)*(z1 - z) + f01*(x - x0)*(z1 - z) + f10*(x1 - x)*(z - z0) + f11*(x - x0)*(z - z0) );
/*
        cout << "x0   " << x0 << endl;
        cout << "x1   " << x1 << endl;
        cout << "x   " << x << endl;

        cout << "z0   " << z0 << endl;
        cout << "z1   " << z1 << endl;
        cout << "z   " << z << endl;

        cout << "f00   " << f00 << endl;
        cout << "f01   " << f01 << endl;
        cout << "f10   " << f10 << endl;
        cout << "f11   " << f11 << endl;
        cout << "result   " << result << endl;
*/ 
    

#if DEBUG_INTERP >=1
        if(isnan(result)) {
#if DEBUG_INTERP >= 2
                cout<<"ERROR: interpolation produced a NaN"<<endl;
#endif
                ++p.status;
                status = 1;
                return TYPE2(0);
        } 
        if(isinf(result)) {
#if DEBUG_INTERP >= 2
                cout<<"ERROR: interpolation produced a INF"<<endl;
#endif
                ++p.status;
                status = 1;
                return TYPE2(0);
        }
#endif
		return result;
	//}
}


//// 2D templated interp w/o Cparticle input
template<class TYPE2>
TYPE2 kj_interp ( const C3Vec &Loc, const fieldMeshClass &fieldMesh, const vector< vector< TYPE2 > > &fVec, int &status ) {

    status = 0;
    float x = Loc.c1;
    float z = Loc.c3;
    
	float _x, x0, x1, _z, z0, z1;
	float xTmp = x;
	float zTmp = z;
    
    int ix, iz;

#if _PARTICLE_BOUNDARY == 1
	if(x < fieldMesh.r.front()||x>fieldMesh.r.back() || z < fieldMesh.z.front() || z > fieldMesh.z.back() ) {
			// Particle absorbing walls
#if DEBUG_INTERP >= 1
            if(fieldMesh.r.size()!=fVec.size()) {
#if DEBUG_INTERP >= 2
                cout<<"ERROR: xVec and fVec are not the same size for interpolation!"<<endl;
#endif
                status = 1;
                return TYPE2(0);
            }
#if DEBUG_INTERP >= 2
			cout<<"Particle absorbed at "<<x<<endl;
            cout<<"x:"<<x<<endl;
            cout<<"xVec.front():"<<fieldMesh.r.front()<<endl;
            cout<<"xVec.back():"<<fieldMesh.r.back()<<endl;
#endif
#endif
            status = 1;
			return TYPE2(0);
	}
#elif _PARTICLE_BOUNDARY == 2
			// Periodic 
            float xRange = fieldMesh.r.back() - fieldMesh.r.front();
            float zRange = fieldMesh.z.back() - fieldMesh.z.front();

            int N;
			if(xTmp<fieldMesh.r.front()){
                N = int(floor( (fieldMesh.r.front() - xTmp)/xRange )) + 1;
                xTmp = xTmp + N*xRange;
               // xTmp = xVec.back()-(xVec.front()-xTmp);
             }
			if(xTmp>fieldMesh.r.back()){
                N = int(floor( (xTmp - fieldMesh.r.back())/xRange )) + 1;
                xTmp = xTmp - N*xRange;
            }
    
			if(zTmp<fieldMesh.z.front()){
                N = int(floor( (fieldMesh.z.front() - zTmp)/zRange )) + 1;
                zTmp = zTmp + N*zRange;
               // xTmp = xVec.back()-(xVec.front()-xTmp);
             }
			if(zTmp>fieldMesh.z.back()){
                N = int(floor( (zTmp - fieldMesh.z.back())/zRange )) + 1;
                zTmp = zTmp - N*zRange;
            }
    
#elif _PARTICLE_BOUNDARY == 3
			// Particle reflecting walls
			if(xTmp<fieldMesh.r.front()) xTmp = fieldMesh.r.front()+(fieldMesh.r.front()-xTmp);
			if(xTmp>fieldMesh.r.back()) xTmp = fieldMesh.r.back()-(xTmp-fieldMesh.r.back());
    
			if(zTmp<fieldMesh.z.front()) zTmp = fieldMesh.z.front()+(fieldMesh.z.front()-zTmp);
			if(zTmp>fieldMesh.z.back()) zTmp = fieldMesh.z.back()-(zTmp-fieldMesh.z.back());
#endif

#if DEBUG_INTERP >= 1    
#if DEBUG_INTERP >= 2
			cout<<"ERROR: Should never get here with _PARTICLE_BOUNDARY ==2|3"<<endl;
#endif
            status = 1;
			return TYPE2(0);
	}
#endif
	//else
	//{
		_x = (xTmp-fieldMesh.r.front())/(fieldMesh.r.back()-fieldMesh.r.front())*(fieldMesh.r.size()-1);
		_z = (zTmp-fieldMesh.z.front())/(fieldMesh.z.back()-fieldMesh.z.front())*(fieldMesh.z.size()-1);
	//}

	//x0 = floor(_x);
	ix = (int) floor(_x);
    //x1 = ceil(_x);

	//z0 = floor(_z);
	iz = (int) floor(_z);
	//z1 = ceil(_z);

	// Catch for particle at point
	//if(x0==x1 || z0==z1) {
//#if DEBUG_INTERP >= 2
//        cout << "Particle version of kj_interp" << endl;
//		cout << "x0: " << x0 << " x1: " <<x1<< " _x: "<<_x << endl;
//		cout << "Particle at point catch: " << x0/x1 << "  "  << abs(1.0-x0/x1) << endl;
//#endif
//		return fVec[ix][iz];
//	}
	//else {

		TYPE2 f00 = fVec[ix][iz];
		TYPE2 f01 = fVec[ix][iz + 1];
		TYPE2 f10 = fVec[ix + 1][iz];
		TYPE2 f11 = fVec[ix + 1][iz + 1];
        
#if DEBUG_INTERP >=2
        //cout << "kj_interp: " << endl;
        //if(typeid(TYPE)==typeid(C3Vec)) cout << "Type is C3Vec" << endl;
        //if(typeid(TYPE)==typeid(float)) cout << "Type is float" << endl;
        //kj_print(x0,"x0");
        //kj_print(x1,"x1");
        //kj_print(y0,"y0");
        //kj_print(y1,"y1");
        //cout << endl;
        if(x0>fVec.size()-1||x0<0||x1>fVec.size()-1||x1<1) {
                cout<<"ERROR: interpolation point off the end of array"<<endl;
                cout<<"x.front: "<<fieldMesh.r.front()<<endl;
                cout<<"x.back: "<<fieldMesh.r.back()<<endl;
                cout<<"x: "<<x<<endl;
                cout<<"xVec.size(): "<<fieldMesh.r.size()<<endl;
                cout<<"fVec.size(): "<<fVec.size()<<endl;
                status = 1;
                return TYPE2(0);
        }
#endif
        x0 = fieldMesh.r[ix];
        x1 = fieldMesh.r[ix + 1];

        z0 = fieldMesh.z[iz];
        z1 = fieldMesh.z[iz + 1];
        
       // TYPE2 result = (1.0/( (x1 - x0)*(z1 - z0) ))*(f00*(x1 - _x)*(z1 - _z) + f01*(_x - x0)*(z1 - _z) + f10*(x1 - _x)*(_z - z0) + f11*(_x - x0)*(z - z0) );
        TYPE2 result = (1.0/( (x1 - x0)*(z1 - z0) ))*(f00*(x1 - x)*(z1 - z) + f01*(x - x0)*(z1 - z) + f10*(x1 - x)*(z - z0) + f11*(x - x0)*(z - z0) );

#if DEBUG_INTERP >=1
        if(isnan(result)) {
#if DEBUG_INTERP >= 2
                cout<<"ERROR: interpolation produced a NaN"<<endl;
#endif
                status = 1;
                return TYPE2(0);
        } 
        if(isinf(result)) {
#if DEBUG_INTERP >= 2
                cout<<"ERROR: interpolation produced a INF"<<endl;
#endif
                status = 1;
                return TYPE2(0);
        }
#endif
		return result;
	//}
}


C3Vec operator* ( const float A[][3], const C3Vec x ) {
        C3Vec B;
        B.c1 = A[0][0]*x.c1 + A[0][1]*x.c2 + A[0][2]*x.c3;
        B.c2 = A[1][0]*x.c1 + A[1][1]*x.c2 + A[1][2]*x.c3;
        B.c3 = A[2][0]*x.c1 + A[2][1]*x.c2 + A[2][2]*x.c3;
        return B;
}

void transpose ( float A[][3] ) {

        float B[3][3];
        
        B[0][0] = A[0][0];
        B[1][0] = A[0][1];
        B[2][0] = A[0][2];

        B[0][1] = A[1][0];
        B[1][1] = A[1][1];
        B[2][1] = A[1][2];

        B[0][2] = A[2][0];
        B[1][2] = A[2][1];
        B[2][2] = A[2][2];

        A = B;
}

C3Vec rot_CYL_to_XYZ ( const float t, const C3Vec vec, const int direction ) {

    // t here is the the cylindrical angle position in rtz (radians)        

    float rot[3][3];

    rot[0][0] =  cos(t);
    rot[0][1] = -sin(t);
    rot[0][2] = 0;

    rot[1][0] =  sin(t);
    rot[1][1] =  cos(t);
    rot[1][2] = 0;

    rot[2][0] = 0;
    rot[2][1] = 0;
    rot[2][2] = 1;

    if(direction<0) {
        transpose(rot);
    }
    
    return rot * vec;
}

C3Vec rot_XYZ_to_abp ( const C3Vec A_XYZ, C3Vec bUnit_XYZ, const int direction ) {

        // If direction<1 then the inverse rotation is applied, i.e., abp_to_XYZ

        C3Vec A_abp;

        C3Vec xu_xyz (1,0,0);
        C3Vec yu_xyz (0,1,0);
        C3Vec zu_xyz (0,0,1);

        /// normalize bUnit_XYZ just in case:
        C3Vec pu_xyz = bUnit_XYZ/mag(bUnit_XYZ);

        // alp is mostly in the +/- x / r direction depending on B toroidal direction
        // bet is mostly z direction

        C3Vec a_xyz = cross(zu_xyz,pu_xyz); 
        C3Vec au_xyz = a_xyz/mag(a_xyz);

        C3Vec b_xyz = cross(pu_xyz,au_xyz);
        C3Vec bu_xyz = b_xyz/mag(b_xyz);

#if DEBUG_ROTATION >=1 
        C3Vec au_xyz2 = au_xyz;
        C3Vec bu_xyz2 = bu_xyz;
        C3Vec pu_xyz2 = pu_xyz;

        cout<<"au_xyz: "<<au_xyz<<endl;
        cout<<"bu_xyz: "<<bu_xyz<<endl;
        cout<<"pu_xyz: "<<pu_xyz<<endl;
#endif

    // Rotation 1

    float theta = acos( dot(xu_xyz,au_xyz) );

#if DEBUG_ROTATION >=1 
        cout<<"theta: "<< theta*180.0/_pi<<endl;
#endif

    float q0  = cos ( theta / 2.0 );
    float q1  = sin ( theta / 2.0 ) * (-zu_xyz.c1); 
    float q2  = sin ( theta / 2.0 ) * (-zu_xyz.c2);  
    float q3  = sin ( theta / 2.0 ) * (-zu_xyz.c3); 

    // Construct the rotation matrix

    float rot1[3][3];

    rot1[0][0] = pow(q0,2)+pow(q1,2)-pow(q2,2)-pow(q3,2);
    rot1[0][1] = 2*(q1*q2-q0*q3);
    rot1[0][2] = 2*(q1*q3+q0*q2);
    rot1[1][0] = 2*(q2*q1+q0*q3);
    rot1[1][1] = pow(q0,2)-pow(q1,2)+pow(q2,2)-pow(q3,2);
    rot1[1][2] = 2*(q2*q3-q0*q1);
    rot1[2][0] = 2*(q3*q1-q0*q2);
    rot1[2][1] = 2*(q3*q2+q0*q1);
    rot1[2][2] = pow(q0,2)-pow(q1,2)-pow(q2,2)+pow(q3,2);

    if(direction<0) {
        transpose(rot1);
    }

    au_xyz = rot1 * au_xyz;
    bu_xyz = rot1 * bu_xyz;
    pu_xyz = rot1 * pu_xyz;

#if DEBUG_ROTATION >=1 
    cout<<"au_rtz 1: "<<au_xyz<<endl;
    cout<<"bu_rtz 1: "<<bu_xyz<<endl;
    cout<<"pu_rtz 1: "<<pu_xyz<<endl;
#endif

    // Rotation 2

    theta = acos( dot(zu_xyz,pu_xyz) );

#if DEBUG_ROTATION >=1
    cout<<"theta: "<<theta * 180.0/_pi << endl;
#endif
    q0  = cos ( theta / 2.0 );
    q1  = sin ( theta / 2.0 ) * (-xu_xyz.c1); 
    q2  = sin ( theta / 2.0 ) * (-xu_xyz.c2);  
    q3  = sin ( theta / 2.0 ) * (-xu_xyz.c3); 

    // Construct the rotation matrix

    float rot2[3][3];

    rot2[0][0] = pow(q0,2)+pow(q1,2)-pow(q2,2)-pow(q3,2);
    rot2[0][1] = 2*(q1*q2-q0*q3);
    rot2[0][2] = 2*(q1*q3+q0*q2);
    rot2[1][0] = 2*(q2*q1+q0*q3);
    rot2[1][1] = pow(q0,2)-pow(q1,2)+pow(q2,2)-pow(q3,2);
    rot2[1][2] = 2*(q2*q3-q0*q1);
    rot2[2][0] = 2*(q3*q1-q0*q2);
    rot2[2][1] = 2*(q3*q2+q0*q1);
    rot2[2][2] = pow(q0,2)-pow(q1,2)-pow(q2,2)+pow(q3,2);

    if(direction<0) {
        transpose(rot2);
    }

    au_xyz = rot2 * au_xyz;
    bu_xyz = rot2 * bu_xyz;
    pu_xyz = rot2 * pu_xyz;

#if DEBUG_ROTATION >=1 
    cout<<"au_xyz 2: "<<au_xyz<<endl;
    cout<<"bu_xyz 2: "<<bu_xyz<<endl;
    cout<<"pu_xyz 2: "<<pu_xyz<<endl;
#endif

    A_abp = rot2 * ( rot1 * A_XYZ );
    
    return A_abp;
}

float GetGyroPhase ( const C3Vec v_abp ) {

        // alp is mostly in the x / r direction
        // bet is mostly z direction

        float alp = v_abp.c1;
        float bet = v_abp.c2;

        return atan2(alp,bet);
}

float GetAlpComp ( const float vPer, const float phs ) {

        return vPer * sin(phs);
}

float GetBetComp ( const float vPer, const float phs ) {

        return vPer * cos(phs);
}


template<class b0TYPE>
C3Vec GetFb0( CParticle &p,const C3Vec &v_XYZ,
				const C3Vec &x, const b0TYPE &b0Vec_CYL,
			  	const fieldMeshClass &fieldMesh, int &status ) {

	float _r = sqrt ( pow(x.c1,2) + pow(x.c2,2) );
	float _p = atan2 ( x.c2, x.c1 );
    
    ///_p = 0; // Really need to do something about this (DG). ///// WHY??? (DB)

#if DEBUGLEVEL >= 3
	cout << "\t\t\tx: " << x.c1 << " y: " << x.c2 << " z: " << x.c3 << endl;
	cout << "\t\t\tr: " << _r << " p: " << _p << endl;
	cout << "\t\t\trVec.front(): " << rVec.front() << endl;
	cout << "\t\t\tv_XYZ: " << v_XYZ.c1 << "  " << v_XYZ.c2 << "  " << v_XYZ.c3 << endl;
#endif

	C3Vec b0_CYL, b0_XYZ;

	b0_CYL = kj_interp(C3Vec(_r,_p,x.c3) , fieldMesh, b0Vec_CYL, p, status);

	b0_XYZ = C3Vec(cos(_p)*b0_CYL.c1-sin(_p)*b0_CYL.c2+0,
					sin(_p)*b0_CYL.c1+cos(_p)*b0_CYL.c2+0,
					0+0+1*b0_CYL.c3);

    C3Vec v_x_b0 = cross(v_XYZ,b0_XYZ);

#if DEBUGLEVEL >= 3
	cout << "\tvxb0: " << v_x_b0.c1 << "  " << v_x_b0.c2 << "  " << v_x_b0.c3 << endl;
	cout << "\tp.q/p.m: " << p.q/p.m << endl;
#endif
	return v_x_b0*(p.q/p.m);
}

// Zero-order orbits
template<class b0TYPE>
C3Vec rk4_evalf ( CParticle &p, const float &t, 
				const C3Vec &v_XYZ, const C3Vec &x, const b0TYPE &b0Vec_CYL,
			  	const fieldMeshClass &fieldMesh, int &status ) {
    
        C3Vec Fb0 = GetFb0(p, v_XYZ, x, b0Vec_CYL, fieldMesh, status);

        return Fb0;
        }

// Zero-order orbits
template<class b0TYPE>
int rk4_move (CParticle &p, const float &dt, const float &t0,
				const b0TYPE &b0, const fieldMeshClass &fieldMesh) {

        int status = 0;
		C3Vec yn0(p.v_c1,p.v_c2,p.v_c3), xn0(p.c1, p.c2, p.c3);
		C3Vec k1, k2, k3, k4, yn1, x1, x2, x3, x4, xn1; 

		k1 = rk4_evalf ( p, t0 + 0.0*dt, yn0         , xn0         , b0, fieldMesh, status ) * dt;
		x1 = yn0 * dt;

		k2 = rk4_evalf ( p, t0 + 0.5*dt, yn0 + 0.5*k1, xn0 + 0.5*x1, b0, fieldMesh, status ) * dt;
		x2 = (yn0 + 0.5*k1) * dt;

		k3 = rk4_evalf ( p, t0 + 0.5*dt, yn0 + 0.5*k2, xn0 + 0.5*x2, b0, fieldMesh, status ) * dt;
		x3 = (yn0 + 0.5*k2) * dt;

		k4 = rk4_evalf ( p, t0 + 1.0*dt, yn0 + 1.0*k3, xn0 + 1.0*x3, b0, fieldMesh, status ) * dt;
		x4 = (yn0 + 1.0*k3) * dt;

		yn1 = yn0 + 1.0/6.0 * (k1+2.0*k2+2.0*k3+k4);
		xn1 = xn0 + 1.0/6.0 * (x1+2.0*x2+2.0*x3+x4);

		p.c1 = xn1.c1;
		p.c2 = xn1.c2;
		p.c3 = xn1.c3;
		p.v_c1 = yn1.c1;
		p.v_c2 = yn1.c2;
		p.v_c3 = yn1.c3;


#if _PARTICLE_BOUNDARY == 1
		// Particle absorbing walls
#elif _PARTICLE_BOUNDARY == 2
		// Periodic 
		if(p.c1<r.front())
		{	
#if DEBUGLEVEL >= 1
			cout<<"Particle went left"<<endl;
#endif
			p.c1 = fieldMesh.r.back()-(fieldMesh.r.front()-p.c1);
		}
		if(p.c1>fieldMesh.r.back())
		{
#if DEBUGLEVEL >= 1
			cout<<"Particle went right"<<endl;
#endif
			p.c1 = fieldMesh.r.front()+(p.c1-fieldMesh.r.back());
		}
#elif _PARTICLE_BOUNDARY == 3
		// Particle reflecting walls
		if(p.c1<fieldMesh.r.front())
		{
			cout<<"Particle hit the left wall"<<endl;
			cout<<"r.front(): "<<fieldMesh.r.front()<<endl;
			p.c1 = fieldMesh.r.front()+(fieldMesh.r.front()-p.c1);
			p.v_c1 = -p.v_c1;
		}
		if(p.c1>r.back())
		{
			cout<<"Particle hit the right wall"<<endl;
			cout<<"r.back(): "<<fieldMesh.r.back()<<endl;
			p.c1 = fieldMesh.r.back()-(p.c1-r.back());
			p.v_c1 = -p.v_c1;
		}
#endif

#if DEBUGLEVEL >= 3
        cout << "\tdt: " << dt << endl;
        cout << "\tr.front(): " << r.front() << endl;
        cout << "\tr.back(): " << r.back() << endl;
		cout << "\tx0_XYZ: " << xn0.c1 << "  " << xn0.c2 << "  " << xn0.c3 << endl;
		cout << "\tv0_XYZ: " << yn0.c1 << "  " << yn0.c2 << "  " << yn0.c3 << endl;
		cout << "\tx1_XYZ: " << xn1.c1 << "  " << xn1.c2 << "  " << xn1.c3 << endl;
		cout << "\tv1_XYZ: " << yn1.c1 << "  " << yn1.c2 << "  " << yn1.c3 << endl;
		cout << "\tE: " << 0.5 * p.m * sqrt (pow(p.v_c1,2)+pow(p.v_c2,2)+pow(p.v_c3,2))/_e << endl;
#endif
        return status;
}

// First-order orbits
template<class b0TYPE, class e1TYPE>
C3Vec rk4_evalf ( CParticle &p, const float &t, const C3Vec &v_XYZ, const C3Vec &x,
				const b0TYPE &b0Vec_CYL,const fieldMeshClass &fieldMesh,
                 const e1TYPE &e1REVec_XYZ,const e1TYPE &e1IMVec_XYZ,
                  const float wrf, int &status ) {

    C3Vec Fb0 = GetFb0(p, v_XYZ, x,b0Vec_CYL,fieldMesh, status);

	float _r = sqrt ( pow(x.c1,2) + pow(x.c2,2) );
	float _p = atan2 ( x.c2, x.c1 );
    _p = 0; // Really need to do something about this. (DG) /// WHY? (DB)

    C3Vec Fe1, e1RE_XYZ, e1IM_XYZ, e1_XYZ;
    
    float ex_a, ex_p, ey_a, ey_p, ez_a, ez_p;
    
    e1RE_XYZ = kj_interp (C3Vec(_r,_p,x.c3), fieldMesh, e1REVec_XYZ, p, status );
	e1IM_XYZ = kj_interp (C3Vec(_r,_p,x.c3), fieldMesh, e1IMVec_XYZ, p, status );
    
    ex_a = sqrt( pow( e1RE_XYZ.c1,2) + pow( e1IM_XYZ.c1,2));
    ey_a = sqrt( pow( e1RE_XYZ.c2,2) + pow( e1IM_XYZ.c2,2));
    ez_a = sqrt( pow( e1RE_XYZ.c3,2) + pow( e1IM_XYZ.c3,2));
    
    ex_p = atan2(e1IM_XYZ.c1,e1RE_XYZ.c1);
    ey_p = atan2(e1IM_XYZ.c2,e1RE_XYZ.c2);
    ez_p = atan2(e1IM_XYZ.c3,e1RE_XYZ.c3);

	e1_XYZ = C3Vec(ex_a*cos(ex_p + wrf*t + p.rfPhase) , ey_a*cos(ey_p + wrf*t + p.rfPhase), ez_a*cos(ez_p + wrf*t + p.rfPhase) );

    Fe1 = (p.q/p.m)*e1_XYZ;
    return Fb0 + (1.0e3)*Fe1;
    
}

// First-order orbits
template<class b0TYPE, class e1TYPE>
int rk4_move ( CParticle &p, float dt, float t0,
				const b0TYPE&b0, const fieldMeshClass &fieldMesh,
                 const e1TYPE &e1Re,  const e1TYPE &e1Im, const float wrf ) {

        int status = 0;
		C3Vec yn0(p.v_c1,p.v_c2,p.v_c3), xn0(p.c1, p.c2, p.c3);
		C3Vec k1, k2, k3, k4, yn1, x1, x2, x3, x4, xn1; 

		k1 = rk4_evalf ( p, t0 + 0.0*dt, yn0 + 0.*yn0, xn0         , b0, fieldMesh, e1Re, e1Im, wrf, status ) * dt;
		x1 = yn0 * dt;

		k2 = rk4_evalf ( p, t0 + 0.5*dt, yn0 + 0.5*k1, xn0 + 0.5*x1, b0, fieldMesh, e1Re, e1Im, wrf, status ) * dt;
		x2 = (yn0 + 0.5*k1) * dt;

		k3 = rk4_evalf ( p, t0 + 0.5*dt, yn0 + 0.5*k2, xn0 + 0.5*x2, b0, fieldMesh, e1Re, e1Im, wrf, status ) * dt;
		x3 = (yn0 + 0.5*k2) * dt;

		k4 = rk4_evalf ( p, t0 + 1.0*dt, yn0 + 1.0*k3, xn0 + 1.0*x3, b0, fieldMesh, e1Re, e1Im, wrf, status ) * dt;
		x4 = (yn0 + 1.0*k3) * dt;

		yn1 = yn0 + 1.0/6.0 * (k1+2.0*k2+2.0*k3+k4);
		xn1 = xn0 + 1.0/6.0 * (x1+2.0*x2+2.0*x3+x4);
    
        p.c1 = xn1.c1;
        p.c2 = xn1.c2;
		p.c3 = xn1.c3;
		p.v_c1 = yn1.c1;
		p.v_c2 = yn1.c2;
		p.v_c3 = yn1.c3;

        return status;
}

///////////////////////////////////// NEED TO CHANGE FIRST ARGUMENT OF kj_interp TO A C3Vec FOR GUIDING CENTER
//////////////////////////////////////////////

// Parallel acceleration
float eval_aPar ( CParticle &p, const C3Vec r, const vector<float> &r_GC, const vector<float> &bDotGradB, int &status ) {

    float This_bDotGradB = kj_interp ( r.c1, r_GC, bDotGradB, p, status );
#if DEBUG_EVAL_APAR >= 1
    if(status>0) {
            cout<<"ERROR 1 in eval_aPar"<<endl;
            exit(1);
    }
#endif
    float aPar = -p.u / p.m * This_bDotGradB;
#if DEBUG_EVAL_APRA >= 1
    if(isnan(aPar)||isinf(aPar)) {
            status = 1;
            cout<<"ERROR 2 in eval_aPar"<<endl;
            exit(1);
    }
#endif
    return aPar;
}

// Perpendicular velocity
float eval_vPer ( CParticle &p, const C3Vec r, const vector<float> &r_b0, const vector<C3Vec> &b0_CYL, int &status ) {

	C3Vec This_b0_CYL = kj_interp ( r.c1, r_b0, b0_CYL, p, status );
    return sqrt ( 2.0 * p.u * mag(This_b0_CYL) / p.m );
}

// Guiding center veclocity
C3Vec eval_vGC ( CParticle &p, const C3Vec r, const float vPer, const float vPar, 
                const vector<float> &r_b0, const vector<C3Vec> &b0_CYL, 
                const vector<float> &r_GC, const vector<C3Vec> &curv_CYL, const vector<C3Vec> &grad_CYL, int &status ) {

	C3Vec This_b0_CYL = kj_interp ( r.c1, r_b0, b0_CYL, p, status );
#if DEBUG_EVAL_VGC >= 1
    if(status>0) {
            cout<<"ERROR 1 in eval_vGC"<<endl;
            exit(1);
    }
#endif

	C3Vec This_curv_CYL = kj_interp ( r.c1, r_GC, curv_CYL, p, status );
#if DEBUG_EVAL_VGC >= 1
    if(status>0) {
            cout<<"ERROR 2 in eval_vGC"<<endl;
            exit(1);
    }
#endif

	C3Vec This_grad_CYL = kj_interp ( r.c1, r_GC, grad_CYL, p, status );
#if DEBUG_EVAL_VGC >= 1
    if(status>0) {
            cout<<"ERROR 3 in eval_vGC"<<endl;
            exit(1);
    }
#endif

#if DEBUG_EVAL_VGC >= 1

    cout << "r.c1: " << r.c1 << endl;
    cout << "p.c1: " << p.c1 << endl;
    cout << "vPar: " << vPar << endl;
    cout << "vPer: " << vPer << endl;
    cout << "b0_CYL: " << This_b0_CYL.c1 << "  " << This_b0_CYL.c2 << "  " << This_b0_CYL.c3 << endl;
    cout << "curv_CYL: " << This_curv_CYL.c1 << "  " << This_curv_CYL.c2 << "  " << This_curv_CYL.c3 << endl;
    cout << "grad_CYL: " << This_grad_CYL.c1 << "  " << This_grad_CYL.c2 << "  " << This_grad_CYL.c3 << endl << endl;
    cout << "max(grad_CYL): " << maxC3VecAbs(grad_CYL) << endl;

#endif

    C3Vec UnitB_CYL = This_b0_CYL / mag(This_b0_CYL);

    C3Vec vGC = vPar * UnitB_CYL  +  pow(vPer,2) * This_grad_CYL  +  pow(vPar,2) * This_curv_CYL;
    return vGC; 
}

// Guiding center orbit
int rk4_move_gc ( CParticle &p, const float &dt, const float &t0, 
				const vector<float> &r_b0, const vector<C3Vec> &b0_CYL, const vector<float> &r_GC, 
                const vector<C3Vec> &curv_CYL, const vector<C3Vec> &grad_CYL, 
                const vector<float> &bDotGradB, const float wrf ) {

                int status=0;
                C3Vec xn0_XYZ(p.c1, p.c2, p.c3);
                C3Vec xn0 = XYZ_to_CYL(xn0_XYZ);

		    	float This_vPer = eval_vPer ( p, xn0, r_b0, b0_CYL, status );
#if DEBUG_GC >= 2 
                cout << "p.vPer: " << p.vPer << endl;
                cout << "p.vPar: " << p.vPar << endl;
                cout << "This_vPer: " << This_vPer << endl;
#endif
		    	C3Vec This_vGC  = eval_vGC  ( p, xn0, This_vPer, p.vPar + 0, r_b0, b0_CYL, r_GC, curv_CYL, grad_CYL, status );
		    	float k1_vPar = dt * eval_aPar ( p, xn0, r_GC, bDotGradB, status );
		    	C3Vec k1_vgc  = dt * This_vGC;
#if DEBUG_GC >= 2
                kj_print(k1_vgc,"k1_vgc");
                kj_print(xn0,"xn0");
                cout<<"Status: "<<status<<endl;
                if(isnan(k1_vgc)||isinf(k1_vgc)||isnan(xn0)||isinf(xn0)||status>0) {
                        status = 1;
                        return status;
                }
#endif    
		    	This_vPer = eval_vPer ( p, xn0 + k1_vgc / 2.0, r_b0, b0_CYL, status );
		    	This_vGC  = eval_vGC  ( p, xn0 + k1_vgc / 2.0, This_vPer, p.vPar + k1_vPar / 2.0, r_b0, b0_CYL, r_GC, curv_CYL, grad_CYL, status );
		    	float k2_vPar = dt * eval_aPar ( p, xn0 + k1_vgc / 2.0, r_GC, bDotGradB, status ); 
		    	C3Vec k2_vgc  = dt * This_vGC;
#if DEBUG_GC >= 2
                kj_print(k2_vgc,"k2_vgc");
                if(isnan(k2_vgc)||isinf(k2_vgc)||isnan(xn0)||isinf(xn0)||status>0) {
                        status = 1;
                        return status;
                }
#endif 
		    	This_vPer = eval_vPer ( p, xn0 + k2_vgc / 2.0, r_b0, b0_CYL, status ); 
		    	This_vGC  = eval_vGC  ( p, xn0 + k2_vgc / 2.0, This_vPer, p.vPar + k2_vPar / 2.0, r_b0, b0_CYL, r_GC, curv_CYL, grad_CYL, status );
		    	float k3_vPar = dt * eval_aPar ( p, xn0 + k2_vgc / 2.0, r_GC, bDotGradB, status ); 
		    	C3Vec k3_vgc  = dt * This_vGC;
#if DEBUG_GC >= 2
                kj_print(k3_vgc,"k3_vgc");
                if(isnan(k3_vgc)||isinf(k3_vgc)||isnan(xn0)||isinf(xn0)||status>0) {
                        status = 1;
                        return status;
                }
#endif 
		    	This_vPer = eval_vPer ( p, xn0 + k3_vgc, r_b0, b0_CYL, status ); 
		    	This_vGC  = eval_vGC  ( p, xn0 + k3_vgc, This_vPer, p.vPar + k3_vPar, r_b0, b0_CYL, r_GC, curv_CYL, grad_CYL, status );
		    	float k4_vPar = dt * eval_aPar ( p, xn0 + k3_vgc, r_GC, bDotGradB, status );
		    	C3Vec k4_vgc  = dt * This_vGC; 
#if DEBUG_GC >= 2
                kj_print(k4_vgc,"k4_vgc");
                if(isnan(k4_vgc)||isinf(k4_vgc)||isnan(xn0)||isinf(xn0)||status>0) {
                        status = 1;
                        return status;
                }
#endif 	
		    	float vPar1 = p.vPar + ( k1_vPar + 2.0 * k2_vPar + 2.0 * k3_vPar + k4_vPar ) / 6.0;
		    	C3Vec xn1 = xn0 + ( k1_vgc + 2.0 * k2_vgc + 2.0 * k3_vgc + k4_vgc ) / 6.0;

#if DEBUG_GC >=1 
                if(isnan(xn1)||isinf(xn1)) {
                        status = 1;
                        return status;
                }
#endif

                // Update particle with moved position and new vPar & vPer

		    	float vPer1 = eval_vPer ( p, xn1, r_b0, b0_CYL, status ); 

                p.vPar = vPar1;
                p.vPer = vPer1;

                C3Vec xn1_XYZ = CYL_to_XYZ(xn1);

                p.c1 = xn1_XYZ.c1; 
		        p.c2 = xn1_XYZ.c2;
		        p.c3 = xn1_XYZ.c3;

                // Update the XYZ velocity also

                C3Vec this_b0_CYL = kj_interp ( xn1.c1, r_b0, b0_CYL, status );
                C3Vec this_b0_XYZ = rot_CYL_to_XYZ ( xn1.c2, this_b0_CYL, 1 );
 
                C3Vec v_abp;

                float this_wc = p.q * mag(this_b0_CYL) / p.m;
                p.phs = this_wc*t0+p.gyroPhase;
                v_abp.c1 = GetAlpComp(vPer1,p.phs); 
                v_abp.c2 = GetBetComp(vPer1,p.phs);
                v_abp.c3 = vPar1;

                p.vAlp = v_abp.c1;
                p.vBet = v_abp.c2;

                C3Vec this_v_XYZ = rot_XYZ_to_abp ( v_abp, this_b0_XYZ, -1 ); 

                p.v_c1 = this_v_XYZ.c1;
                p.v_c2 = this_v_XYZ.c2;
                p.v_c3 = this_v_XYZ.c3;

                return status;
}

float maxwellian ( float vx, float vy, float vz, float vTh ) {

    float weight_x = 1.0 / (vTh*sqrt(_pi)) * exp ( -pow(vx,2) / pow(vTh,2) );
    float weight_y = 1.0 / (vTh*sqrt(_pi)) * exp ( -pow(vy,2) / pow(vTh,2) );
    float weight_z = 1.0 / (vTh*sqrt(_pi)) * exp ( -pow(vz,2) / pow(vTh,2) );
    
    return weight_x * weight_y * weight_z;

}

float get_vTh ( const float _amu, const float _Z, const float _T_keV ) {

    float m = _amu * _mi;
    float q = _Z * _e;
    float kT_joule = _T_keV * 1e3 * _e; // This may actually be E_keV so may need a 3/2 somewhere
    float vTh = sqrt ( 2.0*kT_joule / m );

	return vTh;

}

C3Vec maxwellian_df0_dv (const C3Vec _v, const float _T_keV, const float _n_m3, const float _amu, const float _Z ) {

	C3Vec df0_dv;

	float vTh = get_vTh ( _amu, _Z, _T_keV );

	float _vx = _v.c1;
	float _vy = _v.c2;
	float _vz = _v.c3;

    // Get the 3 components of df0_dv at this point in velocity space

    float h = vTh/1000.0;
    float vxL = _vx-h;
    float vxR = _vx+h;
    float fL = maxwellian(vxL,_vy,_vz,vTh);
    float fR = maxwellian(vxR,_vy,_vz,vTh);
    float _df0_dv = (-fL + fR)/(2*h); 

    df0_dv.c1 = _df0_dv * _n_m3;

    float vyL = _vy-h;
    float vyR = _vy+h;
    fL = maxwellian(_vx,vyL,_vz,vTh);
    fR = maxwellian(_vx,vyR,_vz,vTh);
    _df0_dv = (-fL + fR)/(2*h); 

    df0_dv.c2 = _df0_dv * _n_m3;

    float vzL = _vz-h;
    float vzR = _vz+h;
    fL = maxwellian(_vx,_vy,vzL,vTh);
    fR = maxwellian(_vx,_vy,vzR,vTh);
    _df0_dv = (-fL + fR)/(2*h); 

    df0_dv.c3 = _df0_dv * _n_m3;

	return df0_dv;
}

vector<CParticle> create_particles ( float x, float amu, float Z, float T_keV, float n_m3, 
                int nPx, int nPy, int nPz, int nThermal, float &dv, C3Vec b0_XYZ, float rfPhase) {

        vector<CParticle> pList;

        int nP = nPx * nPy * nPz;
        pList.resize(nP);

        float m = amu * _mi;
		float vTh = get_vTh ( amu, Z, T_keV );

# if DEBUG_MAXWELLIAN >= 1
        cout <<"amu: "<< amu<<endl;
        cout <<"Z: "<< Z<<endl;
        //cout <<"m: "<< m<<endl;
        //cout <<"q: "<< q<<endl;
        cout <<"vTh: "<< vTh<<endl;
#endif

        float vxMin, vyMin, vzMin, dvx, dvy, dvz;

        if(nPx == 1){
        
            vxMin=0.0;
            dvx=1.0;
        }
    
        else{
        
            float vxRange = vTh * nThermal * 2;
             vxMin	=  -vxRange / 2.0;
             dvx = vxRange / (nPx-1);
        }

       if(nPy == 1){
            vyMin=0.0;
            dvy=1.0;
        }
    
        else{
     
            float vyRange = vTh * nThermal * 2;
             vyMin	= -vyRange / 2.0;
             dvy = vyRange / (nPy-1);
        }

       if(nPz == 1){
            vzMin=0.0;
            dvz=1.0;
        }
    
        else{
            float vzRange = vTh * nThermal * 2;
             vzMin	= -vzRange / 2.0;
             dvz = vzRange / (nPz-1);
        }
    
        dv = dvx*dvy*dvz; // Return the Jacobian (volume element for integration later)


        float TestIntegratedValue = 0;

        int cnt = 0;
        for(int i=0;i<nPx;i++) {
            for(int j=0;j<nPy;j++) {
                for(int k=0;k<nPz;k++) {

                    float thisvx = vxMin+i*dvx;
                    float thisvy = vyMin+j*dvy;
                    float thisvz = vzMin+k*dvz;

                    float weight = maxwellian(thisvx,thisvy,thisvz,vTh) * n_m3;

                    TestIntegratedValue += weight * dv;

    	            CParticle p (x,0,0,thisvx,thisvy,thisvz,amu,Z,weight,rfPhase);
                    pList[cnt] = p;
                    pList[cnt].number = cnt;
                    pList[cnt].vTh = vTh;

                    pList[cnt].d3v = dv;

                    //// Get the 3 components of df0_dv at this point in velocity space

					//C3Vec thisVel(thisvx,thisvy,thisvz);
					//C3Vec _df0_dv = maxwellian_df0_dv ( thisVel, T_keV, amu, Z );

                    // Get vPar, vPer and mu for guiding center integration

                    C3Vec thisV_XYZ(thisvx,thisvy,thisvz); 
                    float bMag = mag (b0_XYZ);
                    float vMag = mag (thisV_XYZ);

                    //float vPar = (thisvx*b0_XYZ.c1 + thisvy*b0_XYZ.c2 + thisvz*b0_XYZ.c3) / bMag;
		            //float vPer = sqrt ( pow(vMag,2) - pow(vPar,2) );

                    C3Vec thisV_abp = rot_XYZ_to_abp ( thisV_XYZ, b0_XYZ, 0 );

                    pList[cnt].vPar = thisV_abp.c3;
                    
                    pList[cnt].vPer = sqrt(pow(thisV_abp.c1,2)+pow(thisV_abp.c2,2));
                    pList[cnt].gyroPhase = GetGyroPhase(thisV_abp); 
                    pList[cnt].u = pList[cnt].m * pow(pList[cnt].vPer,2) / ( 2.0 * bMag );

#if DEBUG_MAXWELLIAN >=2 
                    cout<<"ThisVx: "<<thisvx<<endl;
                    cout<<"ThisVy: "<<thisvy<<endl;
                    cout<<"ThisVz: "<<thisvz<<endl;
                    cout<<"vMag: "<<vMag<<endl;
                    cout<<"vPer: "<<pList[cnt].vPer<<endl;
                    cout<<"vPar: "<<pList[cnt].vPar<<endl;
                    cout<<"u: "<<pList[cnt].u<<endl<<endl;
                    if(isnan(pList[cnt].u)) exit(1);
#endif                    
                    cnt++;
                }
            }
        } 

#if DEBUG_MAXWELLIAN >=1 
        cout << "TestIntegratedValue: " << TestIntegratedValue << endl;
#endif
        return pList;
}


// Blob of particles to integrate at each point in velocity space
// Initially created as the phase space particles of Larmor radius about the configuration space point
vector<CParticle> create_particle_blob ( CParticle P, float amu, float Z, float T_keV,
                                        float n_m3, int nPblob, int nThermal, float &dv, C3Vec b0_XYZ) {

        vector<CParticle> pList;

        pList.resize(nPblob);

        float m = amu * _mi;
		float vTh = get_vTh ( amu, Z, T_keV );
    
        float TestIntegratedValue = 0;
    
        C3Vec thisV_XYZ(P.v_c1,P.v_c2,P.v_c3);

        C3Vec b0_perp1_XYZ = cross(b0_XYZ, C3Vec(0.0, 0.0, 1.0));
        b0_perp1_XYZ = b0_perp1_XYZ/mag(b0_perp1_XYZ);
    
    
        C3Vec b0_perp2_XYZ = cross(b0_XYZ, b0_perp1_XYZ);
        b0_perp2_XYZ = b0_perp2_XYZ/mag(b0_perp2_XYZ);
    
    
        float rL = (m*P.v_c2)/(abs(P.q)*mag(b0_XYZ));
        float magb0 = mag(b0_XYZ);
        C3Vec P_CYL = C3Vec(P.c1, P.c2, P.c2);
        C3Vec P_XYZ = CYL_to_XYZ(P_CYL);
    
        int cnt = 0;
//i=0 case:
        float weight = maxwellian(P.v_c1,P.v_c2,P.v_c3,vTh) * n_m3; // n_m3 = density? (why "_3" ?)

        cnt = 0;
        for(int i=0;i<nPblob;i++) {

 	            CParticle p_tmp (P_XYZ.c1 + rL*(cos(2.0*_pi*i/nPblob)*b0_perp1_XYZ.c1 + sin(2.0*_pi*i/nPblob)*b0_perp2_XYZ.c1),
                                P_XYZ.c2 + rL*(cos(2.0*_pi*i/nPblob)*b0_perp1_XYZ.c2 + sin(2.0*_pi*i/nPblob)*b0_perp2_XYZ.c2),
                                P_XYZ.c3 + rL*(cos(2.0*_pi*i/nPblob)*b0_perp1_XYZ.c3 + sin(2.0*_pi*i/nPblob)*b0_perp2_XYZ.c3),
                                -P.v_c2*(-sin(2.0*_pi*i/nPblob)*b0_perp1_XYZ.c1 + cos(2.0*_pi*i/nPblob)*b0_perp2_XYZ.c1) + P.v_c1*b0_XYZ.c1/magb0,
                                -P.v_c2*(-sin(2.0*_pi*i/nPblob)*b0_perp1_XYZ.c2 + cos(2.0*_pi*i/nPblob)*b0_perp2_XYZ.c2) + P.v_c1*b0_XYZ.c2/magb0,
                                -P.v_c2*(-sin(2.0*_pi*i/nPblob)*b0_perp1_XYZ.c3 + cos(2.0*_pi*i/nPblob)*b0_perp2_XYZ.c3) + P.v_c1*b0_XYZ.c3/magb0,
                                amu,Z,weight,P.rfPhase);
            
                pList[cnt] = p_tmp;
                pList[cnt].number = cnt;
                pList[cnt].vTh = vTh;

                pList[cnt].d3v = dv;

                // Get vPar, vPer and mu for guiding center integration

                float bMag = mag (b0_XYZ);
                float vMag = mag (thisV_XYZ);

#if DEBUG_MAXWELLIAN >=2 
                    cout<<"ThisVx: "<<thisvx<<endl;
                    cout<<"ThisVy: "<<thisvy<<endl;
                    cout<<"ThisVz: "<<thisvz<<endl;
                    cout<<"vMag: "<<vMag<<endl;
                    cout<<"vPer: "<<pList[cnt].vPer<<endl;
                    cout<<"vPar: "<<pList[cnt].vPar<<endl;
                    cout<<"u: "<<pList[cnt].u<<endl<<endl;
                    if(isnan(pList[cnt].u)) exit(1);
#endif                    
                cnt++;
        }

#if DEBUG_MAXWELLIAN >=1 
        cout << "TestIntegratedValue: " << TestIntegratedValue << endl;
#endif
        return pList;
}


/// overlaoded create_particle_blob version to average over both gyrophase and toroidal average
/// blob is still a 1D vector to make consistent with rest of code
vector<CParticle> create_particle_blob ( CParticle P, float amu, float Z, float T_keV,
                                        float n_m3, int nPGyro, int nPTor, int nThermal, float &dv, C3Vec b0_XYZ) {

        vector<CParticle> pList;
        int nPblob = nPGyro*nPTor;
        pList.resize(nPblob);

        float m = amu * _mi;
		float vTh = get_vTh ( amu, Z, T_keV );
    
        float TestIntegratedValue = 0;
    
        C3Vec thisV_XYZ(P.v_c1,P.v_c2,P.v_c3);

        C3Vec b0_perp1_XYZ = cross(b0_XYZ, C3Vec(0.0, 0.0, 1.0));
        b0_perp1_XYZ = b0_perp1_XYZ/mag(b0_perp1_XYZ);
    
    
        C3Vec b0_perp2_XYZ = cross(b0_XYZ, b0_perp1_XYZ);
        b0_perp2_XYZ = b0_perp2_XYZ/mag(b0_perp2_XYZ);
    
    
        float rL = (m*P.v_c2)/(abs(P.q)*mag(b0_XYZ));
        float magb0 = mag(b0_XYZ);
        C3Vec P_CYL = C3Vec(P.c1, P.c2, P.c2);
        C3Vec P_XYZ = CYL_to_XYZ(P_CYL);
    
        int cnt = 0;
//i=0 case:
        float weight = maxwellian(P.v_c1,P.v_c2,P.v_c3,vTh) * n_m3; // n_m3 = density? (why "_3" ?)

        cnt = 0;
        for(int i=0;i<nPblob;i++) {
               int iGyro = i % nPGyro;
                int iTor = int(floor(i/nPGyro)) % nPTor;

 	            CParticle p_tmp (P_XYZ.c1*cos(2.0*_pi*iTor/nPTor) - P_XYZ.c2*sin(2.0*_pi*iTor/nPTor) + rL*(cos(2.0*_pi*iGyro/nPGyro)*b0_perp1_XYZ.c1 + sin(2.0*_pi*iGyro/nPGyro)*b0_perp2_XYZ.c1),
                                 P_XYZ.c1*sin(2.0*_pi*iTor/nPTor) + P_XYZ.c2*cos(2.0*_pi*iTor/nPTor) + + rL*(cos(2.0*_pi*iGyro/nPGyro)*b0_perp1_XYZ.c2 + sin(2.0*_pi*iGyro/nPGyro)*b0_perp2_XYZ.c2),
                                P_XYZ.c3 + rL*(cos(2.0*_pi*iGyro/nPGyro)*b0_perp1_XYZ.c3 + sin(2.0*_pi*iGyro/nPGyro)*b0_perp2_XYZ.c3),
                                -P.v_c2*(-sin(2.0*_pi*iGyro/nPGyro)*b0_perp1_XYZ.c1 + cos(2.0*_pi*iGyro/nPGyro)*b0_perp2_XYZ.c1) + P.v_c1*b0_XYZ.c1/magb0,
                                -P.v_c2*(-sin(2.0*_pi*iGyro/nPGyro)*b0_perp1_XYZ.c2 + cos(2.0*_pi*iGyro/nPGyro)*b0_perp2_XYZ.c2) + P.v_c1*b0_XYZ.c2/magb0,
                                -P.v_c2*(-sin(2.0*_pi*iGyro/nPGyro)*b0_perp1_XYZ.c3 + cos(2.0*_pi*iGyro/nPGyro)*b0_perp2_XYZ.c3) + P.v_c1*b0_XYZ.c3/magb0,
                                amu,Z,weight,P.rfPhase);
            
                pList[cnt] = p_tmp;
                pList[cnt].number = cnt;
                pList[cnt].vTh = vTh;

                pList[cnt].d3v = dv;

                // Get vPar, vPer and mu for guiding center integration

                float bMag = mag (b0_XYZ);
                float vMag = mag (thisV_XYZ);

#if DEBUG_MAXWELLIAN >=2 
                    cout<<"ThisVx: "<<thisvx<<endl;
                    cout<<"ThisVy: "<<thisvy<<endl;
                    cout<<"ThisVz: "<<thisvz<<endl;
                    cout<<"vMag: "<<vMag<<endl;
                    cout<<"vPer: "<<pList[cnt].vPer<<endl;
                    cout<<"vPar: "<<pList[cnt].vPar<<endl;
                    cout<<"u: "<<pList[cnt].u<<endl<<endl;
                    if(isnan(pList[cnt].u)) exit(1);
#endif                    
                cnt++;
        }

#if DEBUG_MAXWELLIAN >=1 
        cout << "TestIntegratedValue: " << TestIntegratedValue << endl;
#endif
        return pList;
}




// Calculate the jP given some known E and f(v)

int main ( int argc, char **argv )
{

        // Make sure the "output/" directory exists
		stringstream outputDirName;
		outputDirName << "output/";
		// check directory exists
		struct stat st;
        int dirTest = stat(outputDirName.str().c_str(),&st);
		if( dirTest != 0 ) {
            cout<<"Had to create output/ directory"<<endl;
			int mkDirStat = mkdir(outputDirName.str().c_str(),S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
		}
	
#if CLOCK >= 1
        clock_t ProgramTime = clock();
#endif
	
#if (USEPAPI >= 1 || LOWMEM_USEPAPI >=1 )
		float realTime0, cpuTime0, realTime=0, cpuTime=0, mFlops=0;
		long long flpIns0, flpIns=0;
		int papiReturn;

		cpuTime0=cpuTime;realTime0=realTime;flpIns0=flpIns;
		papiReturn = PAPI_flops ( &realTime, &cpuTime, &flpIns, &mFlops );
        if(papiReturn<0) {
                cout<<"ERROR: PAPI Failed to initialize with error code: "<<papiReturn<<endl;
                cout<<"ERROR: See papi.h for error code explanations "<<endl;
                exit(1);
        }
		printf("Real_time:\t%f\nProc_time:\t%f\nTotal flpins:\t%lld\nMFLOPS:\t\t%f\n",
					   realTime, cpuTime, flpIns, mFlops);

		papiReturn = PAPI_flops ( &realTime, &cpuTime, &flpIns, &mFlops );
        if(papiReturn<0) {
                cout<<"ERROR: PAPI Failed to initialize with error code: "<<papiReturn<<endl;
                cout<<"ERROR: See papi.h for error code explanations "<<endl;
                exit(1);
        } else {
                cout<<"PAPI called successfully with return code: "<<papiReturn<<endl;
        }
		printf("Real_time:\t%f\nProc_time:\t%f\nTotal flpins:\t%lld\nMFLOPS:\t\t%f\n",
					   realTime, cpuTime, flpIns, mFlops);
#endif

		libconfig::Config cfg;
		string cfgName = "kj.cfg";

		// Write a config file if required
		//libconfig::Setting &root = cfg.getRoot();
		//root.add("xGridMin", libconfig::Setting::TypeFloat) = 98.0;
		//root.add("xGridMax", libconfig::Setting::TypeFloat) = 100.0;
		//root.add("nXGrid", libconfig::Setting::TypeInt) = 20;
		//root.add("nRFCycles", libconfig::Setting::TypeFloat) = 10;
		//root.add("nStepsPerCycle", libconfig::Setting::TypeInt) = 100;
		//root.add("nJpCycles", libconfig::Setting::TypeInt) = 6;
		//root.add("nJpPerCycle", libconfig::Setting::TypeInt) = 20;
		//root.add("eField_fName", libconfig::Setting::TypeString) = "data/kj_aorsa_1d.nc";
		//root.add("particleList_fName", libconfig::Setting::TypeString) = "data/f.nc";
		//cfg.writeFile(cfgName.c_str());
		
		// Open the config file
		cfg.readFile(cfgName.c_str());

	    int species_number = cfg.lookup("species_number");

		// Read E
		string eField_fName = cfg.lookup("eField_fName");	
		cout << "Reading eField data file " << eField_fName << endl;

		// Here we are using the cxx-4 netcdf interface by Lynton Appel
		// This needs netCDF 4.1.1 or later build with
		// ./configure --enable-cxx-4 [plus other options]

		
        //// These are defined for DIM = 1 or 2
    
		float freq;

		ifstream file(eField_fName.c_str());
		if(!file.good()) {
			cout << "ERROR: Cannot find file " << eField_fName << endl;
			exit(1);
		}


#if DIM == 1
            vector<C3Vec> b0_CYL, b0_XYZ;

            vector<float>           r, b0_r, b0_p, b0_z,
                                    e_r_re, e_p_re, e_z_re,
                                    e_r_im, e_p_im, e_z_im, n_m3,
                                    b_r_re, b_p_re, b_z_re,
                                    b_r_im, b_p_im, b_z_im;

            vector<complex<float> >  e_r, e_p, e_z;
            vector<complex<float> >  b_r, b_p, b_z;

		try {
				NcFile dataFile ( eField_fName.c_str(), NcFile::read );
	
				NcDim nc_nR(dataFile.getDim("nR"));
				NcDim nc_nSpec(dataFile.getDim("nSpec"));
				NcDim nc_scalar(dataFile.getDim("scalar"));
            
                int nR = nc_nR.getSize();
                int nZ = 1;
                int nSpec = nc_nSpec.getSize();

                if(species_number>nSpec-1) {
                        cout << "ERROR: Asking for species that does not exist in density data" << endl;
                        exit(1);
                }

                cout << "\tnR: " << nR << endl;

                NcVar nc_r(dataFile.getVar("r"));
            
                NcVar nc_freq(dataFile.getVar("freq"));

                NcVar nc_b0_r(dataFile.getVar("B0_r"));
                NcVar nc_b0_p(dataFile.getVar("B0_p"));
                NcVar nc_b0_z(dataFile.getVar("B0_z"));

                NcVar nc_e_r_re(dataFile.getVar("e_r_re"));
                NcVar nc_e_p_re(dataFile.getVar("e_p_re"));
                NcVar nc_e_z_re(dataFile.getVar("e_z_re"));
                NcVar nc_e_r_im(dataFile.getVar("e_r_im"));
                NcVar nc_e_p_im(dataFile.getVar("e_p_im"));
                NcVar nc_e_z_im(dataFile.getVar("e_z_im"));

                NcVar nc_b_r_re(dataFile.getVar("b_r_re"));
                NcVar nc_b_p_re(dataFile.getVar("b_p_re"));
                NcVar nc_b_z_re(dataFile.getVar("b_z_re"));
                NcVar nc_b_r_im(dataFile.getVar("b_r_im"));
                NcVar nc_b_p_im(dataFile.getVar("b_p_im"));
                NcVar nc_b_z_im(dataFile.getVar("b_z_im"));

                NcVar nc_density(dataFile.getVar("density_m3"));

                r.resize(nR);
            
                b0_r.resize(nR);
                b0_p.resize(nR);
                b0_z.resize(nR);

                e_r_re.resize(nR);
                e_p_re.resize(nR);
                e_z_re.resize(nR);
                e_r_im.resize(nR);
                e_p_im.resize(nR);
                e_z_im.resize(nR);

                b_r_re.resize(nR);
                b_p_re.resize(nR);
                b_z_re.resize(nR);
                b_r_im.resize(nR);
                b_p_im.resize(nR);
                b_z_im.resize(nR);

                n_m3.resize(nR);

                nc_r.getVar(&r[0]);
                nc_freq.getVar(&freq);

                nc_b0_r.getVar(&b0_r[0]);
                nc_b0_p.getVar(&b0_p[0]);
                nc_b0_z.getVar(&b0_z[0]);

                // Here im reading a single species' density from a multi species array, 
                // i.e., density[nSpec,nR] and I only want density[1,*] for example where
                // the species is specified by "species_number" in the cfg file
                vector<size_t> start, count;
                start.resize(2);
                count.resize(2);
                start[1] = 0;
                start[0] = species_number;
                count[1] = nR;
                count[0] = 1;

				nc_density.getVar(start, count, &n_m3[0]);

				b0_CYL.resize(nR);
				b0_XYZ.resize(nR);
           
				for(int i=0; i<nR; i++) {
                        b0_CYL[i] = C3Vec(b0_r[i],b0_p[i],b0_z[i]);
                        b0_XYZ[i] = rot_CYL_to_XYZ(0,b0_CYL[i],1);
				}

				nc_e_r_re.getVar(&e_r_re[0]);
				nc_e_p_re.getVar(&e_p_re[0]);
				nc_e_z_re.getVar(&e_z_re[0]);
				nc_e_r_im.getVar(&e_r_im[0]);
				nc_e_p_im.getVar(&e_p_im[0]);
				nc_e_z_im.getVar(&e_z_im[0]);

				nc_b_r_re.getVar(&b_r_re[0]);
				nc_b_p_re.getVar(&b_p_re[0]);
				nc_b_z_re.getVar(&b_z_re[0]);
				nc_b_r_im.getVar(&b_r_im[0]);
				nc_b_p_im.getVar(&b_p_im[0]);
				nc_b_z_im.getVar(&b_z_im[0]);

				for(int i=0; i<nR; i++){
						e_r.push_back(complex<float>( e_r_re[i], e_r_im[i] ) );
						e_p.push_back(complex<float>( e_p_re[i], e_p_im[i] ) );
						e_z.push_back(complex<float>( e_z_re[i], e_z_im[i] ) );
				}

				for(int i=0; i<nR; i++){
						b_r.push_back(complex<float>( b_r_re[i], b_r_im[i] ) );
						b_p.push_back(complex<float>( b_p_re[i], b_p_im[i] ) );
						b_z.push_back(complex<float>( b_z_re[i], b_z_im[i] ) );
				}

				vector<float>::iterator min = min_element(b0_p.begin(),b0_p.end());
				vector<float>::iterator max = max_element(b0_p.begin(),b0_p.end());
#if DEBUGLEVEL >= 1 
				cout << "\tR[0]: " << r[0] << ", R["<<nR<<"]: " << r[r.size()-1] << endl;
				cout << "\tfreq: " << freq << endl;
				cout << "\tmin(b0_p): " << *min << endl;
				cout << "\tmax(b0_p): " << *max << endl;
				cout << "\tabs(e_r[nR/2]): " << abs(e_r[nR/2]) << endl;
				cout << "\tabs(e_p[nR/2]): " << abs(e_p[nR/2]) << endl;
				cout << "\tabs(e_z[nR/2]): " << abs(e_z[nR/2]) << endl;
#endif
		}
		catch(exceptions::NcException &e) {
				cout << "NetCDF: unknown error" << endl;
				e.what();
				exit(1);
		}
    
        fieldMeshClass fieldMesh(r);

		// Read the guiding center terms from file
		string gc_fName = cfg.lookup("gc_fName");	
		cout << "Reading GC terms data file " << gc_fName << endl;

		vector<float> r_gc, curv_r, curv_p, curv_z,
            grad_r, grad_p, grad_z, bDotGradB;
		vector<C3Vec> curv_CYL, grad_CYL;
		
		ifstream gc_file(gc_fName.c_str());
		if(!gc_file.good()) {
			cout << "ERROR: Cannot find file " << gc_fName << endl;
			exit(1);
		}

		try {
				NcFile dataFile ( gc_fName.c_str(), NcFile::read );
	
				NcDim gc_nc_nR(dataFile.getDim("nR"));
				NcDim gc_nc_scalar(dataFile.getDim("scalar"));
	
				int nR_gc = gc_nc_nR.getSize();
                cout << "nR_gc: " << nR_gc << endl;

				NcVar gc_nc_r(dataFile.getVar("r"));

				NcVar gc_nc_curv_r(dataFile.getVar("curv_r"));
				NcVar gc_nc_curv_p(dataFile.getVar("curv_t"));
				NcVar gc_nc_curv_z(dataFile.getVar("curv_z"));

				NcVar gc_nc_grad_r(dataFile.getVar("grad_r"));
				NcVar gc_nc_grad_p(dataFile.getVar("grad_t"));
				NcVar gc_nc_grad_z(dataFile.getVar("grad_z"));

				NcVar gc_nc_bDotGradB(dataFile.getVar("bDotGradB"));

				r_gc.resize(nR_gc);

				curv_r.resize(nR_gc);
				curv_p.resize(nR_gc);
				curv_z.resize(nR_gc);

				grad_r.resize(nR_gc);
				grad_p.resize(nR_gc);
				grad_z.resize(nR_gc);

				bDotGradB.resize(nR_gc);

				gc_nc_r.getVar(&r_gc[0]);

				gc_nc_curv_r.getVar(&curv_r[0]);
				gc_nc_curv_p.getVar(&curv_p[0]);
				gc_nc_curv_z.getVar(&curv_z[0]);

				gc_nc_grad_r.getVar(&grad_r[0]);
				gc_nc_grad_p.getVar(&grad_p[0]);
				gc_nc_grad_z.getVar(&grad_z[0]);

                gc_nc_bDotGradB.getVar(&bDotGradB[0]);

				curv_CYL.resize(nR_gc);
				grad_CYL.resize(nR_gc);
				for(int i=0; i<nR_gc; i++) {
						curv_CYL[i] = C3Vec(curv_r[i],curv_p[i],curv_z[i]);
		                grad_CYL[i] = C3Vec(grad_r[i],grad_p[i],grad_z[i]);
				}

		}
		catch(exceptions::NcException &e) {
				cout << "NetCDF: unknown error." << endl;
				e.what();
				exit(1);
		}

		// Rotate the e & b fields to XYZ

		vector<C3Vec> e1Re_CYL, e1Im_CYL, b1Re_CYL, b1Im_CYL;
		e1Re_CYL.resize(e_r.size());
		e1Im_CYL.resize(e_r.size());
		b1Re_CYL.resize(e_r.size());
		b1Im_CYL.resize(e_r.size());

		vector<C3Vec> e1Re_XYZ, e1Im_XYZ, b1Re_XYZ, b1Im_XYZ;
		e1Re_XYZ.resize(e_r.size());
		e1Im_XYZ.resize(e_r.size());
		b1Re_XYZ.resize(e_r.size());
		b1Im_XYZ.resize(e_r.size());

		for(int i=0;i<e_r.size();i++) {

			e1Re_CYL[i].c1 = real(e_r[i]);
			e1Re_CYL[i].c2 = real(e_p[i]);
			e1Re_CYL[i].c3 = real(e_z[i]);
			e1Im_CYL[i].c1 = imag(e_r[i]);
			e1Im_CYL[i].c2 = imag(e_p[i]);
			e1Im_CYL[i].c3 = imag(e_z[i]);

			b1Re_CYL[i].c1 = real(b_r[i]);
			b1Re_CYL[i].c2 = real(b_p[i]);
			b1Re_CYL[i].c3 = real(b_z[i]);
			b1Im_CYL[i].c1 = imag(b_r[i]);
			b1Im_CYL[i].c2 = imag(b_p[i]);
			b1Im_CYL[i].c3 = imag(b_z[i]);

			float _p = 0;

            e1Re_XYZ[i] = rot_CYL_to_XYZ ( _p, e1Re_CYL[i], 1);
            e1Im_XYZ[i] = rot_CYL_to_XYZ ( _p, e1Im_CYL[i], 1);
            
            b1Re_XYZ[i] = rot_CYL_to_XYZ ( _p, b1Re_CYL[i], 1);
            b1Im_XYZ[i] = rot_CYL_to_XYZ ( _p, b1Im_CYL[i], 1);

		}
#endif

#if DIM == 2

            vector< vector<C3Vec> > b0_CYL, b0_XYZ;
            vector<float> r;
            vector<float> z;
            vector<vector<float> >    b0_r, b0_p, b0_z,
                                    e_r_re, e_p_re, e_z_re,
                                    e_r_im, e_p_im, e_z_im, n_m3,
                                    b_r_re, b_p_re, b_z_re,
                                    b_r_im, b_p_im, b_z_im;

            vector< vector<complex<float> > > e_r, e_p, e_z;
            vector <vector<complex<float> > > b_r, b_p, b_z;

		try {
				NcFile dataFile ( eField_fName.c_str(), NcFile::read );

	
				NcDim nc_nR(dataFile.getDim("nR"));
				NcDim nc_nSpec(dataFile.getDim("nSpec"));
				NcDim nc_scalar(dataFile.getDim("scalar"));

                NcDim nc_nZ(dataFile.getDim("nZ"));
    
                int nR = nc_nR.getSize();
                int nZ = nc_nZ.getSize();
                int nSpec = nc_nSpec.getSize();

                if(species_number>nSpec-1) {
                        cout << "ERROR: Asking for species that does not exist in density data" << endl;
                        exit(1);
                }

                cout << "\tnR: " << nR << endl;

                NcVar nc_r(dataFile.getVar("r"));
                NcVar nc_z(dataFile.getVar("z"));
            
                NcVar nc_freq(dataFile.getVar("freq"));

                NcVar nc_b0_r(dataFile.getVar("B0_r"));
                NcVar nc_b0_p(dataFile.getVar("B0_p"));
                NcVar nc_b0_z(dataFile.getVar("B0_z"));

                NcVar nc_e_r_re(dataFile.getVar("e_r_re"));
                NcVar nc_e_p_re(dataFile.getVar("e_p_re"));
                NcVar nc_e_z_re(dataFile.getVar("e_z_re"));
                NcVar nc_e_r_im(dataFile.getVar("e_r_im"));
                NcVar nc_e_p_im(dataFile.getVar("e_p_im"));
                NcVar nc_e_z_im(dataFile.getVar("e_z_im"));

                NcVar nc_b_r_re(dataFile.getVar("b_r_re"));
                NcVar nc_b_p_re(dataFile.getVar("b_p_re"));
                NcVar nc_b_z_re(dataFile.getVar("b_z_re"));
                NcVar nc_b_r_im(dataFile.getVar("b_r_im"));
                NcVar nc_b_p_im(dataFile.getVar("b_p_im"));
                NcVar nc_b_z_im(dataFile.getVar("b_z_im"));

                NcVar nc_density(dataFile.getVar("density_m3"));

                r.resize(nR);
                z.resize(nZ);

				for(int i=0; i<nR; i++)
				{
					vector<float> row(nZ);
                    vector<C3Vec> rowC3Vec(nZ);
                    
					b0_r.push_back(row);
					b0_p.push_back(row);
					b0_z.push_back(row);

                    e_r_re.push_back(row);
                    e_p_re.push_back(row);
                    e_z_re.push_back(row);
                    e_r_im.push_back(row);
                    e_p_im.push_back(row);
                    e_z_im.push_back(row);

                    b_r_re.push_back(row);
                    b_p_re.push_back(row);
                    b_z_re.push_back(row);
                    b_r_im.push_back(row);
                    b_p_im.push_back(row);
                    b_z_im.push_back(row);
                    
                    n_m3.push_back(row);

                    b0_CYL.push_back(rowC3Vec);
                    b0_XYZ.push_back(rowC3Vec);
				}
 
                nc_r.getVar(&r[0]);
                nc_z.getVar(&z[0]);

                nc_freq.getVar(&freq);
            
                float* b0_r_2D = new float[nR*nZ];
				float* b0_p_2D = new float[nR*nZ];
				float* b0_z_2D = new float[nR*nZ];

                float* e_r_re_2D = new float[nR*nZ];
				float* e_p_re_2D = new float[nR*nZ];
				float* e_z_re_2D = new float[nR*nZ];
                float* e_r_im_2D = new float[nR*nZ];
				float* e_p_im_2D = new float[nR*nZ];
				float* e_z_im_2D = new float[nR*nZ];

                float* b_r_re_2D = new float[nR*nZ];
				float* b_p_re_2D = new float[nR*nZ];
				float* b_z_re_2D = new float[nR*nZ];
                float* b_r_im_2D = new float[nR*nZ];
				float* b_p_im_2D = new float[nR*nZ];
				float* b_z_im_2D = new float[nR*nZ];

				float* n_m3_2D = new float[nR*nZ*nSpec];
            
                nc_b0_r.getVar(b0_r_2D);
                nc_b0_p.getVar(b0_p_2D);
                nc_b0_z.getVar(b0_z_2D);

				nc_e_r_re.getVar(e_r_re_2D);
				nc_e_p_re.getVar(e_p_re_2D);
				nc_e_z_re.getVar(e_z_re_2D);
				nc_e_r_im.getVar(e_r_im_2D);
				nc_e_p_im.getVar(e_p_im_2D);
				nc_e_z_im.getVar(e_z_im_2D);

				nc_b_r_re.getVar(b_r_re_2D);
				nc_b_p_re.getVar(b_p_re_2D);
				nc_b_z_re.getVar(b_z_re_2D);
				nc_b_r_im.getVar(b_r_im_2D);
				nc_b_p_im.getVar(b_p_im_2D);
				nc_b_z_im.getVar(b_z_im_2D);
            
				for(int i=0; i<nR; i++) {
                    for (int j = 0; j< nZ; j++){
						b0_r[i][j] = b0_r_2D[i*nZ+j];
						b0_p[i][j] = b0_p_2D[i*nZ+j];
						b0_z[i][j] = b0_z_2D[i*nZ+j];

						e_r_re[i][j] = e_r_re_2D[i*nZ+j];
						e_p_re[i][j] = e_p_re_2D[i*nZ+j];
						e_z_re[i][j] = e_z_re_2D[i*nZ+j];
						e_r_im[i][j] = e_r_im_2D[i*nZ+j];
						e_p_im[i][j] = e_p_im_2D[i*nZ+j];
						e_z_im[i][j] = e_z_im_2D[i*nZ+j];

						b_r_re[i][j] = b_r_re_2D[i*nZ+j];
						b_p_re[i][j] = b_p_re_2D[i*nZ+j];
						b_z_re[i][j] = b_z_re_2D[i*nZ+j];
						b_r_im[i][j] = b_r_im_2D[i*nZ+j];
						b_p_im[i][j] = b_p_im_2D[i*nZ+j];
						b_z_im[i][j] = b_z_im_2D[i*nZ+j];
                        
                    }
                }
            
                cout << "hello, got vars from netcdf field and deleted arrays" << endl;
                // Here im reading a single species' density from a multi species array,
                // i.e., density[nSpec,nR] and I only want density[1,*] for example where
                // the species is specified by "species_number" in the cfg file
            
                vector<size_t> start, count;
                start.resize(2);
                count.resize(2);
                start[1] = 0;
                start[0] = species_number;
                count[1] = nR*nZ;
                count[0] = 1;
            
                cout << "species number        " << species_number << endl;
                cout << "about to get n_m3     " << endl;
                //nc_density.getVar(start, count, n_m3_2D);
            
                nc_density.getVar(n_m3_2D);
            
                cout << "got n_m3" << endl;

 				for(int i=0; i<nR; i++) {
                    for (int j = 0; j< nZ; j++){
                        n_m3[i][j] = n_m3_2D[i*nZ + j];
                    }
                //cout << "n_m3[i][j]    " << n_m3[i][0] << endl;

                }

                cout << "got and stored n_m3" << endl;

				for(int i=0; i<nR; i++) {
                    for (int j = 0; j< nZ; j++){
                        b0_CYL[i][j] = C3Vec(b0_r[i][j],b0_p[i][j],b0_z[i][j]);
                        b0_XYZ[i][j] = rot_CYL_to_XYZ(0,b0_CYL[i][j],1);
                        }
				}

                delete [] b0_r_2D;
				delete [] b0_p_2D;
				delete [] b0_z_2D;

                delete [] e_r_re_2D;
				delete [] e_p_re_2D;
				delete [] e_z_re_2D;
                delete [] e_r_im_2D;
				delete [] e_p_im_2D;
				delete [] e_z_im_2D;

                delete [] b_r_re_2D;
				delete [] b_p_re_2D;
				delete [] b_z_re_2D;
                delete [] b_r_im_2D;
				delete [] b_p_im_2D;
				delete [] b_z_im_2D;
            
                vector<complex<float> > e_r_tmp, e_p_tmp, e_z_tmp;

				for(int i=0; i<nR; i++){
                    for(int j=0; j<nZ; j++){
                    
                        e_r_tmp.push_back(complex<float>( e_r_re[i][j], e_r_im[i][j] ) );
                        e_p_tmp.push_back(complex<float>( e_r_re[i][j], e_r_im[i][j] ) );
                        e_z_tmp.push_back(complex<float>( e_r_re[i][j], e_r_im[i][j] ) );
                        
                    }
                    
                    e_r.push_back(e_r_tmp);
                    e_p.push_back(e_r_tmp);
                    e_z.push_back(e_r_tmp);
                    
				}

                vector<complex<float> > b_r_tmp, b_p_tmp, b_z_tmp;
				for(int i=0; i<nR; i++){
                    for(int j=0; j<nZ; j++){
                        b_r_tmp.push_back(complex<float>( b_r_re[i][j], b_r_im[i][j] ) );
                        b_p_tmp.push_back(complex<float>( b_r_re[i][j], b_r_im[i][j] ) );
                        b_z_tmp.push_back(complex<float>( b_r_re[i][j], b_r_im[i][j] ) );
                        
                    }
            
                    b_r.push_back(b_r_tmp);
                    b_p.push_back(b_r_tmp);
                    b_z.push_back(b_r_tmp);
                    
				}
//////////////////////////////  how to change to 2D? Commented out iterators, do not seemed to be used.
            //	vector<float>::iterator min = min_element(b0_p.begin(),b0_p.end());
			//	vector<float>::iterator max = max_element(b0_p.begin(),b0_p.end());
            
//////////////// guiding center data not changed to 2D...............
#if DEBUGLEVEL >= 1 
				cout << "\tR[0]: " << r[0] << ", R["<<nR<<"]: " << r[r.size()-1] << endl;
				cout << "\tfreq: " << freq << endl;
				cout << "\tmin(b0_p): " << *min << endl;
				cout << "\tmax(b0_p): " << *max << endl;
				cout << "\tabs(e_r[nR/2]): " << abs(e_r[nR/2]) << endl;
				cout << "\tabs(e_p[nR/2]): " << abs(e_p[nR/2]) << endl;
				cout << "\tabs(e_z[nR/2]): " << abs(e_z[nR/2]) << endl;
#endif
		}
		catch(exceptions::NcException &e) {
				cout << "NetCDF: unknown error" << endl;
				e.what();
				exit(1);
		}
        fieldMeshClass fieldMesh(r,z);
    
    // make sure z is monotonically increasing
        sort(fieldMesh.z.begin(), fieldMesh.z.end());

		// Read the guiding center terms from file
		string gc_fName = cfg.lookup("gc_fName");	
		cout << "Reading GC terms data file " << gc_fName << endl;

		vector<float> r_gc, curv_r, curv_p, curv_z,
            grad_r, grad_p, grad_z, bDotGradB;
		vector<C3Vec> curv_CYL, grad_CYL;
		
		ifstream gc_file(gc_fName.c_str());
		if(!gc_file.good()) {
			cout << "ERROR: Cannot find file " << gc_fName << endl;
			exit(1);
		}

		try {
				NcFile dataFile ( gc_fName.c_str(), NcFile::read );
	
				NcDim gc_nc_nR(dataFile.getDim("nR"));
				NcDim gc_nc_scalar(dataFile.getDim("scalar"));
	
				int nR_gc = gc_nc_nR.getSize();
                cout << "nR_gc: " << nR_gc << endl;

				NcVar gc_nc_r(dataFile.getVar("r"));

				NcVar gc_nc_curv_r(dataFile.getVar("curv_r"));
				NcVar gc_nc_curv_p(dataFile.getVar("curv_t"));
				NcVar gc_nc_curv_z(dataFile.getVar("curv_z"));

				NcVar gc_nc_grad_r(dataFile.getVar("grad_r"));
				NcVar gc_nc_grad_p(dataFile.getVar("grad_t"));
				NcVar gc_nc_grad_z(dataFile.getVar("grad_z"));

				NcVar gc_nc_bDotGradB(dataFile.getVar("bDotGradB"));

				r_gc.resize(nR_gc);

				curv_r.resize(nR_gc);
				curv_p.resize(nR_gc);
				curv_z.resize(nR_gc);

				grad_r.resize(nR_gc);
				grad_p.resize(nR_gc);
				grad_z.resize(nR_gc);

				bDotGradB.resize(nR_gc);

				gc_nc_r.getVar(&r_gc[0]);

				gc_nc_curv_r.getVar(&curv_r[0]);
				gc_nc_curv_p.getVar(&curv_p[0]);
				gc_nc_curv_z.getVar(&curv_z[0]);

				gc_nc_grad_r.getVar(&grad_r[0]);
				gc_nc_grad_p.getVar(&grad_p[0]);
				gc_nc_grad_z.getVar(&grad_z[0]);

                gc_nc_bDotGradB.getVar(&bDotGradB[0]);

				curv_CYL.resize(nR_gc);
				grad_CYL.resize(nR_gc);
				for(int i=0; i<nR_gc; i++) {
						curv_CYL[i] = C3Vec(curv_r[i],curv_p[i],curv_z[i]);
		                grad_CYL[i] = C3Vec(grad_r[i],grad_p[i],grad_z[i]);
				}

		}
		catch(exceptions::NcException &e) {
				cout << "NetCDF: unknown error." << endl;
				e.what();
				exit(1);
		}


		// Rotate the e & b fields to XYZ
		vector< vector<C3Vec> > e1Re_CYL, e1Im_CYL, b1Re_CYL, b1Im_CYL;
        vector< vector<C3Vec> > e1Re_XYZ, e1Im_XYZ, b1Re_XYZ, b1Im_XYZ;
    
		for (int i = 0; i<e_r.size(); i++){
            vector<C3Vec> rowC3Vec(e_r[0].size());
            
            b0_CYL.push_back(rowC3Vec);
            
            e1Re_CYL.push_back(rowC3Vec);
            e1Im_CYL.push_back(rowC3Vec);
            b1Re_CYL.push_back(rowC3Vec);
            b1Im_CYL.push_back(rowC3Vec);

            e1Re_XYZ.push_back(rowC3Vec);
            e1Im_XYZ.push_back(rowC3Vec);
            b1Re_XYZ.push_back(rowC3Vec);
            b1Im_XYZ.push_back(rowC3Vec);
        }

		for(int i=0;i<e_r.size();i++) {
            for(int j=0;j<e_r[0].size();j++) {
			e1Re_CYL[i][j].c1 = real(e_r[i][j]);
			e1Re_CYL[i][j].c2 = real(e_p[i][j]);
			e1Re_CYL[i][j].c3 = real(e_z[i][j]);
			e1Im_CYL[i][j].c1 = imag(e_r[i][j]);
			e1Im_CYL[i][j].c2 = imag(e_p[i][j]);
			e1Im_CYL[i][j].c3 = imag(e_z[i][j]);

			b1Re_CYL[i][j].c1 = real(b_r[i][j]);
			b1Re_CYL[i][j].c2 = real(b_p[i][j]);
			b1Re_CYL[i][j].c3 = real(b_z[i][j]);
			b1Im_CYL[i][j].c1 = imag(b_r[i][j]);
			b1Im_CYL[i][j].c2 = imag(b_p[i][j]);
			b1Im_CYL[i][j].c3 = imag(b_z[i][j]);

			float _p = 0;

            e1Re_XYZ[i][j] = rot_CYL_to_XYZ ( _p, e1Re_CYL[i][j], 1);
            e1Im_XYZ[i][j] = rot_CYL_to_XYZ ( _p, e1Im_CYL[i][j], 1);
            
            b1Re_XYZ[i][j] = rot_CYL_to_XYZ ( _p, b1Re_CYL[i][j], 1);
            b1Im_XYZ[i][j] = rot_CYL_to_XYZ ( _p, b1Im_CYL[i][j], 1);
            }
		}

#endif
    
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////// NOTHING BELOW HERE SEEMS TO DEPEND ON WHETHER DOING A 1D OR A 2D CALCULATION
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

	float wrf = freq * 2 * _pi;

    cout << "Loaded all netcdf data!!! " << endl;

  //  cout << " x  " << x << "z  " << z << " fieldMesh.r.front  " << endl;
    
	//string googlePerfFileName = "/home/dg6/code/kineticj/googlep";
	//ProfilerStart(googlePerfFileName.c_str());
    
#if USEPAPI >= 1
	cpuTime0=cpuTime;realTime0=realTime;flpIns0=flpIns;
	papiReturn = PAPI_flops ( &realTime, &cpuTime, &flpIns, &mFlops );
	printf("\nStartup performance:\n");
	printf("Real_time:\t%f\nProc_time:\t%f\nTotal flpins:\t%lld\nMFLOPS:\t\t%f\n",
					   realTime-realTime0, cpuTime-cpuTime0, flpIns-flpIns0, mFlops);
#endif

	float nRFCycles 		= cfg.lookup("nRFCycles");
	float nStepsPerCycle 	= cfg.lookup("nStepsPerCycle"); 
	float tRF 			= (2*_pi)/wrf;
	int nJpCycles 		= cfg.lookup("nJpCycles");
	int nJpPerCycle 	= cfg.lookup("nJpPerCycle");
    int nPhi            = cfg.lookup("nPhi");
	int nJp 			= nJpCycles * nJpPerCycle;
	float dtJp 			= tRF / nJpPerCycle;
	int istat = 0;
    int nPx = cfg.lookup("nPx");
    int nPy = cfg.lookup("nPy");
    int nPz = cfg.lookup("nPz");
    int nPVpar = cfg.lookup("nPVpar");
    int nPVper = cfg.lookup("nPVper");
    
    int nPGyro = cfg.lookup("nPGyro");
    int nPTor =  cfg.lookup("nPTor");
    
    float amu = cfg.lookup("species_amu");
    float Z = cfg.lookup("species_Z");
    int nThermal = cfg.lookup("nThermal");
	int nP;
    //= nPx*nPy*nPz;
    int nPblob = cfg.lookup("nPblob");
    int nSavePerRFCycle = cfg.lookup("nSavePerRFCycle");
    int nList;
    
#if defined(_OPENMP)
        int nThreads;
#endif
        float dv;

//////////// Initialize primary phase space worklist:
	float rGridMin = cfg.lookup("rGridMin");
	float rGridMax = cfg.lookup("rGridMax");

	float zGridMin = cfg.lookup("zGridMin");
	float zGridMax = cfg.lookup("zGridMax");

	float phiGridMin = cfg.lookup("phiGridMin");
	float phiGridMax = cfg.lookup("phiGridMax");
    
	int nRGrid = cfg.lookup("nRGrid");
	int nZGrid = cfg.lookup("nZGrid");
	int nPhiGrid = cfg.lookup("nPhiGrid");

	int nVparGrid = cfg.lookup("nVparGrid");
	int nVperGrid = cfg.lookup("nVperGrid");

	int nPhaseGrid = cfg.lookup("nPhaseGrid");

    nList = nRGrid*nVparGrid*nVperGrid*nPhaseGrid;
    cout << "nList     " << nList << endl;

	vector<float> rGrid(nRGrid), zGrid(nZGrid), phiGrid(nPhiGrid);
    vector<float> density_m3(nList), T_keV(nList), wrf_wc(nList);
    
    vector<C3Vec> b0_XYZ_T_at_List(nList);
    
    vector<float> VparGrid(nVparGrid), VperGrid(nVperGrid);
    vector<float> PhaseGrid(nPhaseGrid);
    
    vector<CParticle> PrimaryWorkList;
    PrimaryWorkList.resize(nList);
    
    vector<float> Dvperp(nList), Dvpar(nList);
    
    float rGridRng = 0;
	float rGridStep = 0;

    float zGridRng = 0;
	float zGridStep = 0;

    float phiGridRng = 0;
	float phiGridStep = 0;
	
	if(nRGrid>1) {
		rGridRng = rGridMax-rGridMin;
		rGridStep = rGridRng/(nRGrid-1);
	}

	if(nZGrid>1) {
		zGridRng = zGridMax-zGridMin;
		zGridStep = zGridRng/(nZGrid-1);
	}

	if(nPhiGrid>1) {
		phiGridRng = phiGridMax-phiGridMin;
		phiGridStep = phiGridRng/(nPhiGrid-1);
	}

    for (int iList=0;iList<nList;iList++){
        T_keV[iList] = 1.0;
    }
    
    float vth_part = get_vTh(amu, Z,T_keV[0]);

    /// initialize (rGrid, VparGrid, VperGrid)
	for(int iR=0;iR<nRGrid;iR++) {
		rGrid[iR] = rGridMin+iR*rGridStep;
    }

	for(int iPhi=0;iPhi<nPhiGrid;iPhi++) {
		phiGrid[iPhi] = phiGridMin+iPhi*phiGridStep;
    }

	for(int iZ=0;iZ<nZGrid;iZ++) {
		zGrid[iZ] = zGridMin+iZ*zGridStep;
    }
    
	for(int iVpar=0;iVpar<nVparGrid;iVpar++) {
		VparGrid[iVpar] = 0.3*vth_part + iVpar*(vth_part/nVparGrid);
    }
    
	for(int iVper=0;iVper<nVperGrid;iVper++) {
		VperGrid[iVper] = 0.3*vth_part + iVper*(vth_part/nVperGrid);
    }
    
    for(int iPhase=0;iPhase<nPhaseGrid;iPhase++){
        PhaseGrid[iPhase] = iPhase*2.0*_pi/nPhaseGrid;
    }
    
    // initialize primary worklist (assuming y,z,gyrophase = 0)
    for (int iList=0;iList<nList;iList++){
            int iR = iList % nRGrid;
            int iPhi = int(floor(iList/nRGrid)) % nPhiGrid;
            int iZ = int(floor(iList/nRGrid*nPhiGrid)) % nZGrid;
        
            int iVpar = int(floor( iList/nRGrid*nZGrid*nPhiGrid)) % nVparGrid;
            int iVper = int(floor(iList/(nRGrid*nZGrid*nPhiGrid*nVparGrid))) % nVperGrid;
            int iPhase = int(floor(iList/(nRGrid*nZGrid*nPhiGrid*nVparGrid*nVperGrid)));
        
            cout << "iList   = " << iList << "    iPhase = " << iPhase << "   PhaseGrid[iPhase] = " << PhaseGrid[iPhase]
            << ",          iR    =  "     << iR << "     " << "r[iR] = " << rGrid[iR]
            << ",          iZ    =  "     << iZ << "     " << "z[iZ] = " << zGrid[iZ]
            << ",          iPhi    =  "     << iPhi << "     " << "phi[iPhi] = " << phiGrid[iPhi]
            << ",          iVpar   =  "     << iR << "     " << "Vpar[iVpar] = " << VparGrid[iVpar]
            << ",          iVper    =  "     << iR << "     " << "Vper[iVper] = " << VperGrid[iVper] << endl;
        
           PrimaryWorkList[iList] = CParticle(rGrid[iR], phiGrid[iPhi],zGrid[iZ], VparGrid[iVpar],VperGrid[iVper],0.0,amu,Z,0.0, PhaseGrid[iPhase]);
    }

/*
    for (int i = 0; i< fieldMesh.r.size(); i++){
            cout << fieldMesh.r[i] << endl;
    }

    for (int i = 0; i< fieldMesh.z.size(); i++){
            cout << fieldMesh.z[i] << endl;
    }
*/
    
    // Initialize other variables on worklist, density, Bfield, etc
    for (int iList=0;iList<nList;iList++){
        int iStat;
        density_m3[iList] = kj_interp(C3Vec(PrimaryWorkList[iList].c1,PrimaryWorkList[iList].c2,PrimaryWorkList[iList].c3),fieldMesh,n_m3,iStat);
        cout << "density_m3 from kj_interp    " << density_m3[iList] << endl;
        b0_XYZ_T_at_List[iList] = kj_interp(C3Vec(PrimaryWorkList[iList].c1,PrimaryWorkList[iList].c2,PrimaryWorkList[iList].c3),fieldMesh,b0_XYZ,iStat);
        cout << "b0 from kj_interp    " << b0_XYZ_T_at_List[iList] << endl;
    }

    for (int iList=0;iList<nList;iList++){
    float this_wc =	Z*_e*mag(b0_XYZ_T_at_List[iList])/(amu*_mi);
		wrf_wc[iList] =  wrf / this_wc;
	}

//	for(int iR=0;iR<nRGrid;iR++) {
//	}

    float MaxB0 = maxC3VecAbs(b0_XYZ_T_at_List);
    float wc = Z*_e*MaxB0/(amu*_mi);
    float cyclotronPeriod = 2*_pi / wc;
	//float dtMin 	= -tRF/nStepsPerCycle;
    float dtMin 	= -cyclotronPeriod/nStepsPerCycle;
	int nSteps 	= nRFCycles*nStepsPerCycle;
	//int nSteps 		= nRFCycles*tRF/abs(dtMin)+1;
    bool SaveSingleOrbit = cfg.lookup("SaveSingleOribt");


#if PRINT_INFO >= 1
    cout << "dtMin [s]: " << dtMin << endl;
    cout << "Cyclotron Period: "<< cyclotronPeriod<<endl;
    cout << "RF Period: " << tRF << endl;
    cout << "nSteps: " << nSteps << endl;
    cout << "freq: " << freq << endl;
#endif
	
#if !(LOWMEM >= 1)
	vector<float> tJp(nJp,0);
	for(int jt=0;jt<nJp;jt++) {
		tJp[jt] = jt*dtJp;//-0.33*dtJp;
	}
#endif
	vector<float> thisT(nSteps);
	for(int i=0;i<nSteps;i++) {	
		thisT[i]=i*dtMin;//+1.5*dtMin;
	}

	vector<float> hanningWeight(nSteps);
	vector<float> expWeight(nSteps);
	vector<float> linearWeight(nSteps);
	for(int i=0;i<nSteps;i++) {
		//linearWeight[i]=thisT[i]*1.0/(tRF*nRFCycles)+1.0;
		hanningWeight[i]=0.5*(1-cos(2*_pi*i/(nSteps-1))); //Regular
		//hanningWeight[i]=0.5*(1-cos(2*_pi*i/(nSteps*0.25-1))); //Sharper

		//hanningWeight[i] = linearWeight[i];
		if(i<nSteps/2) hanningWeight[i]=1; //Regular
		//if(i<nSteps*7.0/8.0) hanningWeight[i]=1; //Sharper

		complex<float> _i (0.0,1.0);	
		complex<float> wrf_c (wrf,wrf*0.0025);
		expWeight[i] = abs(exp(-_i*wrf_c*thisT[i]));
		hanningWeight[i] = hanningWeight[i] * expWeight[i];
	}

	vector<vector<float> > j1x(nList), j1y(nList), j1z(nList);
#if LOWMEM >= 1
	vector<complex<float> > j1xc(nList), j1yc(nList), j1zc(nList);
#else
	vector<vector<complex<float> > >j1xc(nList), j1yc(nList), j1zc(nList);
#endif


	#pragma omp parallel for private(istat)
	for(int iList=0;iList<nList;iList++) {

#if defined(_OPENMP)
        nThreads = omp_get_num_threads();
#endif

#if !(LOWMEM >= 1)
		j1x[iList].resize(nJp);
		j1xc[iList].resize(nJp);

		j1y[iList].resize(nJp);
		j1yc[iList].resize(nJp);

		j1z[iList].resize(nJp);
		j1zc[iList].resize(nJp);
#endif
#if CLOCK >= 1
        clock_t startTime = clock();
#endif
#if LOWMEM >= 1
        j1xc[iList] = complex<float>(0,0);
        j1yc[iList] = complex<float>(0,0);
        j1zc[iList] = complex<float>(0,0);
#else
		for(int jt=0;jt<nJp;jt++) {
            j1xc[iList][jt] = complex<float>(0,0);
            j1yc[iList][jt] = complex<float>(0,0);
            j1zc[iList][jt] = complex<float>(0,0);
        }
#endif

#if LOWMEM_USEPAPI >= 1
		cpuTime0=cpuTime;realTime0=realTime;flpIns0=flpIns;
		papiReturn = PAPI_flops ( &realTime, &cpuTime, &flpIns, &mFlops );
#endif	

#if LOWMEM >= 1

 // START OF THE LOWMEM CODING vvv


//// kineticj: 1-D create nPx*nPy*nPZ phase space particles
//        vector<CParticle> ThisParticleList(create_particles(PrimaryWorkList[iList].c1,amu,Z,T_keV[iList],density_m3[iList],
//                                nPx,nPy,nPz,nThermal,dv,b0_XYZ_T_at_rGrid[iList]));

/////// dc: create blob, Larmor orbit only
//        nP = nPblob;
//        vector<CParticle> ThisParticleList(create_particle_blob(PrimaryWorkList[iList],amu,Z,T_keV[iList],density_m3[iList],
//                                nP,nThermal,dv,b0_XYZ_T_at_List[iList]));

/////// dc: create blob, e.g. gyrophase and toroidal average
        nP = nPGyro*nPTor;
        vector<CParticle> ThisParticleList(create_particle_blob(PrimaryWorkList[iList],amu,Z,T_keV[iList],density_m3[iList],
                                nPGyro, nPTor,nThermal,dv,b0_XYZ_T_at_List[iList]));
        
        cout << "nP       " << nP << endl;
        cout << "iList       " << iList << endl;

		vector<float> f1(nP);
		//vector<complex<float> > f1xc(nP), f1yc(nP), f1zc(nP);
		vector<complex<float> > f1c(nP);
        
/// Save particles to file to test
//
        if (iList == 0){
            ofstream ParticleFile;
            ParticleFile.open("output/particles.txt", ios::out | ios::trunc);
 //           ParticleFile << rGrid[iList] <<"    "<< 0.0 <<"    "<< 0.0 <<    <<     <<     << endl;
            
            for(int iP=0;iP <nP;iP++) {
            
                ParticleFile <<  std::fixed << std::setprecision(8) <<
                     ThisParticleList[iP].c1
                        <<"    "<<  ThisParticleList[iP].c2
                        <<"    "<<  ThisParticleList[iP].c3
                        <<"    "<<  ThisParticleList[iP].v_c1
                        <<"    "<<  ThisParticleList[iP].v_c2
                        <<"    "<<  ThisParticleList[iP].v_c3
                        <<"    "<<  sqrt(ThisParticleList[iP].v_c1*ThisParticleList[iP].v_c1 + ThisParticleList[iP].v_c3*ThisParticleList[iP].v_c3)
                        << endl;
                
                //cout << ThisParticleList <<
            }
        }
      //  ParticleFile.close();

		for(int iP=0;iP <nP;iP++) {
			vector<C3Vec> thisOrbitE1_re_XYZ(nSteps,C3Vec(0,0,0));
			vector<C3Vec> thisOrbitE1_im_XYZ(nSteps,C3Vec(0,0,0));

			vector<C3Vec> thisOrbitB1_re_XYZ(nSteps,C3Vec(0,0,0));
			vector<C3Vec> thisOrbitB1_im_XYZ(nSteps,C3Vec(0,0,0));

			CParticle thisParticle_XYZ(ThisParticleList[iP]);
            
			float qOverm =  thisParticle_XYZ.q/thisParticle_XYZ.m;
            
			float Ze = thisParticle_XYZ.q;

#if LOWMEM_ORBIT_WRITE >= 1
            ofstream OrbitFile;
            ofstream v1File;
			ofstream e1_dot_grad_File;
			ofstream df0dv_File;

            stringstream OrbitFileName;
            OrbitFileName << "output/unperturbed/orbit";
            OrbitFileName << '_' << setw(3) << setfill('0') <<  iList << '_' << setw(3) << setfill('0') << iP;
            OrbitFileName << ".txt";
            
            int write_iX;
            int write_iP;
            
            if(SaveSingleOrbit){
                 write_iX = 0;
                 write_iP = 5;
               
                if(iList==write_iX && iP==write_iP) {
                    cout<<"Write Particle Properties:"<<endl;
                    cout<<" vTh: "<<thisParticle_XYZ.vTh<<endl;
                    cout<<" v1: "<<thisParticle_XYZ.v_c1<<endl;
                    cout<<" v2: "<<thisParticle_XYZ.v_c2<<endl;
                    cout<<" v3: "<<thisParticle_XYZ.v_c3<<endl;

                    OrbitFile.open("output/orbit.txt", ios::out | ios::trunc);
                    OrbitFile<<"wc / wrf: "<< wrf_wc[iList]<<endl;
                    OrbitFile<<" t  x  y  z  re(e1)  im(e1)  re(e2)  im(e2)  re(e3)  im(e3)  re(b1)  im(b1)  re(b2)  im(b2)  re(b3)  im(b3)"<<endl;
                    v1File.open("output/orbit_v1.txt", ios::out | ios::trunc);
                    v1File<<" t  re(v11)  im(v11)  re(v12)  im(v12)  re(v13)  im(v13)"<<endl;
                    e1_dot_grad_File.open("output/orbit_e1_dot_grad_df0_dv.txt", ios::out | ios::trunc);
                    e1_dot_grad_File<<" t  re(v1xb01)  im(v1xb01)  re(v1xb02)  im(v1xb02)  re(v1xb03)  im(v1xb03)"<<endl;
                    df0dv_File.open("output/df0dv.txt", ios::out | ios::trunc);
                    df0dv_File<<" t  vx  vy  vz  valp  vbet  vpar  vper  gyroAngle  df0dv_x  df0dv_y  df0dv_z"<<endl;
                }
             }
            
             else{
                OrbitFile.open(OrbitFileName.str().c_str(), ios::out | ios::trunc);
                OrbitFile<<"wc / wrf: "<< wrf_wc[iList]<<endl;
                OrbitFile<<" t  x  y  z vpar vper vx vy vz re(e1)  im(e1)  re(e2)  im(e2)  re(e3)  im(e3)  re(b1)  im(b1)  re(b2)  im(b2)  re(b3)  im(b3)"<<endl;
            }
#endif    
			// generate  orbit and get time-harmonic e along it

			vector<C3Vec> thisOrbit_XYZ(nSteps);

			vector<C3VecI> thisE1c_XYZ(nSteps,C3VecI());
			vector<C3VecI> thisB1c_XYZ(nSteps,C3VecI());
			C3VecI thisV1c_(0,0,0), dVc(0,0,0), crossTerm(0,0,0);
            vector<complex<float> > this_e1_dot_gradvf0(nSteps);

	 		for(int i=0;i<nSteps;i++) {	

				thisOrbit_XYZ[i] = C3Vec(thisParticle_XYZ.c1,thisParticle_XYZ.c2,thisParticle_XYZ.c3);
#if GC_ORBITS >=1 
                int MoveStatus = rk4_move_gc ( thisParticle_XYZ, dtMin, thisT[i], 
                                r, b0_CYL, r_gc, curv_CYL, grad_CYL, bDotGradB, wrf );
#else
				 int MoveStatus = rk4_move ( thisParticle_XYZ, dtMin, thisT[i], b0_CYL, fieldMesh );

               // int MoveStatus = rk4_move ( thisParticle_XYZ, dtMin, thisT[i], b0_CYL,r, e1Re_XYZ, e1Im_XYZ, wrf );
                
#endif
                int OverallStatus = max(thisParticle_XYZ.status,MoveStatus);
#if DEBUG_MOVE >=1 
                cout << "Position After Move: " << thisParticle_XYZ.c1 << "  " << thisParticle_XYZ.c2 << "  " << thisParticle_XYZ.c3 << endl;
                if(MoveStatus>0) {
                    cout<<"ERROR: rk4_move_gc threw an error"<<endl;
                    cout<<"MoveStatus: "<<MoveStatus<<endl;
                    exit(1); 
                }
#endif
            
                C3Vec thisPos(thisParticle_XYZ.c1,thisParticle_XYZ.c2,thisParticle_XYZ.c3);
                C3Vec thisVel_XYZ(thisParticle_XYZ.v_c1,thisParticle_XYZ.v_c2,thisParticle_XYZ.v_c3);
				C3Vec thisB0 = kj_interp (XYZ_to_CYL(C3Vec(thisOrbit_XYZ[i].c1,thisOrbit_XYZ[i].c2,thisOrbit_XYZ[i].c3)), fieldMesh, b0_CYL, istat );

                float this_Theta = sqrt(pow(thisParticle_XYZ.c1,2)+pow(thisParticle_XYZ.c2,2));

				C3Vec gradv_f0_XYZ = maxwellian_df0_dv ( thisVel_XYZ, T_keV[iList], density_m3[iList], thisParticle_XYZ.amu, thisParticle_XYZ.Z );

				C3Vec e1ReTmp_XYZ = kj_interp ( XYZ_to_CYL(C3Vec(thisOrbit_XYZ[i].c1,thisOrbit_XYZ[i].c2,thisOrbit_XYZ[i].c3)), fieldMesh, e1Re_XYZ, istat );
				C3Vec e1ImTmp_XYZ = kj_interp ( XYZ_to_CYL(C3Vec(thisOrbit_XYZ[i].c1,thisOrbit_XYZ[i].c2,thisOrbit_XYZ[i].c3)), fieldMesh, e1Im_XYZ, istat );
				C3Vec b1ReTmp_XYZ = kj_interp ( XYZ_to_CYL(C3Vec(thisOrbit_XYZ[i].c1,thisOrbit_XYZ[i].c2,thisOrbit_XYZ[i].c3)), fieldMesh, b1Re_XYZ, istat );
				C3Vec b1ImTmp_XYZ = kj_interp ( XYZ_to_CYL(C3Vec(thisOrbit_XYZ[i].c1,thisOrbit_XYZ[i].c2,thisOrbit_XYZ[i].c3)), fieldMesh, b1Im_XYZ, istat );
	
				thisOrbitE1_re_XYZ[i] = e1ReTmp_XYZ*(1-OverallStatus);
				thisOrbitE1_im_XYZ[i] = e1ImTmp_XYZ*(1-OverallStatus);
				thisOrbitB1_re_XYZ[i] = b1ReTmp_XYZ*(1-OverallStatus);
				thisOrbitB1_im_XYZ[i] = b1ImTmp_XYZ*(1-OverallStatus);

				float tTmp = thisT[i];
				float weight = hanningWeight[i];
                float phs = -(wrf * tTmp); 
                complex<float> _i(0,1);
                complex<float> exp_inphi = exp(_i*(float)nPhi*this_Theta);

                weight = 1.0;
                phs = 0.0;
                exp_inphi = 1.0;

				thisE1c_XYZ[i] = exp_inphi * C3VecI(
								weight*complex<float>(
										thisOrbitE1_re_XYZ[i].c1*cos(phs)-thisOrbitE1_im_XYZ[i].c1*sin(phs),
										thisOrbitE1_im_XYZ[i].c1*cos(phs)+thisOrbitE1_re_XYZ[i].c1*sin(phs)),
								weight*complex<float>(
										thisOrbitE1_re_XYZ[i].c2*cos(phs)-thisOrbitE1_im_XYZ[i].c2*sin(phs),
										thisOrbitE1_im_XYZ[i].c2*cos(phs)+thisOrbitE1_re_XYZ[i].c2*sin(phs)),
								weight*complex<float>(
										thisOrbitE1_re_XYZ[i].c3*cos(phs)-thisOrbitE1_im_XYZ[i].c3*sin(phs),
										thisOrbitE1_im_XYZ[i].c3*cos(phs)+thisOrbitE1_re_XYZ[i].c3*sin(phs))
								);	

				thisB1c_XYZ[i] = exp_inphi * C3VecI(
								weight*complex<float>(
										thisOrbitB1_re_XYZ[i].c1*cos(phs)-thisOrbitB1_im_XYZ[i].c1*sin(phs),
										thisOrbitB1_im_XYZ[i].c1*cos(phs)+thisOrbitB1_re_XYZ[i].c1*sin(phs)),
								weight*complex<float>(
										thisOrbitB1_re_XYZ[i].c2*cos(phs)-thisOrbitB1_im_XYZ[i].c2*sin(phs),
										thisOrbitB1_im_XYZ[i].c2*cos(phs)+thisOrbitB1_re_XYZ[i].c2*sin(phs)),
								weight*complex<float>(
										thisOrbitB1_re_XYZ[i].c3*cos(phs)-thisOrbitB1_im_XYZ[i].c3*sin(phs),
										thisOrbitB1_im_XYZ[i].c3*cos(phs)+thisOrbitB1_re_XYZ[i].c3*sin(phs))
								);	

                //complex<float> _i(0,1);
				//complex<float> phs = -(_i*wrf*tTmp);
                //complex<float> amp1(thisOrbitE_re_XYZ[i].c1,thisOrbitE_im_XYZ[i].c1);
                //complex<float> amp2(thisOrbitE_re_XYZ[i].c2,thisOrbitE_im_XYZ[i].c2);
                //complex<float> amp3(thisOrbitE_re_XYZ[i].c3,thisOrbitE_im_XYZ[i].c3);
                //complex<float> expPart = exp(phs);
                //thisE1c[i] = weight*C3VecI(amp1*expPart,amp2*expPart,amp3*expPart); 

#if DEBUG_MOVE >= 2
                cout << "thisE1c[i].c1: "<<thisE1c_XYZ[i].c1<<endl;
                cout << "thisE1c[i].c2: "<<thisE1c_XYZ[i].c2<<endl;
                cout << "thisE1c[i].c3: "<<thisE1c_XYZ[i].c3<<endl;

                cout << "thisB1c[i].c1: "<<thisB1c_XYZ[i].c1<<endl;
                cout << "thisB1c[i].c2: "<<thisB1c_XYZ[i].c2<<endl;
                cout << "thisB1c[i].c3: "<<thisB1c_XYZ[i].c3<<endl;
#endif
#if DEBUG_FORCE_TERM >= 1
                cout << "thisE1c[i].c1: "<<thisE1c_XYZ[i].c1<<endl;
                cout << "thisE1c[i].c2: "<<thisE1c_XYZ[i].c2<<endl;
                cout << "thisE1c[i].c3: "<<thisE1c_XYZ[i].c3<<endl;

                cout << "thisB1c[i].c1: "<<thisB1c_XYZ[i].c1<<endl;
                cout << "thisB1c[i].c2: "<<thisB1c_XYZ[i].c2<<endl;
                cout << "thisB1c[i].c3: "<<thisB1c_XYZ[i].c3<<endl;

                cout << "thisVel_XYZ.c1: "<<thisVel_XYZ.c1<<endl;
                cout << "thisVel_XYZ.c2: "<<thisVel_XYZ.c2<<endl;
                cout << "thisVel_XYZ.c3: "<<thisVel_XYZ.c3<<endl;

#endif
                this_e1_dot_gradvf0[i] = dot(thisE1c_XYZ[i], gradv_f0_XYZ);

				//this_e1_dot_df0dv[i].c1 = thisE1c[i].c1 * thisdf0_dv.c1;  
				//this_e1_dot_df0dv[i].c2 = thisE1c[i].c2 * thisdf0_dv.c2;  
				//this_e1_dot_df0dv[i].c3 = thisE1c[i].c3 * thisdf0_dv.c3;  

#if LOWMEM_ORBIT_WRITE >= 1

                if(SaveSingleOrbit){
                
                    if(iList==write_iX && iP==write_iP) {
                        df0dv_File<<scientific;
                        df0dv_File<< thisT[i]
                                <<"    "<< thisVel_XYZ.c1
                                <<"    "<< thisVel_XYZ.c2
                                <<"    "<< thisVel_XYZ.c3 
                                <<"    "<< thisParticle_XYZ.vAlp
                                <<"    "<< thisParticle_XYZ.vBet
                                <<"    "<< thisParticle_XYZ.vPar
                                <<"    "<< thisParticle_XYZ.vPer
                                <<"    "<< thisParticle_XYZ.phs
                                <<"    "<< gradv_f0_XYZ.c1 
                                <<"    "<< gradv_f0_XYZ.c2 
                                <<"    "<< gradv_f0_XYZ.c3 
                                << endl;
                    }
     
                    if(iList==write_iX && iP==write_iP) {
                        OrbitFile<<scientific;
                        OrbitFile<< thisT[i]
                                <<"    "<< thisPos.c1
                                <<"    "<< thisPos.c2
                                <<"    "<< thisPos.c3 
                                <<"    "<< real(thisE1c_XYZ[i].c1)
                                <<"    "<< imag(thisE1c_XYZ[i].c1)
                                <<"    "<< real(thisE1c_XYZ[i].c2)
                                <<"    "<< imag(thisE1c_XYZ[i].c2)
                                <<"    "<< real(thisE1c_XYZ[i].c3)
                                <<"    "<< imag(thisE1c_XYZ[i].c3)
                                <<"    "<< real(thisB1c_XYZ[i].c1)
                                <<"    "<< imag(thisB1c_XYZ[i].c1)
                                <<"    "<< real(thisB1c_XYZ[i].c2)
                                <<"    "<< imag(thisB1c_XYZ[i].c2)
                                <<"    "<< real(thisB1c_XYZ[i].c3)
                                <<"    "<< imag(thisB1c_XYZ[i].c3)
                                << endl;
                    }
                    if(iList==write_iX && iP==write_iP) {
                        e1_dot_grad_File<<scientific;
                        e1_dot_grad_File<< thisT[i]
                                <<"    "<< real(this_e1_dot_gradvf0[i])
                                <<"    "<< imag(this_e1_dot_gradvf0[i])
                                << endl;
                    }
                    
                }
                else{
                        if ( i % (int)floor(nStepsPerCycle/nSavePerRFCycle) == 0){
    
                            int tmp_Stat;
                            C3Vec b0_XYZ_T_at_ThisPos = kj_interp(XYZ_to_CYL(C3Vec(thisPos.c1,thisPos.c2,thisPos.c3)),fieldMesh,b0_XYZ,tmp_Stat);
                            
                            C3Vec thisV_abp = rot_XYZ_to_abp (thisVel_XYZ,b0_XYZ_T_at_ThisPos, 0 );
                            float vPar = thisV_abp.c3;
                            float vPer = sqrt(pow(thisV_abp.c1,2)+pow(thisV_abp.c2,2));

                            //OrbitFile<<scientific;
                            OrbitFile <<  std::fixed << std::setprecision(12) << thisT[i]
                                <<"    "<< thisPos.c1
                                <<"    "<< thisPos.c2
                                <<"    "<< thisPos.c3
                                <<"    "<< vPar
                                <<"    "<< vPer
                                <<"    "<< thisVel_XYZ.c1
                                <<"    "<< thisVel_XYZ.c2
                                <<"    "<< thisVel_XYZ.c3
                                <<"    "<< real(thisE1c_XYZ[i].c1)
                                <<"    "<< imag(thisE1c_XYZ[i].c1)
                                <<"    "<< real(thisE1c_XYZ[i].c2)
                                <<"    "<< imag(thisE1c_XYZ[i].c2)
                                <<"    "<< real(thisE1c_XYZ[i].c3)
                                <<"    "<< imag(thisE1c_XYZ[i].c3)
                                <<"    "<< real(thisB1c_XYZ[i].c1)
                                <<"    "<< imag(thisB1c_XYZ[i].c1)
                                <<"    "<< real(thisB1c_XYZ[i].c2)
                                <<"    "<< imag(thisB1c_XYZ[i].c2)
                                <<"    "<< real(thisB1c_XYZ[i].c3)
                                <<"    "<< imag(thisB1c_XYZ[i].c3)
                                << endl;
                                }
                }
                
#endif
			}
#if LOWMEM_ORBIT_WRITE >= 1
            if(iList==write_iX && iP==write_iP) {
                OrbitFile.close();
            }
#endif
			complex<float> this_f1c = -qOverm * intVecArray ( thisT, this_e1_dot_gradvf0 );

#if LOWMEM_ORBIT_WRITE >= 1
            if(iList==write_iX && iP==write_iP) {

                complex<float> tmp = 0.0;
                for(int i=0;i<nSteps;i++) {
						//float tmp_vTh = get_vTh ( thisParticle_XYZ.amu, thisParticle_XYZ.Z, T_keV[iX] );
                        tmp += -qOverm * this_e1_dot_gradvf0[i] * dtMin;
                        v1File<<thisT[i]
                                <<"    "<< real(tmp)
                                <<"    "<< imag(tmp)
                                << endl;
                }
				//cout<<"tmp: "<<tmp.c1<<"  "<<tmp.c2<<"  "<<tmp.c3<<endl;
				//cout<<"thisV1c: "<<thisV1c.c1<<"  "<<thisV1c.c2<<"  "<<thisV1c.c3<<endl;
            }
#endif
			f1c[iP] = -this_f1c;
			//f1xc[iP] = -thisV1c.c1;
			//f1yc[iP] = -thisV1c.c2;
			//f1zc[iP] = -thisV1c.c3;

			float v0x_i = ThisParticleList[iP].v_c1;
			float v0y_i = ThisParticleList[iP].v_c2;
			float v0z_i = ThisParticleList[iP].v_c3;

			float h = dv * Ze;

			#pragma omp critical // "atomic" does not work for complex numbers
			{
				//j1xc[iX] += h * ( v0x_i*f1xc[iP] ); 
				//j1yc[iX] += h * ( v0y_i*f1yc[iP] ); 
				//j1zc[iX] += h * ( v0z_i*f1zc[iP] ); 


////////// DEBUGGING!
				//j1xc[iList] += h * ( v0x_i*f1c[iP] );
				//j1yc[iList] += h * ( v0y_i*f1c[iP] );
				//j1zc[iList] += h * ( v0z_i*f1c[iP] );
			}
		}

///////////////////
///////// Begin perturbed orbit optional section
///////////////////
#if COMPUTE_PERTURBED_ORBITS >= 1

		for(int iP=0;iP <nP;iP++) {

			vector<C3Vec> thisOrbitE1_re_XYZ(nSteps,C3Vec(0,0,0));
			vector<C3Vec> thisOrbitE1_im_XYZ(nSteps,C3Vec(0,0,0));

			vector<C3Vec> thisOrbitB1_re_XYZ(nSteps,C3Vec(0,0,0));
			vector<C3Vec> thisOrbitB1_im_XYZ(nSteps,C3Vec(0,0,0));
            CParticle thisParticle_XYZ(ThisParticleList[iP]);

			float qOverm =  thisParticle_XYZ.q/thisParticle_XYZ.m;
            
			float Ze = thisParticle_XYZ.q;

#if LOWMEM_ORBIT_WRITE >= 1
            ofstream OrbitFile;
            ofstream v1File;
			ofstream e1_dot_grad_File;
			ofstream df0dv_File;

            stringstream OrbitFileName;
            OrbitFileName << "output/perturbed/orbit";
            OrbitFileName << '_' << setw(3) << setfill('0') <<  iList << '_' << setw(3) << setfill('0') << iP;
            OrbitFileName << ".txt";
            
            int write_iX;
            int write_iP;
            if(SaveSingleOrbit){
                 write_iX = 0;
                 write_iP = 5;
               
                if(iList==write_iX && iP==write_iP) {
                    cout<<"Write Particle Properties:"<<endl;
                    cout<<" vTh: "<<thisParticle_XYZ.vTh<<endl;
                    cout<<" v1: "<<thisParticle_XYZ.v_c1<<endl;
                    cout<<" v2: "<<thisParticle_XYZ.v_c2<<endl;
                    cout<<" v3: "<<thisParticle_XYZ.v_c3<<endl;

                    OrbitFile.open("output/orbit.txt", ios::out | ios::trunc);
                    OrbitFile<<"wc / wrf: "<< wrf_wc[iList]<<endl;
                    OrbitFile<<" t  x  y  z  re(e1)  im(e1)  re(e2)  im(e2)  re(e3)  im(e3)  re(b1)  im(b1)  re(b2)  im(b2)  re(b3)  im(b3)"<<endl;
                    v1File.open("output/orbit_v1.txt", ios::out | ios::trunc);
                    v1File<<" t  re(v11)  im(v11)  re(v12)  im(v12)  re(v13)  im(v13)"<<endl;
                    e1_dot_grad_File.open("output/orbit_e1_dot_grad_df0_dv.txt", ios::out | ios::trunc);
                    e1_dot_grad_File<<" t  re(v1xb01)  im(v1xb01)  re(v1xb02)  im(v1xb02)  re(v1xb03)  im(v1xb03)"<<endl;
                    df0dv_File.open("output/df0dv.txt", ios::out | ios::trunc);
                    df0dv_File<<" t  vx  vy  vz  valp  vbet  vpar  vper  gyroAngle  df0dv_x  df0dv_y  df0dv_z"<<endl;
                }
             }
            
             else{
                OrbitFile.open(OrbitFileName.str().c_str(), ios::out | ios::trunc);
                OrbitFile<<"wc / wrf: "<< wrf_wc[iList]<<endl;
                OrbitFile<<" t  x  y  z vpar vper vx vy vz re(e1)  im(e1)  re(e2)  im(e2)  re(e3)  im(e3)  re(b1)  im(b1)  re(b2)  im(b2)  re(b3)  im(b3)"<<endl;
                
            }
#endif    
			// generate  orbit and get time-harmonic e along it

			vector<C3Vec> thisOrbit_XYZ(nSteps);

			vector<C3VecI> thisE1c_XYZ(nSteps,C3VecI());
			vector<C3VecI> thisB1c_XYZ(nSteps,C3VecI());
			C3VecI thisV1c_(0,0,0), dVc(0,0,0), crossTerm(0,0,0);
            vector<complex<float> > this_e1_dot_gradvf0(nSteps);

	 		for(int i=0;i<nSteps;i++) {	
				thisOrbit_XYZ[i] = C3Vec(thisParticle_XYZ.c1,thisParticle_XYZ.c2,thisParticle_XYZ.c3);
#if GC_ORBITS >=1 
                int MoveStatus = rk4_move_gc ( thisParticle_XYZ, dtMin, thisT[i], 
                                r, b0_CYL, r_gc, curv_CYL, grad_CYL, bDotGradB, wrf );
#else
				//int MoveStatus = rk4_move ( thisParticle_XYZ, dtMin, thisT[i], b0_CYL, r );

                int MoveStatus = rk4_move ( thisParticle_XYZ, dtMin, thisT[i], b0_CYL,fieldMesh, e1Re_XYZ, e1Im_XYZ, wrf );
                
#endif
                int OverallStatus = max(thisParticle_XYZ.status,MoveStatus);
#if DEBUG_MOVE >=1
                cout << "Position After Move: " << thisParticle_XYZ.c1 << "  " << thisParticle_XYZ.c2 << "  " << thisParticle_XYZ.c3 << endl;
                if(MoveStatus>0) {
                    cout<<"ERROR: rk4_move_gc threw an error"<<endl;
                    cout<<"MoveStatus: "<<MoveStatus<<endl;
                    exit(1); 
                }
#endif
            
                C3Vec thisPos(thisParticle_XYZ.c1,thisParticle_XYZ.c2,thisParticle_XYZ.c3);
                C3Vec thisVel_XYZ(thisParticle_XYZ.v_c1,thisParticle_XYZ.v_c2,thisParticle_XYZ.v_c3);
				C3Vec thisB0 = kj_interp ( XYZ_to_CYL(C3Vec(thisOrbit_XYZ[i].c1,thisOrbit_XYZ[i].c2,thisOrbit_XYZ[i].c3)), fieldMesh, b0_CYL, istat );

                float this_Theta = sqrt(pow(thisParticle_XYZ.c1,2)+pow(thisParticle_XYZ.c2,2));

				C3Vec gradv_f0_XYZ = maxwellian_df0_dv ( thisVel_XYZ, T_keV[iList], density_m3[iList], thisParticle_XYZ.amu, thisParticle_XYZ.Z );

				C3Vec e1ReTmp_XYZ = kj_interp ( XYZ_to_CYL(C3Vec(thisOrbit_XYZ[i].c1,thisOrbit_XYZ[i].c2,thisOrbit_XYZ[i].c3)), fieldMesh, e1Re_XYZ, istat );
				C3Vec e1ImTmp_XYZ = kj_interp ( XYZ_to_CYL(C3Vec(thisOrbit_XYZ[i].c1,thisOrbit_XYZ[i].c2,thisOrbit_XYZ[i].c3)), fieldMesh, e1Im_XYZ, istat );
				C3Vec b1ReTmp_XYZ = kj_interp ( XYZ_to_CYL(C3Vec(thisOrbit_XYZ[i].c1,thisOrbit_XYZ[i].c2,thisOrbit_XYZ[i].c3)), fieldMesh, b1Re_XYZ, istat );
				C3Vec b1ImTmp_XYZ = kj_interp ( XYZ_to_CYL(C3Vec(thisOrbit_XYZ[i].c1,thisOrbit_XYZ[i].c2,thisOrbit_XYZ[i].c3)), fieldMesh, b1Im_XYZ, istat );
	
				thisOrbitE1_re_XYZ[i] = e1ReTmp_XYZ*(1-OverallStatus);
				thisOrbitE1_im_XYZ[i] = e1ImTmp_XYZ*(1-OverallStatus);
				thisOrbitB1_re_XYZ[i] = b1ReTmp_XYZ*(1-OverallStatus);
				thisOrbitB1_im_XYZ[i] = b1ImTmp_XYZ*(1-OverallStatus);

				float tTmp = thisT[i];
				float weight = hanningWeight[i];
                float phs = -(wrf * tTmp);
                complex<float> _i(0,1);
                complex<float> exp_inphi = exp(_i*(float)nPhi*this_Theta);
 
                weight = 1.0;
                exp_inphi = 1.0;
                
				thisE1c_XYZ[i] = exp_inphi * C3VecI(
								weight*complex<float>(
										thisOrbitE1_re_XYZ[i].c1*cos(phs)-thisOrbitE1_im_XYZ[i].c1*sin(phs),
										thisOrbitE1_im_XYZ[i].c1*cos(phs)+thisOrbitE1_re_XYZ[i].c1*sin(phs)),
								weight*complex<float>(
										thisOrbitE1_re_XYZ[i].c2*cos(phs)-thisOrbitE1_im_XYZ[i].c2*sin(phs),
										thisOrbitE1_im_XYZ[i].c2*cos(phs)+thisOrbitE1_re_XYZ[i].c2*sin(phs)),
								weight*complex<float>(
										thisOrbitE1_re_XYZ[i].c3*cos(phs)-thisOrbitE1_im_XYZ[i].c3*sin(phs),
										thisOrbitE1_im_XYZ[i].c3*cos(phs)+thisOrbitE1_re_XYZ[i].c3*sin(phs))
								);	

				thisB1c_XYZ[i] = exp_inphi * C3VecI(
								weight*complex<float>(
										thisOrbitB1_re_XYZ[i].c1*cos(phs)-thisOrbitB1_im_XYZ[i].c1*sin(phs),
										thisOrbitB1_im_XYZ[i].c1*cos(phs)+thisOrbitB1_re_XYZ[i].c1*sin(phs)),
								weight*complex<float>(
										thisOrbitB1_re_XYZ[i].c2*cos(phs)-thisOrbitB1_im_XYZ[i].c2*sin(phs),
										thisOrbitB1_im_XYZ[i].c2*cos(phs)+thisOrbitB1_re_XYZ[i].c2*sin(phs)),
								weight*complex<float>(
										thisOrbitB1_re_XYZ[i].c3*cos(phs)-thisOrbitB1_im_XYZ[i].c3*sin(phs),
										thisOrbitB1_im_XYZ[i].c3*cos(phs)+thisOrbitB1_re_XYZ[i].c3*sin(phs))
								);	

                //complex<float> _i(0,1);
				//complex<float> phs = -(_i*wrf*tTmp);
                //complex<float> amp1(thisOrbitE_re_XYZ[i].c1,thisOrbitE_im_XYZ[i].c1);
                //complex<float> amp2(thisOrbitE_re_XYZ[i].c2,thisOrbitE_im_XYZ[i].c2);
                //complex<float> amp3(thisOrbitE_re_XYZ[i].c3,thisOrbitE_im_XYZ[i].c3);
                //complex<float> expPart = exp(phs);
                //thisE1c[i] = weight*C3VecI(amp1*expPart,amp2*expPart,amp3*expPart); 

#if DEBUG_MOVE >= 2
                cout << "thisE1c[i].c1: "<<thisE1c_XYZ[i].c1<<endl;
                cout << "thisE1c[i].c2: "<<thisE1c_XYZ[i].c2<<endl;
                cout << "thisE1c[i].c3: "<<thisE1c_XYZ[i].c3<<endl;

                cout << "thisB1c[i].c1: "<<thisB1c_XYZ[i].c1<<endl;
                cout << "thisB1c[i].c2: "<<thisB1c_XYZ[i].c2<<endl;
                cout << "thisB1c[i].c3: "<<thisB1c_XYZ[i].c3<<endl;
#endif
#if DEBUG_FORCE_TERM >= 1
                cout << "thisE1c[i].c1: "<<thisE1c_XYZ[i].c1<<endl;
                cout << "thisE1c[i].c2: "<<thisE1c_XYZ[i].c2<<endl;
                cout << "thisE1c[i].c3: "<<thisE1c_XYZ[i].c3<<endl;

                cout << "thisB1c[i].c1: "<<thisB1c_XYZ[i].c1<<endl;
                cout << "thisB1c[i].c2: "<<thisB1c_XYZ[i].c2<<endl;
                cout << "thisB1c[i].c3: "<<thisB1c_XYZ[i].c3<<endl;

                cout << "thisVel_XYZ.c1: "<<thisVel_XYZ.c1<<endl;
                cout << "thisVel_XYZ.c2: "<<thisVel_XYZ.c2<<endl;
                cout << "thisVel_XYZ.c3: "<<thisVel_XYZ.c3<<endl;

#endif
                this_e1_dot_gradvf0[i] = dot(thisE1c_XYZ[i], gradv_f0_XYZ);

				//this_e1_dot_df0dv[i].c1 = thisE1c[i].c1 * thisdf0_dv.c1;  
				//this_e1_dot_df0dv[i].c2 = thisE1c[i].c2 * thisdf0_dv.c2;  
				//this_e1_dot_df0dv[i].c3 = thisE1c[i].c3 * thisdf0_dv.c3;  

#if LOWMEM_ORBIT_WRITE >= 1

                if(SaveSingleOrbit){
                
                    if(iList==write_iX && iP==write_iP) {
                        df0dv_File<<scientific;
                        df0dv_File<< thisT[i]
                                <<"    "<< thisVel_XYZ.c1
                                <<"    "<< thisVel_XYZ.c2
                                <<"    "<< thisVel_XYZ.c3 
                                <<"    "<< thisParticle_XYZ.vAlp
                                <<"    "<< thisParticle_XYZ.vBet
                                <<"    "<< thisParticle_XYZ.vPar
                                <<"    "<< thisParticle_XYZ.vPer
                                <<"    "<< thisParticle_XYZ.phs
                                <<"    "<< gradv_f0_XYZ.c1 
                                <<"    "<< gradv_f0_XYZ.c2 
                                <<"    "<< gradv_f0_XYZ.c3 
                                << endl;
                    }
     
                    if(iList==write_iX && iP==write_iP) {
                        OrbitFile<<scientific;
                        OrbitFile<< thisT[i]
                                <<"    "<< thisPos.c1
                                <<"    "<< thisPos.c2
                                <<"    "<< thisPos.c3 
                                <<"    "<< real(thisE1c_XYZ[i].c1)
                                <<"    "<< imag(thisE1c_XYZ[i].c1)
                                <<"    "<< real(thisE1c_XYZ[i].c2)
                                <<"    "<< imag(thisE1c_XYZ[i].c2)
                                <<"    "<< real(thisE1c_XYZ[i].c3)
                                <<"    "<< imag(thisE1c_XYZ[i].c3)
                                <<"    "<< real(thisB1c_XYZ[i].c1)
                                <<"    "<< imag(thisB1c_XYZ[i].c1)
                                <<"    "<< real(thisB1c_XYZ[i].c2)
                                <<"    "<< imag(thisB1c_XYZ[i].c2)
                                <<"    "<< real(thisB1c_XYZ[i].c3)
                                <<"    "<< imag(thisB1c_XYZ[i].c3)
                                << endl;
                    }
                    if(iList==write_iX && iP==write_iP) {
                        e1_dot_grad_File<<scientific;
                        e1_dot_grad_File<< thisT[i]
                                <<"    "<< real(this_e1_dot_gradvf0[i])
                                <<"    "<< imag(this_e1_dot_gradvf0[i])
                                << endl;
                    }
                    
                }
                else{
                        if ( i % (int)floor(nStepsPerCycle/nSavePerRFCycle) == 0){
    
                            int tmp_Stat;
                            C3Vec b0_XYZ_T_at_ThisPos = kj_interp(XYZ_to_CYL(C3Vec(thisPos.c1,thisPos.c2,thisPos.c3)),fieldMesh,b0_XYZ,tmp_Stat);
                            C3Vec thisV_abp = rot_XYZ_to_abp (thisVel_XYZ,b0_XYZ_T_at_ThisPos, 0 );
                            float vPar = thisV_abp.c3;
                            float vPer = sqrt(pow(thisV_abp.c1,2)+pow(thisV_abp.c2,2));

                            //OrbitFile<<scientific;
                            OrbitFile <<  std::fixed << std::setprecision(12) << thisT[i]
                                <<"    "<< thisPos.c1
                                <<"    "<< thisPos.c2
                                <<"    "<< thisPos.c3
                                <<"    "<< vPar
                                <<"    "<< vPer
                                <<"    "<< thisVel_XYZ.c1
                                <<"    "<< thisVel_XYZ.c2
                                <<"    "<< thisVel_XYZ.c3
                                <<"    "<< real(thisE1c_XYZ[i].c1)
                                <<"    "<< imag(thisE1c_XYZ[i].c1)
                                <<"    "<< real(thisE1c_XYZ[i].c2)
                                <<"    "<< imag(thisE1c_XYZ[i].c2)
                                <<"    "<< real(thisE1c_XYZ[i].c3)
                                <<"    "<< imag(thisE1c_XYZ[i].c3)
                                <<"    "<< real(thisB1c_XYZ[i].c1)
                                <<"    "<< imag(thisB1c_XYZ[i].c1)
                                <<"    "<< real(thisB1c_XYZ[i].c2)
                                <<"    "<< imag(thisB1c_XYZ[i].c2)
                                <<"    "<< real(thisB1c_XYZ[i].c3)
                                <<"    "<< imag(thisB1c_XYZ[i].c3)
                                << endl;
                                }
                }
                
#endif
			}
#if LOWMEM_ORBIT_WRITE >= 1
            if(iList==write_iX && iP==write_iP) {
                OrbitFile.close();
            }
#endif
			complex<float> this_f1c = -qOverm * intVecArray ( thisT, this_e1_dot_gradvf0 );

#if LOWMEM_ORBIT_WRITE >= 1
            if(iList==write_iX && iP==write_iP) {

                complex<float> tmp = 0.0;
                for(int i=0;i<nSteps;i++) {
						//float tmp_vTh = get_vTh ( thisParticle_XYZ.amu, thisParticle_XYZ.Z, T_keV[iX] );
                        tmp += -qOverm * this_e1_dot_gradvf0[i] * dtMin;
                        v1File<<thisT[i]
                                <<"    "<< real(tmp)
                                <<"    "<< imag(tmp)
                                << endl;
                }
				//cout<<"tmp: "<<tmp.c1<<"  "<<tmp.c2<<"  "<<tmp.c3<<endl;
				//cout<<"thisV1c: "<<thisV1c.c1<<"  "<<thisV1c.c2<<"  "<<thisV1c.c3<<endl;
            }
#endif
			f1c[iP] = -this_f1c;
			//f1xc[iP] = -thisV1c.c1;
			//f1yc[iP] = -thisV1c.c2;
			//f1zc[iP] = -thisV1c.c3;

			float v0x_i = ThisParticleList[iP].v_c1;
			float v0y_i = ThisParticleList[iP].v_c2;
			float v0z_i = ThisParticleList[iP].v_c3;

			float h = dv * Ze;

			#pragma omp critical // "atomic" does not work for complex numbers
			{
				//j1xc[iX] += h * ( v0x_i*f1xc[iP] ); 
				//j1yc[iX] += h * ( v0y_i*f1yc[iP] ); 
				//j1zc[iX] += h * ( v0z_i*f1zc[iP] ); 

				j1xc[iList] += h * ( v0x_i*f1c[iP] );
				j1yc[iList] += h * ( v0y_i*f1c[iP] );
				j1zc[iList] += h * ( v0z_i*f1c[iP] );
			}
		}

#endif
///////////////////
///////// End perturbed orbit optional section
///////////////////








#if CLOCK >= 1
#if not defined(_OPENMP)
        clock_t endTime = clock();
        double timeInSeconds = (endTime-startTime)/(double) CLOCKS_PER_SEC;
        cout<<"Time for this spatial point: "<<timeInSeconds<<endl;
        cout<<"Time per particle: "<<timeInSeconds/nP<<endl;  
#endif
#endif

#if LOWMEM_USEPAPI >= 1
		    cpuTime0=cpuTime;realTime0=realTime;flpIns0=flpIns;
		    papiReturn = PAPI_flops ( &realTime, &cpuTime, &flpIns, &mFlops );
		    printf("\nLOWMEM Oribit calculation performance ...\n");
		    printf("Real_time:\t%f\nProc_time:\t%f\nTotal flpins:\t%lld\nMFLOPS:\t\t%f\n",
					   realTime-realTime0, cpuTime-cpuTime0, flpIns-flpIns0, mFlops);
#endif

//#else // END OF LOWMEM CODING ^^^
//////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////

#endif


#if USEPAPI >= 1
		printf("\nGet e(t) and integrate performance ...\n");
		printf("Real_time:\t%f\nProc_time:\t%f\nTotal flpins:\t%lld\nMFLOPS:\t\t%f\n",
					   eT_realTime, eT_cpuTime, eT_flpIns, eT_mFlops/(nJp-1));
		printf("\nGet v(t) and integrate performance ...\n");
		printf("Real_time:\t%f\nProc_time:\t%f\nTotal flpins:\t%lld\nMFLOPS:\t\t%f\n",
					   vT_realTime, vT_cpuTime, vT_flpIns, vT_mFlops/(nJp-1));

		cpuTime0=cpuTime;realTime0=realTime;flpIns0=flpIns;
		papiReturn = PAPI_flops ( &realTime, &cpuTime, &flpIns, &mFlops );
		printf("\nj(t) performance ...\n");
		printf("Real_time:\t%f\nProc_time:\t%f\nTotal flpins:\t%lld\nMFLOPS:\t\t%f\n",
					   realTime-realTime0, cpuTime-cpuTime0, flpIns-flpIns0, mFlops);
#endif

#if __SAVE_ORBITS__>=1
		// Write orbits to file
	
		cout << "Writing orbits to file ... " << endl;

		stringstream ncOrbitsFileName;
		ncOrbitsFileName << "output/orbits_";
		ncOrbitsFileName << setw(3) << setfill('0') << iX;
	   	ncOrbitsFileName << ".nc"; 	

		try {
				// Really need to fix this but I don't know how to 
				// write a vector of structures using netCDF yet.

				NcFile ncOrbitsFile (ncOrbitsFileName.str().c_str(), NcFile::replace);
		
				NcDim nc_nP = ncOrbitsFile.addDim("nP", this_particles_XYZ.size());
				NcDim nc_nSteps = ncOrbitsFile.addDim("nSteps", nSteps);
				NcDim nc_nJp = ncOrbitsFile.addDim("nJp", nJp);
	
				vector<NcDim> nc_nPxnSteps(2);
				nc_nPxnSteps[0]=nc_nP;
				nc_nPxnSteps[1]=nc_nSteps;

				vector<NcDim> nc_nPxnJpxnSteps(3);
				nc_nPxnJpxnSteps[0]=nc_nP;
				nc_nPxnJpxnSteps[1]=nc_nJp;
				nc_nPxnJpxnSteps[2]=nc_nSteps;
		
				NcVar nc_t = ncOrbitsFile.addVar("t",ncFloat,nc_nSteps);
		
				NcVar nc_x = ncOrbitsFile.addVar("x",ncFloat,nc_nPxnSteps);
				NcVar nc_y = ncOrbitsFile.addVar("y",ncFloat,nc_nPxnSteps);
				NcVar nc_z = ncOrbitsFile.addVar("z",ncFloat,nc_nPxnSteps);

				NcVar nc_vx = ncOrbitsFile.addVar("vx",ncFloat,nc_nPxnSteps);
				NcVar nc_vy = ncOrbitsFile.addVar("vy",ncFloat,nc_nPxnSteps);
				NcVar nc_vz = ncOrbitsFile.addVar("vz",ncFloat,nc_nPxnSteps);
		
				NcVar nc_e1_x = ncOrbitsFile.addVar("e1_x",ncFloat,nc_nPxnSteps);
				NcVar nc_e1_y = ncOrbitsFile.addVar("e1_y",ncFloat,nc_nPxnSteps);
				NcVar nc_e1_z = ncOrbitsFile.addVar("e1_z",ncFloat,nc_nPxnSteps);

				NcVar nc_e1_x_re = ncOrbitsFile.addVar("e1_x_re",ncFloat,nc_nPxnSteps);
				NcVar nc_e1_y_re = ncOrbitsFile.addVar("e1_y_re",ncFloat,nc_nPxnSteps);
				NcVar nc_e1_z_re = ncOrbitsFile.addVar("e1_z_re",ncFloat,nc_nPxnSteps);

				NcVar nc_e1_x_im = ncOrbitsFile.addVar("e1_x_im",ncFloat,nc_nPxnSteps);
				NcVar nc_e1_y_im = ncOrbitsFile.addVar("e1_y_im",ncFloat,nc_nPxnSteps);
				NcVar nc_e1_z_im = ncOrbitsFile.addVar("e1_z_im",ncFloat,nc_nPxnSteps);

				NcVar nc_v1_x = ncOrbitsFile.addVar("v1x",ncFloat,nc_nPxnJpxnSteps);
				NcVar nc_v1_y = ncOrbitsFile.addVar("v1y",ncFloat,nc_nPxnJpxnSteps);
				NcVar nc_v1_z = ncOrbitsFile.addVar("v1z",ncFloat,nc_nPxnJpxnSteps);

				NcVar nc_v1_x_re = ncOrbitsFile.addVar("v1x_re",ncFloat,nc_nPxnJpxnSteps);
				NcVar nc_v1_y_re = ncOrbitsFile.addVar("v1y_re",ncFloat,nc_nPxnJpxnSteps);
				NcVar nc_v1_z_re = ncOrbitsFile.addVar("v1z_re",ncFloat,nc_nPxnJpxnSteps);

				NcVar nc_v1_x_im = ncOrbitsFile.addVar("v1x_im",ncFloat,nc_nPxnJpxnSteps);
				NcVar nc_v1_y_im = ncOrbitsFile.addVar("v1y_im",ncFloat,nc_nPxnJpxnSteps);
				NcVar nc_v1_z_im = ncOrbitsFile.addVar("v1z_im",ncFloat,nc_nPxnJpxnSteps);
	
				vector<size_t> startpA(2);
				vector<size_t> countpA(2);
				for(int iP=0;iP<this_particles_XYZ.size();iP++) {
		
						startpA[0]=iP;
						startpA[1]=0;
						countpA[0]=1;
						countpA[1]=nSteps;
		
						vector<float> tmpData (nSteps,0);
						for(int iS=0;iS<nSteps;iS++){tmpData[iS] = orbits_XYZ[iP][iS].c1;}
						nc_x.putVar(startpA,countpA,&tmpData[0]);
						for(int iS=0;iS<nSteps;iS++){tmpData[iS] = orbits_XYZ[iP][iS].c2;}
						nc_y.putVar(startpA,countpA,&tmpData[0]);
						for(int iS=0;iS<nSteps;iS++){tmpData[iS] = orbits_XYZ[iP][iS].c3;}
						nc_z.putVar(startpA,countpA,&tmpData[0]);

						for(int iS=0;iS<nSteps;iS++){tmpData[iS] = orbits_v_XYZ[iP][iS].c1;}
						nc_vx.putVar(startpA,countpA,&tmpData[0]);
						for(int iS=0;iS<nSteps;iS++){tmpData[iS] = orbits_v_XYZ[iP][iS].c2;}
						nc_vy.putVar(startpA,countpA,&tmpData[0]);
						for(int iS=0;iS<nSteps;iS++){tmpData[iS] = orbits_v_XYZ[iP][iS].c3;}
						nc_vz.putVar(startpA,countpA,&tmpData[0]);

						for(int iS=0;iS<nSteps;iS++){tmpData[iS] = e1[iP][iS].c1;}
						nc_e1_x.putVar(startpA,countpA,&tmpData[0]);
						for(int iS=0;iS<nSteps;iS++){tmpData[iS] = e1[iP][iS].c2;}
						nc_e1_y.putVar(startpA,countpA,&tmpData[0]);
						for(int iS=0;iS<nSteps;iS++){tmpData[iS] = e1[iP][iS].c3;}
						nc_e1_z.putVar(startpA,countpA,&tmpData[0]);

						for(int iS=0;iS<nSteps;iS++){tmpData[iS] = real(e1c[iP][iS].c1);}
						nc_e1_x_re.putVar(startpA,countpA,&tmpData[0]);
						for(int iS=0;iS<nSteps;iS++){tmpData[iS] = real(e1c[iP][iS].c2);}
						nc_e1_y_re.putVar(startpA,countpA,&tmpData[0]);
						for(int iS=0;iS<nSteps;iS++){tmpData[iS] = real(e1c[iP][iS].c3);}
						nc_e1_z_re.putVar(startpA,countpA,&tmpData[0]);

						for(int iS=0;iS<nSteps;iS++){tmpData[iS] = imag(e1c[iP][iS].c1);}
						nc_e1_x_im.putVar(startpA,countpA,&tmpData[0]);
						for(int iS=0;iS<nSteps;iS++){tmpData[iS] = imag(e1c[iP][iS].c2);}
						nc_e1_y_im.putVar(startpA,countpA,&tmpData[0]);
						for(int iS=0;iS<nSteps;iS++){tmpData[iS] = imag(e1c[iP][iS].c3);}
						nc_e1_z_im.putVar(startpA,countpA,&tmpData[0]);

				}

				vector<size_t> startpB(3);
				vector<size_t> countpB(3);
				for(int iP=0;iP<this_particles_XYZ.size();iP++) {
						for(int iJ=0;iJ<nJp;iJ++) {
		
							startpB[0]=iP;
							startpB[1]=iJ;
							startpB[2]=0;
							countpB[0]=1;
							countpB[1]=1;
							countpB[2]=nSteps;
		
							vector<float> tmpData (nSteps,0);

							for(int iS=0;iS<nSteps;iS++){tmpData[iS] = v1[iP][iJ][iS].c1;}
							nc_v1_x.putVar(startpB,countpB,&tmpData[0]);
							for(int iS=0;iS<nSteps;iS++){tmpData[iS] = v1[iP][iJ][iS].c2;}
							nc_v1_y.putVar(startpB,countpB,&tmpData[0]);
							for(int iS=0;iS<nSteps;iS++){tmpData[iS] = v1[iP][iJ][iS].c3;}
							nc_v1_z.putVar(startpB,countpB,&tmpData[0]);

							for(int iS=0;iS<nSteps;iS++){tmpData[iS] = real(v1c[iP][iJ][iS].c1);}
							nc_v1_x_re.putVar(startpB,countpB,&tmpData[0]);
							for(int iS=0;iS<nSteps;iS++){tmpData[iS] = real(v1c[iP][iJ][iS].c2);}
							nc_v1_y_re.putVar(startpB,countpB,&tmpData[0]);
							for(int iS=0;iS<nSteps;iS++){tmpData[iS] = real(v1c[iP][iJ][iS].c3);}
							nc_v1_z_re.putVar(startpB,countpB,&tmpData[0]);

							for(int iS=0;iS<nSteps;iS++){tmpData[iS] = imag(v1c[iP][iJ][iS].c1);}
							nc_v1_x_im.putVar(startpB,countpB,&tmpData[0]);
							for(int iS=0;iS<nSteps;iS++){tmpData[iS] = imag(v1c[iP][iJ][iS].c2);}
							nc_v1_y_im.putVar(startpB,countpB,&tmpData[0]);
							for(int iS=0;iS<nSteps;iS++){tmpData[iS] = imag(v1c[iP][iJ][iS].c3);}
							nc_v1_z_im.putVar(startpB,countpB,&tmpData[0]);

						}
				}
		
				vector<size_t> startp (1,0);
				vector<size_t> countp (1,nSteps);
		
				nc_t.putVar(startp,countp,&thisT[0]);
		
		
		}
				catch(exceptions::NcException &e) {
						cout << "NetCDF: unknown error" << endl;
						e.what();
						exit(1);
		}

		//cout << "DONE" << endl;
#endif

	} // End of xGrid loop

	// Write current(s) to file
	
	//cout << "Writing jP to file ... ";


////////////////////////////////////////////////////////////////////// COMMENTED OUT COMPUTING OF J
/*
	for(int iList=0;iList<nList;iList++) {

		stringstream ncjPFileName;
		ncjPFileName << "output/";
		// check directory exists
		struct stat st;
		if(stat(ncjPFileName.str().c_str(),&st) != 1) {
			int mkDirStat = mkdir(ncjPFileName.str().c_str(),S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
		}
		ncjPFileName << "/jP_";
		ncjPFileName << setw(3) << setfill('0') << iList;
		ncjPFileName << ".nc"; 	
#if DEBUGLEVEL >= 1
		cout<<ncjPFileName.str().c_str()<<endl;
#endif

		NcFile ncjPFile (ncjPFileName.str().c_str(), NcFile::replace);

		NcDim nc_nJp = ncjPFile.addDim("nJp", nJp);
		NcDim nc_scalar = ncjPFile.addDim("scalar", 1);

		NcVar nc_t = ncjPFile.addVar("t",ncFloat,nc_nJp);

        cout << "HELLO!"<< endl;

		NcVar nc_x = ncjPFile.addVar("x",ncFloat,nc_scalar);
		NcVar nc_freq = ncjPFile.addVar("freq",ncFloat,nc_scalar);

		NcVar nc_j1x = ncjPFile.addVar("j1x",ncFloat,nc_nJp);
		NcVar nc_j1y = ncjPFile.addVar("j1y",ncFloat,nc_nJp);
		NcVar nc_j1z = ncjPFile.addVar("j1z",ncFloat,nc_nJp);

		NcVar nc_j1xc_re = ncjPFile.addVar("j1xc_re",ncFloat,nc_scalar);
		NcVar nc_j1xc_im = ncjPFile.addVar("j1xc_im",ncFloat,nc_scalar);

		NcVar nc_j1yc_re = ncjPFile.addVar("j1yc_re",ncFloat,nc_scalar);
		NcVar nc_j1yc_im = ncjPFile.addVar("j1yc_im",ncFloat,nc_scalar);

		NcVar nc_j1zc_re = ncjPFile.addVar("j1zc_re",ncFloat,nc_scalar);
		NcVar nc_j1zc_im = ncjPFile.addVar("j1zc_im",ncFloat,nc_scalar);


/////////////////////////////////  NEED TO CHANGE THIS..... SAVES BASED ON iList////////////////////////////////

		nc_x.putVar(&rGrid[iList]);
		nc_freq.putVar(&freq);

		vector<size_t> startp (1,0);
		vector<size_t> countp (1,nJp);

#if LOWMEM >= 1
		float tmpJxRe = real(j1xc[iList]);
		float tmpJxIm = imag(j1xc[iList]);
		nc_j1xc_re.putVar(&tmpJxRe);
		nc_j1xc_im.putVar(&tmpJxIm);

		float tmpJyRe = real(j1yc[iList]);
		float tmpJyIm = imag(j1yc[iList]);
		nc_j1yc_re.putVar(&tmpJyRe);
		nc_j1yc_im.putVar(&tmpJyIm);

	    float tmpJzRe = real(j1zc[iList]);
		float tmpJzIm = imag(j1zc[iList]);
		nc_j1zc_re.putVar(&tmpJzRe);
		nc_j1zc_im.putVar(&tmpJzIm);
#else 
		nc_t.putVar(startp,countp,&tJp[0]);

		nc_j1x.putVar(startp,countp,&j1x[iX][0]);
		nc_j1y.putVar(startp,countp,&j1y[iX][0]);
		nc_j1z.putVar(startp,countp,&j1z[iX][0]);

		float tmpJxRe = real(j1xc[iX][0]);
		float tmpJxIm = imag(j1xc[iX][0]);
		nc_j1xc_re.putVar(&tmpJxRe);
		nc_j1xc_im.putVar(&tmpJxIm);

		float tmpJyRe = real(j1yc[iX][0]);
		float tmpJyIm = imag(j1yc[iX][0]);
		nc_j1yc_re.putVar(&tmpJyRe);
		nc_j1yc_im.putVar(&tmpJyIm);

	    float tmpJzRe = real(j1zc[iX][0]);
		float tmpJzIm = imag(j1zc[iX][0]);
		nc_j1zc_re.putVar(&tmpJzRe);
		nc_j1zc_im.putVar(&tmpJzIm);
#endif	
	}

	//ProfilerStop();

	cout << "DONE" << endl;
    */

#if CLOCK >= 1
        clock_t ProgramTime_ = clock();
        double ProgramTimeInSeconds = (ProgramTime_-ProgramTime)/(double) CLOCKS_PER_SEC;
#if defined(_OPENMP)
        ProgramTimeInSeconds = ProgramTimeInSeconds / nThreads;
        cout<<"nThreads: "<<nThreads<<endl;
#endif
        cout<<"Total Time [s]: "<<ProgramTimeInSeconds<<endl;
        cout<<"Total Time [m]: "<<ProgramTimeInSeconds/60.0<<endl;
        cout<<"Total Time [h]: "<<ProgramTimeInSeconds/3600.0<<endl;
#endif
	return EXIT_SUCCESS;
}
