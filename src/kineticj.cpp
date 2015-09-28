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
#include "read_gc_file.hpp"
#include "c3vec.hpp"
#include "read_e_field.hpp"
#include <new> // for stl::bad_alloc

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
                float dvx, dvy, dvz, d3v;
                float vPar, vPer, gyroPhase, u, vTh;
                float vAlp, vBet, phs;

				CParticle ();
				CParticle ( double _amu, int _Z);
				CParticle (float c1, float c2, float c3, 
								float v_c1, float v_c2, float v_c3, 
								double _amu, int _Z, float _weight );
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
	 double _amu, int _Z, float _weight ) :
	CSpecies (_amu, _Z) {

	c1 = _c1;
	c2 = _c2;
	c3 = _c3;

	v_c1 = _v_c1;
	v_c2 = _v_c2;
	v_c3 = _v_c3;

	weight = _weight;
	
	status = 0;
}

CParticle::CParticle ( CSpecies _species ) : CSpecies(_species) {
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

// This is not really a usefule magnitude
// and is only used to reduce a complex valued
// vector to a number for checking for Inf & NaNs
float mag ( const C3VecI &in ) {
        float c1 = abs(in.c1);
        float c2 = abs(in.c2);
        float c3 = abs(in.c3);
		return sqrt(pow(c1,2)+pow(c2,2)+pow(c3,2));
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

template <typename C3T1, typename C3T2>
C3VecI cross ( const C3T1 A, const C3T2 B ) {

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
    if(isnan(mag(arg))) answer = 1;
    return answer;
}

int isnan ( const C3VecI arg ) {
    int answer = 0;
    if( isnan( mag(arg) ) ) answer = 1;
    return answer;
}

int isinf ( const C3Vec arg ) {
    int answer = 0;
    if(isinf(mag(arg))) answer = 1;
    return answer;
}

int isinf ( const C3VecI arg ) {
    int answer = 0;
    if(isinf(mag(arg))) answer = 1;
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

// First-order orbits
C3Vec rk4_evalf ( CParticle &p, const float &t, const C3Vec &v, const C3Vec &x,
				const vector<C3Vec> &b0Vec, const vector<C3VecI> &e1, const float wrf ) {

	C3Vec b0(0,0,0), F;

	C3Vec v_x_b0 ( v.c2*b0.c3-v.c3*b0.c2, -1.0*(v.c1*b0.c3-v.c3*b0.c1), v.c1*b0.c2-v.c2*b0.c1); 

	return F*(p.q/p.m);	
}

template<class TYPE>
TYPE kj_interp1D ( const float &x, const vector<float> &xVec, const vector<TYPE> &yVec, int &status ) {

	float _x, x0, x1;
	float xTmp;
	xTmp = x;

#if _PARTICLE_BOUNDARY == 1
	if(x<xVec.front()||x>xVec.back()) {
			// Particle absorbing walls
#if DEBUG_INTERP >= 1
            if(xVec.size()!=yVec.size()) {
#if DEBUG_INTERP >= 2
                cout<<"ERROR: xVec and yVec are not the same size for interpolation!"<<endl;
#endif
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
			return TYPE(0);
	}
#elif _PARTICLE_BOUNDARY == 2
			// Periodic 
			if(xTmp<xVec.front()) xTmp = xVec.back()-(xVec.front()-xTmp);			
			if(xTmp>xVec.back()) xTmp = xVec.front()+(xTmp-xVec.back());			
#elif _PARTICLE_BOUNDARY == 3
			// Particle reflecting walls
			if(xTmp<xVec.front()) xTmp = xVec.front()+(xVec.front()-xTmp);			
			if(xTmp>xVec.back()) xTmp = xVec.back()-(xTmp-xVec.back());			
#endif

#if DEBUG_INTERP >= 1    
	if(status>0){
#if DEBUG_INTERP >= 2
			cout<<"ERROR: Should never get here with _PARTICLE_BOUNDARY ==2|3"<<endl;
#endif
			return TYPE(0);
	}
#endif
	_x = (xTmp-xVec.front())/(xVec.back()-xVec.front())*(xVec.size()-1);

	x0 = floor(_x);
	x1 = ceil(_x);
	
	// Catch for particle at point
	if(x0==x1) {
#if DEBUG_INTERP >= 2
        cout << "Particle version of kj_interp1D" << endl;
		cout << "x0: " << x0 << " x1: " <<x1<< " _x: "<<_x << endl;
		cout << "Particle at point catch: " << x0/x1 << "  "  << abs(1.0-x0/x1) << endl;
#endif
		return yVec[x0];
	}
	else {
		TYPE y0 = yVec[x0];
		TYPE y1 = yVec[x1];
#if DEBUG_INTERP >=2
        //cout << "kj_interp1D: " << endl;
        //if(typeid(TYPE)==typeid(C3Vec)) cout << "Type is C3Vec" << endl;
        //if(typeid(TYPE)==typeid(float)) cout << "Type is float" << endl;
        //kj_print(x0,"x0");
        //kj_print(x1,"x1");
        //kj_print(y0,"y0");
        //kj_print(y1,"y1");
        //cout << endl;
        if(x0>yVec.size()-1||x0<0||x1>yVec.size()-1||x1<1) {
                cout<<"ERROR: interpolation point off the end of array"<<endl;
                cout<<"x.front: "<<xVec.front()<<endl;
                cout<<"x.back: "<<xVec.back()<<endl;
                cout<<"x: "<<x<<endl;
                cout<<"xVec.size(): "<<xVec.size()<<endl;
                cout<<"yVec.size(): "<<yVec.size()<<endl;
                status=1;
                return TYPE(0);
        }
#endif
        TYPE result = y0+(_x-x0)*(y1-y0)/(x1-x0);

#if DEBUG_INTERP >=1
        if(isnan(result)) {
#if DEBUG_INTERP >= 2
                cout<<"ERROR: interpolation produced a NaN"<<endl;
#endif
                status=1;
                return TYPE(0);
        } 
        if(isinf(result)) {
#if DEBUG_INTERP >= 2
                cout<<"ERROR: interpolation produced a INF"<<endl;
#endif
                status = 1;
                return TYPE(0);
        }
#endif
		return result;
	}

}

C3Vec operator* ( const float A[][3], const C3Vec x ) {
        C3Vec B;
        B.c1 = A[0][0]*x.c1 + A[0][1]*x.c2 + A[0][2]*x.c3;
        B.c2 = A[1][0]*x.c1 + A[1][1]*x.c2 + A[1][2]*x.c3;
        B.c3 = A[2][0]*x.c1 + A[2][1]*x.c2 + A[2][2]*x.c3;
        return B;
}

C3VecI operator* ( const float A[][3], const C3VecI x ) {
        C3VecI B;
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

template <typename T>
T rot_CYL_to_XYZ ( const float t, const T vec, const int direction ) {

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

C3Vec rot_XYZ_to_abp ( const C3Vec A_XYZ, const C3Vec bUnit_XYZ, const int direction ) {

        // If direction<1 then the inverse rotation is applied, i.e., abp_to_XYZ

        C3Vec A_abp;

        C3Vec xu_xyz (1,0,0);
        C3Vec yu_xyz (0,1,0);
        C3Vec zu_xyz (0,0,1);

        C3Vec pu_xyz = bUnit_XYZ;

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


C3Vec getB_XYZ ( CParticle &p_XYZ, const vector<float> &rVec, const vector<C3Vec> &b0Vec_CYL ) {

	float _r = sqrt ( pow(p_XYZ.c1,2) + pow(p_XYZ.c2,2) );
	float _p = atan2 ( p_XYZ.c2, p_XYZ.c1 );

	C3Vec b0_CYL, b0_XYZ;

    int status = 0;
	b0_CYL = kj_interp1D ( _r, rVec, b0Vec_CYL, status );
    p_XYZ.status = max(p_XYZ.status,status);

    b0_XYZ = rot_CYL_to_XYZ ( _p, b0_CYL, 1 );

    return b0_XYZ;
}

C3VecI getE1orB1_XYZ ( CParticle &p_XYZ, const vector<float> &rVec, const vector<C3VecI> &E1Vec_CYL, int nPhi ) {

	float _r = sqrt ( pow(p_XYZ.c1,2) + pow(p_XYZ.c2,2) );
	float _p = atan2 ( p_XYZ.c2, p_XYZ.c1 );

	C3VecI E1_CYL, E1_XYZ;

    int status = 0;
	E1_CYL = kj_interp1D ( _r, rVec, E1Vec_CYL, status );
    p_XYZ.status = max(p_XYZ.status,status);

    complex<float> ii(0,1);

    E1_XYZ = exp(ii*float(nPhi*_p)) * rot_CYL_to_XYZ ( _p, E1_CYL, 1 ) ;

    return E1_XYZ;
}

// Zero-order orbits
C3Vec rk4_evalf ( CParticle &p, const float &t, 
				const C3Vec &v_XYZ, const C3Vec &x, const vector<C3Vec> &b0Vec_CYL,
			  	const vector<float> &rVec ) {

    C3Vec b0_XYZ;
    b0_XYZ = getB_XYZ(p,rVec,b0Vec_CYL);

    C3Vec v_x_b0 = cross(v_XYZ,b0_XYZ);

	return v_x_b0*(p.q/p.m);	
}

// Zero-order orbits
int rk4_move ( CParticle &p, const float &dt, const float &t0, 
				const vector<float> &r, const vector<C3Vec> &b0 ) {

		C3Vec yn0(p.v_c1,p.v_c2,p.v_c3), xn0(p.c1, p.c2, p.c3);
		C3Vec k1, k2, k3, k4, yn1, x1, x2, x3, x4, xn1; 

		k1 = rk4_evalf ( p, t0 + 0.0*dt, yn0         , xn0         , b0, r ) * dt;	
		x1 = yn0 * dt;

		k2 = rk4_evalf ( p, t0 + 0.5*dt, yn0 + 0.5*k1, xn0 + 0.5*x1, b0, r ) * dt;	
		x2 = (yn0 + 0.5*k1) * dt;

		k3 = rk4_evalf ( p, t0 + 0.5*dt, yn0 + 0.5*k2, xn0 + 0.5*x2, b0, r ) * dt;	
		x3 = (yn0 + 0.5*k2) * dt;

		k4 = rk4_evalf ( p, t0 + 1.0*dt, yn0 + 1.0*k3, xn0 + 1.0*x3, b0, r ) * dt;	
		x4 = (yn0 + 1.0*k3) * dt;

		yn1 = yn0 + 1.0/6.0 * (k1+2.0*k2+2.0*k3+k4) * (1-p.status); // the * (1-p.status) sets the move to zero for dead particles;
		xn1 = xn0 + 1.0/6.0 * (x1+2.0*x2+2.0*x3+x4) * (1-p.status);

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
			p.c1 = r.back()-(r.front()-p.c1);
		}
		if(p.c1>r.back())
		{
#if DEBUGLEVEL >= 1
			cout<<"Particle went right"<<endl;
#endif
			p.c1 = r.front()+(p.c1-r.back());
		}
#elif _PARTICLE_BOUNDARY == 3
		// Particle reflecting walls
		if(p.c1<r.front())
		{
			cout<<"Particle hit the left wall"<<endl;
			cout<<"r.front(): "<<r.front()<<endl;
			p.c1 = r.front()+(r.front()-p.c1);
			p.v_c1 = -p.v_c1;
		}
		if(p.c1>r.back())
		{
			cout<<"Particle hit the right wall"<<endl;
			cout<<"r.back(): "<<r.back()<<endl;
			p.c1 = r.back()-(p.c1-r.back());
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
        return p.status;
}

// Parallel acceleration
float eval_aPar ( CParticle &p, const C3Vec r, const vector<float> &r_GC, const vector<float> &bDotGradB ) {

    int status = 0;
    float This_bDotGradB = kj_interp1D ( r.c1, r_GC, bDotGradB, status );
    p.status = max(p.status,status);
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
float eval_vPer ( CParticle &p, const C3Vec r, const vector<float> &r_b0, const vector<C3Vec> &b0_CYL ) {

    int status = 0;         
	C3Vec This_b0_CYL = kj_interp1D ( r.c1, r_b0, b0_CYL, status );
    p.status = max(p.status,status);
    return sqrt ( 2.0 * p.u * mag(This_b0_CYL) / p.m );
}

// Guiding center veclocity
C3Vec eval_vGC ( CParticle &p, const C3Vec r, const float vPer, const float vPar, 
                const vector<float> &r_b0, const vector<C3Vec> &b0_CYL, 
                const vector<float> &r_GC, const vector<C3Vec> &curv_CYL, const vector<C3Vec> &grad_CYL ) {

    int status = 0; 
	C3Vec This_b0_CYL = kj_interp1D ( r.c1, r_b0, b0_CYL, status );
    p.status = max(p.status,status);
#if DEBUG_EVAL_VGC >= 1
    if(status>0) {
            cout<<"ERROR 1 in eval_vGC"<<endl;
            exit(1);
    }
#endif

    status = 0;
	C3Vec This_curv_CYL = kj_interp1D ( r.c1, r_GC, curv_CYL, status );
    p.status = max(p.status,status);

#if DEBUG_EVAL_VGC >= 1
    if(status>0) {
            cout<<"ERROR 2 in eval_vGC"<<endl;
            exit(1);
    }
#endif

    status = 0;
	C3Vec This_grad_CYL = kj_interp1D ( r.c1, r_GC, grad_CYL, status );
    p.status = max(p.status,status);
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

                C3Vec xn0_XYZ(p.c1, p.c2, p.c3);
                C3Vec xn0 = XYZ_to_CYL(xn0_XYZ);

		    	float This_vPer = eval_vPer ( p, xn0, r_b0, b0_CYL );
#if DEBUG_GC >= 2 
                cout << "p.vPer: " << p.vPer << endl;
                cout << "p.vPar: " << p.vPar << endl;
                cout << "This_vPer: " << This_vPer << endl;
                if(isnan(p.vPer)) exit(1);
#endif
		    	C3Vec This_vGC  = eval_vGC  ( p, xn0, This_vPer, p.vPar + 0, r_b0, b0_CYL, r_GC, curv_CYL, grad_CYL );
		    	float k1_vPar = dt * eval_aPar ( p, xn0, r_GC, bDotGradB );
		    	C3Vec k1_vgc  = dt * This_vGC;
#if DEBUG_GC >= 2
                kj_print(k1_vgc,"k1_vgc");
                kj_print(xn0,"xn0");
                cout<<"Status: "<<p.status<<endl;
                if(isnan(k1_vgc)||isinf(k1_vgc)||isnan(xn0)||isinf(xn0)||p.status>0) {
                        p.status = 1;
                        return p.status;
                }
#endif    
		    	This_vPer = eval_vPer ( p, xn0 + k1_vgc / 2.0, r_b0, b0_CYL );
		    	This_vGC  = eval_vGC  ( p, xn0 + k1_vgc / 2.0, This_vPer, p.vPar + k1_vPar / 2.0, r_b0, b0_CYL, r_GC, curv_CYL, grad_CYL );
		    	float k2_vPar = dt * eval_aPar ( p, xn0 + k1_vgc / 2.0, r_GC, bDotGradB ); 
		    	C3Vec k2_vgc  = dt * This_vGC;
#if DEBUG_GC >= 2
                kj_print(k2_vgc,"k2_vgc");
                if(isnan(k2_vgc)||isinf(k2_vgc)||isnan(xn0)||isinf(xn0)||p.status>0) {
                        p.status = 1;
                        return p.status;
                }
#endif 
		    	This_vPer = eval_vPer ( p, xn0 + k2_vgc / 2.0, r_b0, b0_CYL ); 
		    	This_vGC  = eval_vGC  ( p, xn0 + k2_vgc / 2.0, This_vPer, p.vPar + k2_vPar / 2.0, r_b0, b0_CYL, r_GC, curv_CYL, grad_CYL );
		    	float k3_vPar = dt * eval_aPar ( p, xn0 + k2_vgc / 2.0, r_GC, bDotGradB ); 
		    	C3Vec k3_vgc  = dt * This_vGC;
#if DEBUG_GC >= 2
                kj_print(k3_vgc,"k3_vgc");
                if(isnan(k3_vgc)||isinf(k3_vgc)||isnan(xn0)||isinf(xn0)||p.status>0) {
                        p.status = 1;
                        return p.status;
                }
#endif 
		    	This_vPer = eval_vPer ( p, xn0 + k3_vgc, r_b0, b0_CYL ); 
		    	This_vGC  = eval_vGC  ( p, xn0 + k3_vgc, This_vPer, p.vPar + k3_vPar, r_b0, b0_CYL, r_GC, curv_CYL, grad_CYL );
		    	float k4_vPar = dt * eval_aPar ( p, xn0 + k3_vgc, r_GC, bDotGradB );
		    	C3Vec k4_vgc  = dt * This_vGC; 
#if DEBUG_GC >= 2
                kj_print(k4_vgc,"k4_vgc");
                if(isnan(k4_vgc)||isinf(k4_vgc)||isnan(xn0)||isinf(xn0)||p.status>0) {
                        p.status = 1;
                        return p.status;
                }
#endif 	
		    	float vPar1 = p.vPar + ( k1_vPar + 2.0 * k2_vPar + 2.0 * k3_vPar + k4_vPar ) / 6.0 * (1-p.status);
		    	C3Vec xn1 = xn0 + ( k1_vgc + 2.0 * k2_vgc + 2.0 * k3_vgc + k4_vgc ) / 6.0* (1-p.status);

#if DEBUG_GC >=1 
                if(isnan(xn1)||isinf(xn1)) {
                        p.status = 1;
                        return p.status;
                }
#endif

                // Update particle with moved position and new vPar & vPer

		    	float vPer1 = eval_vPer ( p, xn1, r_b0, b0_CYL ); 

                p.vPar = vPar1;
                p.vPer = vPer1;

                C3Vec xn1_XYZ = CYL_to_XYZ(xn1);

                p.c1 = xn1_XYZ.c1; 
		        p.c2 = xn1_XYZ.c2;
		        p.c3 = xn1_XYZ.c3;

                // Update the XYZ velocity also

                int status = 0;
                C3Vec this_b0_CYL = kj_interp1D ( xn1.c1, r_b0, b0_CYL, status );
                p.status = max(p.status,status);

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

                return p.status;
}


// First-order orbits
int rk4_move ( CParticle &p, float dt, float t0, 
				const vector<C3Vec> &b0, const vector<C3VecI> &e1, const float wrf ) {

		C3Vec yn0(p.v_c1,p.v_c2,p.v_c3), xn0(p.c1, p.c2, p.c3);
		C3Vec k1, k2, k3, k4, yn1, x1, x2, x3, x4, xn1; 

		k1 = rk4_evalf ( p, t0 + 0.0*dt, yn0 + 0.*yn0, xn0         , b0, e1, wrf ) * dt;	
		x1 = k1 * dt;                                               
		k2 = rk4_evalf ( p, t0 + 0.5*dt, yn0 + 0.5*k1, xn0 + 0.5*x1, b0, e1, wrf ) * dt;	
		x2 = k2 * dt;                                               
		k3 = rk4_evalf ( p, t0 + 0.5*dt, yn0 + 0.5*k2, xn0 + 0.5*x2, b0, e1, wrf ) * dt;	
		x3 = k3 * dt;                                               
		k4 = rk4_evalf ( p, t0 + 1.0*dt, yn0 + 1.0*k3, xn0 + 1.0*x3, b0, e1, wrf ) * dt;	
		x4 = k4 * dt;

		yn1 = yn0 + 1.0/6.0 * (k1+2.0*k2+2.0*k3+k4) * (1-p.status); // the * (1-p.status) sets the move to zero for dead particles
		xn1 = xn0 + 1.0/6.0 * (x1+2.0*x2+2.0*x3+x4) * (1-p.status);

		p.c1 = xn1.c1;
		p.c2 = xn1.c2;
		p.c3 = xn1.c3;
		p.v_c1 = yn1.c1;
		p.v_c2 = yn1.c2;
		p.v_c3 = yn1.c3;

#if DEBUG_RK4 >= 1
        if(p.status!=0){
            cout<<"ERROR: p.status != 0"<<endl;
            //exit(1)
        }
#endif
        return p.status;
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
                int nPx, int nPy, int nPz, int nThermal, float &dv, vector<float> &r, vector<C3Vec> &b0_CYL) {

        vector<CParticle> pList;

        int nP = nPx * nPy * nPz;
        pList.resize(nP);

        float m = amu * _mi;
		float vTh = get_vTh ( amu, Z, T_keV );

# if DEBUG_MAXWELLIAN >= 1
        cout <<"amu: "<< amu<<endl;
        cout <<"Z: "<< Z<<endl;
        cout <<"vTh: "<< vTh<<endl;
#endif

		float vxRange = vTh * nThermal * 2;
		float vxMin	= -vxRange / 2.0;
		float dvx = vxRange / (nPx-1); 

		float vyRange = vTh * nThermal * 2;
		float vyMin	= -vyRange / 2.0;
		float dvy = vyRange / (nPy-1);

		float vzRange = vTh * nThermal * 2;
		float vzMin	= -vzRange / 2.0;
		float dvz = vzRange / (nPz-1);

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

    	            CParticle p (x,0,0,thisvx,thisvy,thisvz,amu,Z,weight);
                    pList[cnt] = p;
                    pList[cnt].number = cnt;
                    pList[cnt].vTh = vTh;

                    pList[cnt].d3v = dv;

                    // Get vPar, vPer and mu for guiding center integration

                    C3Vec thisV_XYZ(thisvx,thisvy,thisvz); 
                    int iStat = 0;
                    C3Vec this_b0_CYL = kj_interp1D(x,r,b0_CYL,iStat);
                    C3Vec this_b0_XYZ = rot_CYL_to_XYZ ( 0, this_b0_CYL, 1 );
                    float bMag = mag (this_b0_XYZ);
                    float vMag = mag (thisV_XYZ);

                    C3Vec thisV_abp = rot_XYZ_to_abp ( thisV_XYZ, this_b0_XYZ, 0 );

                    pList[cnt].vPar = thisV_abp.c3;
                    pList[cnt].vPer = sqrt(pow(thisV_abp.c1,2)+pow(thisV_abp.c2,2));
                    pList[cnt].gyroPhase = GetGyroPhase(thisV_abp); 
                    pList[cnt].u = pList[cnt].m * pow(pList[cnt].vPer,2) / ( 2.0 * bMag );

#if DEBUG_MAXWELLIAN >=2 
                    cout<<"ThisVx: "<<thisvx<<endl;
                    cout<<"ThisVy: "<<thisvy<<endl;
                    cout<<"ThisVz: "<<thisvz<<endl;
                    cout<<"b0_XYZ: "<<this_b0_XYZ.c1<<", "<<this_b0_XYZ.c2<<", "<<this_b0_XYZ.c3<<endl;
                    cout<<"vMag: "<<vMag<<endl;
                    cout<<"vPer: "<<pList[cnt].vPer<<endl;
                    cout<<"vPar: "<<pList[cnt].vPar<<endl;
                    cout<<"u: "<<pList[cnt].u<<endl<<endl;
                    if(isnan(pList[cnt].u)) exit(1);
                    if(vMag>3e8) exit(1);
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


// Calculate the jP given some know E and f(v)

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

        // Read config file

		libconfig::Config cfg;
		string cfgName = "kj.cfg";

        try {
            cfg.readFile(cfgName.c_str());
        }
        catch(const libconfig::FileIOException &fioex) {
            std::cerr << "I/O error while reading file." << std::endl;
            return(EXIT_FAILURE);
        }
        catch(const libconfig::ParseException &pex) {
            std::cerr << "Parse error at " << pex.getFile() << ":" << pex.getLine()
                << " - " << pex.getError() << std::endl;
            return(EXIT_FAILURE);
        }

	    int species_number = cfg.lookup("species_number");

		// Read E
		string eField_fName = cfg.lookup("eField_fName");	
		vector<C3Vec> e1Re_CYL, e1Im_CYL, b1Re_CYL, b1Im_CYL;
        vector<C3VecI> e1_CYL, b1_CYL;
		vector<C3Vec> b0_CYL, b0_XYZ;
        vector<float> r, n_m3;
		float freq;
        int eReadStat = read_e_field( eField_fName, species_number, freq, r, n_m3,  
                        e1_CYL, b1_CYL, e1Re_CYL, e1Im_CYL, b1Re_CYL, b1Im_CYL, b0_CYL);

        // Read GC terms
        string gc_fName = cfg.lookup("gc_fName");	
        vector<C3Vec> curv_CYL, grad_CYL;
        std::vector<float> r_gc, bDotGradB;
        int gcReadStat = read_gc_file( gc_fName, r_gc, curv_CYL, grad_CYL, bDotGradB );


	float wrf = freq * 2 * _pi;
	float xGridMin = cfg.lookup("xGridMin");
	float xGridMax = cfg.lookup("xGridMax");
	int nXGrid = cfg.lookup("nXGrid");
    cout<<"nXGrid: "<<nXGrid<<endl;

	vector<float> xGrid(nXGrid);
    vector<float> density_m3(nXGrid);
    vector<float> T_keV(nXGrid);
    vector<float> wrf_wc(nXGrid);
    vector<float> bMag_kjGrid(nXGrid);

	float xGridRng = 0;
	float xGridStep = 0;

	if(nXGrid>1) {
		xGridRng = xGridMax-xGridMin;
		xGridStep = xGridRng/(nXGrid-1);
	}

	for(int iX=0;iX<nXGrid;iX++) {
		xGrid[iX] = xGridMin+iX*xGridStep;
        int iStat;
        density_m3[iX] = kj_interp1D(xGrid[iX],r,n_m3,iStat);
        C3Vec this_b0= kj_interp1D(xGrid[iX],r,b0_CYL,iStat);
        bMag_kjGrid[iX] = mag(this_b0);
        T_keV[iX] = 2.0;//kj_interp1D(xGrid[iX],r,n_m3);
	}

    float MaxB0 = *max_element(bMag_kjGrid.begin(),bMag_kjGrid.end());

	//string googlePerfFileName = "/home/dg6/code/kineticj/googlep";
	//ProfilerStart(googlePerfFileName.c_str());
	//
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
    float amu = cfg.lookup("species_amu");
    float Z = cfg.lookup("species_Z");
    int nThermal = cfg.lookup("nThermal");
	int nP = nPx*nPy*nPz;
    float wc = Z*_e*MaxB0/(amu*_mi);
    float cyclotronPeriod = 2*_pi / wc;
    float dtMin 	= -cyclotronPeriod/nStepsPerCycle;
	int nSteps 		= nRFCycles*tRF/abs(dtMin)+1;

	for(int iX=0;iX<nXGrid;iX++) {
		float this_wc =	Z*_e*bMag_kjGrid[iX]/(amu*_mi);
		wrf_wc[iX] =  wrf / this_wc;
	}


#if PRINT_INFO >= 1
    cout << "dtMin [s]: " << dtMin << endl;
    cout << "Cyclotron Period: "<< cyclotronPeriod<<endl;
    cout << "RF Period: " << tRF << endl;
    cout << "nSteps: " << nSteps << endl;
    cout << "nStepsPerCycle: " << nStepsPerCycle << endl;
    cout << "freq: " << freq << endl;
    cout << "Max B0: " << MaxB0 << endl;
#endif
	
    vector<float> thisT;
    try {
	    thisT.resize(nSteps);
    }
    catch ( const std::bad_alloc &error ) {
        cout<<"Allocation error at "<<__FILE__<<__LINE__<<endl;
        cout<<error.what();
    }

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



	vector<vector<float> > j1x(nXGrid), j1y(nXGrid), j1z(nXGrid);
	vector<complex<float> > j1xc(nXGrid), j1yc(nXGrid), j1zc(nXGrid);

#if defined(_OPENMP)
        int nThreads, tid, spoken = 0;
#endif
 
	#pragma omp parallel for private(istat, tid, spoken)
	for(int iX=0;iX<nXGrid;iX++) {

#if defined(_OPENMP)
        nThreads = omp_get_num_threads();
        tid = omp_get_thread_num();
        if( tid == 0 && spoken == 0 ) {
            cout<<"tid : "<<tid<<endl;
            cout<<"OMP_NUM_THREADS: "<<nThreads<<endl;
            spoken = 1;
        }
#endif
        float dv;
        vector<CParticle> ThisParticleList(create_particles(xGrid[iX],amu,Z,T_keV[iX],density_m3[iX],
                                nPx,nPy,nPz,nThermal,dv,r,b0_CYL));

#if CLOCK >= 1
        clock_t startTime = clock();
#endif
        j1xc[iX] = complex<float>(0,0);
        j1yc[iX] = complex<float>(0,0);
        j1zc[iX] = complex<float>(0,0);

#if LOWMEM_USEPAPI >= 1
		cpuTime0=cpuTime;realTime0=realTime;flpIns0=flpIns;
		papiReturn = PAPI_flops ( &realTime, &cpuTime, &flpIns, &mFlops );
#endif	

		vector<float> f1(nP);
		//vector<complex<float> > f1xc(nP), f1yc(nP), f1zc(nP);
		vector<complex<float> > f1c(nP);

		for(int iP=0;iP<nP;iP++) {

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

            int write_iX = 75;
            int write_iP = 31;
            if(iX==write_iX && iP==write_iP) {
                cout<<"Write Particle Properties:"<<endl;
                cout<<" vTh: "<<thisParticle_XYZ.vTh<<endl;
                cout<<" v1: "<<thisParticle_XYZ.v_c1<<endl;
                cout<<" v2: "<<thisParticle_XYZ.v_c2<<endl;
                cout<<" v3: "<<thisParticle_XYZ.v_c3<<endl;

                OrbitFile.open("output/orbit.txt", ios::out | ios::trunc);
				OrbitFile<<"wc / wrf: "<< wrf_wc[iX]<<endl;
                OrbitFile<<" t  x  y  z  re(e1)  im(e1)  re(e2)  im(e2)  re(e3)  im(e3)  re(b1)  im(b1)  re(b2)  im(b2)  re(b3)  im(b3) status"<<endl;
                v1File.open("output/orbit_v1.txt", ios::out | ios::trunc);
                v1File<<" t  re(v11)  im(v11)  re(v12)  im(v12)  re(v13)  im(v13)"<<endl;
                e1_dot_grad_File.open("output/orbit_e1_dot_grad_df0_dv.txt", ios::out | ios::trunc);
                e1_dot_grad_File<<" t  re(v1xb01)  im(v1xb01)  re(v1xb02)  im(v1xb02)  re(v1xb03)  im(v1xb03)"<<endl;
                df0dv_File.open("output/df0dv.txt", ios::out | ios::trunc);
                df0dv_File<<" t  vx  vy  vz  valp  vbet  vpar  vper  gyroAngle  df0dv_x  df0dv_y  df0dv_z"<<endl;
            }
#endif    
			// generate orbit and get time-harmonic e along it

			vector<C3Vec> thisOrbit_XYZ(nSteps);
			vector<C3VecI> thisE1c_XYZ(nSteps,C3VecI());
			vector<C3VecI> thisB1c_XYZ(nSteps,C3VecI());
			C3VecI thisV1c_(0,0,0), dVc(0,0,0), crossTerm(0,0,0);
            vector<complex<float> > this_e1_dot_gradvf0(nSteps);
            vector<C3VecI> this_vCrossB1(nSteps);

	 		for(int i=0;i<nSteps;i++) {	
#if DEBUG_MOVE >=1 
                cout << "Position Before Move: " << thisParticle_XYZ.c1 << "  " << thisParticle_XYZ.c2 << "  " << thisParticle_XYZ.c3 << endl;
                cout << "p.status: "<<thisParticle_XYZ.status << endl;
#endif
				thisOrbit_XYZ[i] = C3Vec(thisParticle_XYZ.c1,thisParticle_XYZ.c2,thisParticle_XYZ.c3);
#if GC_ORBITS >=1 
                int MoveStatus = rk4_move_gc ( thisParticle_XYZ, dtMin, thisT[i], 
                                r, b0_CYL, r_gc, curv_CYL, grad_CYL, bDotGradB, wrf );
#else
				int MoveStatus = rk4_move ( thisParticle_XYZ, dtMin, thisT[i], r, b0_CYL );
#endif
                int OverallStatus = max(thisParticle_XYZ.status,MoveStatus);
#if DEBUG_MOVE >=1 
                if(MoveStatus>0) {
                    cout << "Position After Move: " << thisParticle_XYZ.c1 << "  " << thisParticle_XYZ.c2 << "  " << thisParticle_XYZ.c3 << endl;
                    cout<<"ERROR: rk4_move* threw an error"<<endl;
                    cout<<"MoveStatus: "<<MoveStatus<<endl;
                    exit(1); 
                }
#endif
            
                C3Vec thisPos(thisParticle_XYZ.c1,thisParticle_XYZ.c2,thisParticle_XYZ.c3);
                C3Vec thisVel_XYZ(thisParticle_XYZ.v_c1,thisParticle_XYZ.v_c2,thisParticle_XYZ.v_c3);
				C3Vec thisB0 = kj_interp1D ( thisOrbit_XYZ[i].c1, r, b0_CYL, istat );
#if GC_ORBITS >= 1
                thisVel_XYZ = thisB0 / mag(thisB0) * thisParticle_XYZ.vPar; // vPar vector in XYZ
                cout<<thisParticle_XYZ.vPar<<"  "<<thisParticle_XYZ.vPer<<endl;
                kj_print(thisVel_XYZ,"thisVel_XYZ");
				C3Vec gradv_f0_XYZ = maxwellian_df0_dv ( thisVel_XYZ, T_keV[iX], density_m3[iX], thisParticle_XYZ.amu, thisParticle_XYZ.Z );
#else
				C3Vec gradv_f0_XYZ = maxwellian_df0_dv ( thisVel_XYZ, T_keV[iX], density_m3[iX], thisParticle_XYZ.amu, thisParticle_XYZ.Z );
#endif

                C3VecI E1_XYZ;
                complex<float> _i(0,1);
                E1_XYZ = hanningWeight[i] * exp(-_i*wrf*thisT[i]) * getE1orB1_XYZ(thisParticle_XYZ,r,e1_CYL,nPhi);
                thisE1c_XYZ[i] = E1_XYZ * (1-thisParticle_XYZ.status);

	            C3VecI B1_XYZ;
                B1_XYZ = hanningWeight[i] * exp(-_i*wrf*thisT[i]) * getE1orB1_XYZ(thisParticle_XYZ,r,b1_CYL,nPhi);
                thisB1c_XYZ[i] = B1_XYZ * (1-thisParticle_XYZ.status);
	
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

#if GC_ORBITS >= 1
                // For GC orbits (electrons) use only the orbit parallel piece, 
                // since the perp peice will cancel due to E not varying within
                // a cyclotron period.
                //
                // NO - WE DO NOT HAVE thisVel_XYZ for GC!!!!!
                C3Vec orbitParallelUnitVector_XYZ = thisB0 / mag(thisB0);
                //kj_print(orbitParallelUnitVector_XYZ,"unit");
                complex<float> orbitParallel_E = dot(thisE1c_XYZ[i], orbitParallelUnitVector_XYZ);
                float orbitParallel_gradv_f0 = dot(gradv_f0_XYZ, orbitParallelUnitVector_XYZ);
                //cout<<"E : "<<orbitParallel_E<<" gf0: "<<orbitParallel_gradv_f0<<endl;
                this_e1_dot_gradvf0[i] = orbitParallel_E * orbitParallel_gradv_f0;
                //cout<<this_e1_dot_gradvf0[i]<<" e1dotgrad"<<endl;
                complex<float> _full = dot(thisE1c_XYZ[i], gradv_f0_XYZ);
                //cout<<"_full : "<<_full<<endl;
#else
                this_vCrossB1[i] = cross(thisVel_XYZ,thisB1c_XYZ[i]);
                C3VecI this_force = this_vCrossB1[i] + thisE1c_XYZ[i];
                this_e1_dot_gradvf0[i] = dot(this_force, gradv_f0_XYZ);
#endif

#if LOWMEM_ORBIT_WRITE >= 1
                if(iX==write_iX && iP==write_iP) {
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
 
                if(iX==write_iX && iP==write_iP) {
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
                            <<"    "<< real(this_vCrossB1[i].c1)
                            <<"    "<< imag(this_vCrossB1[i].c1)
                            <<"    "<< real(this_vCrossB1[i].c2)
                            <<"    "<< imag(this_vCrossB1[i].c2)
                            <<"    "<< real(this_vCrossB1[i].c3)
                            <<"    "<< imag(this_vCrossB1[i].c3)
                            <<"    "<< thisParticle_XYZ.status 
                            << endl;
                }
                if(iX==write_iX && iP==write_iP) {
                    e1_dot_grad_File<<scientific;
                    e1_dot_grad_File<< thisT[i]
                            <<"    "<< real(this_e1_dot_gradvf0[i])
                            <<"    "<< imag(this_e1_dot_gradvf0[i])
                            << endl;
                }
#endif
			}
#if LOWMEM_ORBIT_WRITE >= 1
            if(iX==write_iX && iP==write_iP) {
                OrbitFile.close();
            }
#endif
			complex<float> this_f1c = -qOverm * intVecArray ( thisT, this_e1_dot_gradvf0 );

#if LOWMEM_ORBIT_WRITE >= 1
            if(iX==write_iX && iP==write_iP) {

                complex<float> tmp = 0.0;
                for(int i=0;i<nSteps;i++) {
                        tmp += -qOverm * this_e1_dot_gradvf0[i] * dtMin;
                        v1File<<thisT[i]
                                <<"    "<< real(tmp)
                                <<"    "<< imag(tmp)
                                << endl;
                }
            }
#endif
			f1c[iP] = -this_f1c;

			float v0x_i = ThisParticleList[iP].v_c1;
			float v0y_i = ThisParticleList[iP].v_c2;
			float v0z_i = ThisParticleList[iP].v_c3;

			float h = dv * Ze;

			#pragma omp critical // "atomic" does not work for complex numbers
			{
				j1xc[iX] += h * ( v0x_i*f1c[iP] ); 
				j1yc[iX] += h * ( v0y_i*f1c[iP] ); 
				j1zc[iX] += h * ( v0z_i*f1c[iP] ); 
			}
		}

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

	for(int iX=0;iX<nXGrid;iX++) {

		stringstream ncjPFileName;
		ncjPFileName << "output/";
		// check directory exists
		struct stat st;
		if(stat(ncjPFileName.str().c_str(),&st) != 1) {
			int mkDirStat = mkdir(ncjPFileName.str().c_str(),S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
		}
		ncjPFileName << "/jP_";
		ncjPFileName << setw(3) << setfill('0') << iX;
		ncjPFileName << ".nc"; 	
#if DEBUGLEVEL >= 1
		cout<<ncjPFileName.str().c_str()<<endl;
#endif

		NcFile ncjPFile (ncjPFileName.str().c_str(), NcFile::replace);

		NcDim nc_nJp = ncjPFile.addDim("nJp", nJp);
		NcDim nc_scalar = ncjPFile.addDim("scalar", 1);

		NcVar nc_t = ncjPFile.addVar("t",ncFloat,nc_nJp);

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

		nc_x.putVar(&xGrid[iX]);
		nc_freq.putVar(&freq);

		vector<size_t> startp (1,0);
		vector<size_t> countp (1,nJp);

		float tmpJxRe = real(j1xc[iX]);
		float tmpJxIm = imag(j1xc[iX]);
		nc_j1xc_re.putVar(&tmpJxRe);
		nc_j1xc_im.putVar(&tmpJxIm);

		float tmpJyRe = real(j1yc[iX]);
		float tmpJyIm = imag(j1yc[iX]);
		nc_j1yc_re.putVar(&tmpJyRe);
		nc_j1yc_im.putVar(&tmpJyIm);

	    float tmpJzRe = real(j1zc[iX]);
		float tmpJzIm = imag(j1zc[iX]);
		nc_j1zc_re.putVar(&tmpJzRe);
		nc_j1zc_im.putVar(&tmpJzIm);
	
	}

	//ProfilerStop();

	cout << "DONE" << endl;

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
