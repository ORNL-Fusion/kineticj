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

#if USEPAPI >= 1
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

class C3VecI {
		public:
				complex<float> c1, c2, c3;

				C3VecI () {c1=complex<float>(0.0f,0.0f);c2=complex<float>(0.0f,0.0f);c3=complex<float>(0.0f,0.0f);};
				C3VecI ( complex<float> _c1, complex<float> _c2, complex<float> _c3 ) {c1=_c1;c2=_c2;c3=_c3;};
};

class C3Vec {
		public:
				float c1, c2, c3;

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

C3Vec C3Vec::operator+ (const C3Vec &other) {
		//return C3Vec(*this)+=other;
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

// Global (not member) functions for lhs operators

C3Vec operator* ( const float &other, const C3Vec &rhs ) {
		return C3Vec(rhs)*=other;
}

C3Vec operator+ ( const C3Vec &other, const C3Vec &rhs) {
		return C3Vec(other.c1+rhs.c1,other.c2+rhs.c2,other.c3+rhs.c3);
}

// First-order orbits
C3Vec rk4_evalf ( CParticle &p, const float &t, const C3Vec &v, const C3Vec &x,
				const vector<C3Vec> &b0Vec, const vector<C3VecI> &e1, const float wrf ) {

	C3Vec b0(0,0,0), F;

	C3Vec v_x_b0 ( v.c2*b0.c3-v.c3*b0.c2, -1.0*(v.c1*b0.c3-v.c3*b0.c1), v.c1*b0.c2-v.c2*b0.c1); 
	//C3Vec F ( real(e1.c1) * cos ( wrf * t ) + imag(e1.c1) * sin ( wrf * t ) + v_x_b0.c1,
	//	  	real(e1.c2) * cos ( wrf * t ) + imag(e1.c2) * sin ( wrf * t ) + v_x_b0.c2,
	//	  	real(e1.c3) * cos ( wrf * t ) + imag(e1.c3) * sin ( wrf * t ) + v_x_b0.c3 );

	return F*(p.q/p.m);	
}

C3Vec kj_interp1D ( const float &x, const vector<float> &xVec, const vector<C3Vec> &yVec, int &stat ) {

	if(x<xVec.front()||x>xVec.back()||stat>0) {
			//printf("\t%s line: %i\n",__FILE__,__LINE__);
			//cout << "\tERROR: Interpolation pt off grid." << endl;
			//cout << "\tx: " << x << " x.front(): " << xVec.front() << " x.back(): " << xVec.back() << endl;
			++stat;
			return C3Vec(0,0,0);
	}

	float _x = (x-xVec.front())/(xVec.back()-xVec.front())*(xVec.size()-1);
	float x0 = floor(_x);
	float x1 = ceil(_x);

	// Catch for particle at point
	if(x0==x1) {
		//cout << "x0: " << x0 << " x1: " <<x1<< " _x: "<<_x << endl;
		//cout << "Particle at point catch: " << x0/x1 << "  "  << abs(1.0-x0/x1) << endl;
		return yVec[x0];
	}
	else {
		C3Vec y0 = yVec[x0];
		C3Vec y1 = yVec[x1];

		return y0+(_x-x0)*(y1-y0)/(x1-x0);
	}

}

// Zero-order orbits
C3Vec rk4_evalf ( CParticle &p, const float &t, 
				const C3Vec &v_XYZ, const C3Vec &x, const vector<C3Vec> &b0Vec_CYL,
			  	const vector<float> &rVec ) {

	// Interpolate b0 at location in CYL
	
	float _r = sqrt ( pow(x.c1,2) + pow(x.c2,2) );
	float _p = atan2 ( x.c2, x.c1 );

#if DEBUGLEVEL >= 3
	cout << "\t\t\tx: " << x.c1 << " y: " << x.c2 << " z: " << x.c3 << endl;
	cout << "\t\t\tr: " << _r << " p: " << _p << endl;
	cout << "\t\t\trVec.front(): " << rVec.front() << endl;
	cout << "\t\t\tv_XYZ: " << v_XYZ.c1 << "  " << v_XYZ.c2 << "  " << v_XYZ.c3 << endl;
#endif

	C3Vec b0_CYL, b0_XYZ;

	b0_CYL = kj_interp1D ( _r, rVec, b0Vec_CYL, p.status );

	b0_XYZ = C3Vec( cos(_p)*b0_CYL.c1-sin(_p)*b0_CYL.c2+0,
					sin(_p)*b0_CYL.c1+cos(_p)*b0_CYL.c2+0,
					0+0+1*b0_CYL.c3 );

	C3Vec v_x_b0 ( v_XYZ.c2*b0_XYZ.c3-v_XYZ.c3*b0_XYZ.c2, 
					-1.0*(v_XYZ.c1*b0_XYZ.c3-v_XYZ.c3*b0_XYZ.c1), 
					v_XYZ.c1*b0_XYZ.c2-v_XYZ.c2*b0_XYZ.c1);

#if DEBUGLEVEL >= 3
	cout << "\tvxb0: " << v_x_b0.c1 << "  " << v_x_b0.c2 << "  " << v_x_b0.c3 << endl;
	cout << "\tp.q/p.m: " << p.q/p.m << endl;
#endif

	return v_x_b0*(p.q/p.m);	
}

// Zero-order orbits
void rk4_move ( CParticle &p, const float &dt, const float &t0, 
				const vector<C3Vec> &b0, const vector<float> &r ) {

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

		yn1 = yn0 + 1.0/6.0 * (k1+2.0*k2+2.0*k3+k4);
		xn1 = xn0 + 1.0/6.0 * (x1+2.0*x2+2.0*x3+x4);

		p.c1 = xn1.c1;
		p.c2 = xn1.c2;
		p.c3 = xn1.c3;
		p.v_c1 = yn1.c1;
		p.v_c2 = yn1.c2;
		p.v_c3 = yn1.c3;

#if DEBUGLEVEL >= 3
		cout << "\tx0_XYZ: " << xn0.c1 << "  " << xn0.c2 << "  " << xn0.c3 << endl;
		cout << "\tv0_XYZ: " << yn0.c1 << "  " << yn0.c2 << "  " << yn0.c3 << endl;
		cout << "\tx1_XYZ: " << xn1.c1 << "  " << xn1.c2 << "  " << xn1.c3 << endl;
		cout << "\tv1_XYZ: " << yn1.c1 << "  " << yn1.c2 << "  " << yn1.c3 << endl;
		cout << "\tE: " << 0.5 * p.m * sqrt (pow(p.v_c1,2)+pow(p.v_c2,2)+pow(p.v_c3,2))/_e << endl;
#endif

}

// First-order orbits
void rk4_move ( CParticle &p, float dt, float t0, 
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

		yn1 = yn0 + 1.0/6.0 * (k1+2.0*k2+2.0*k3+k4);
		xn1 = xn0 + 1.0/6.0 * (x1+2.0*x2+2.0*x3+x4);

		p.c1 = xn1.c1;
		p.c2 = xn1.c2;
		p.c3 = xn1.c3;
		p.v_c1 = yn1.c1;
		p.v_c2 = yn1.c2;
		p.v_c3 = yn1.c3;
}


float maxC3VecAbs ( const vector<C3Vec> &input ) {

	vector<float> inputAbs(input.size());
	for(int i=0;i<input.size();i++) {
		inputAbs[i] = sqrt(pow(input[i].c1,2)+pow(input[i].c2,2)+pow(input[i].c3,2));
	}
	return *max_element(inputAbs.begin(),inputAbs.end());
}

C3Vec intC3VecArray ( const vector<float> &x, const vector<C3Vec> &f ) {

	C3Vec result(0,0,0);
	for(int i=1;i<f.size()-1;i++) {
		float h = x[i+1]-x[i];
		result += h/2.0*(f[i]+f[i+1]);
	}

	return result;
}


// Calculate the jP given some know E and f(v)

int main ( int argc, char **argv )
{

#if USEPAPI >= 1
		float realTime0, cpuTime0, realTime=0, cpuTime=0, mFlops;
		long_long flpIns0, flpIns=0;
		int papiReturn;

		cpuTime0=cpuTime;realTime0=realTime;flpIns0=flpIns;
		papiReturn = PAPI_flops ( &realTime, &cpuTime, &flpIns, &mFlops );
		printf("Real_time:\t%f\nProc_time:\t%f\nTotal flpins:\t%lld\nMFLOPS:\t\t%f\n",
					   realTime-realTime0, cpuTime-cpuTime0, flpIns-flpIns0, mFlops);
#endif

		libconfig::Config cfg;
		string cfgName = "kj.cfg";

		// Write a config file if required
		//libconfig::Setting &root = cfg.getRoot();
		//root.add("xGridMin", libconfig::Setting::TypeFloat) = 98.0;
		//root.add("xGridMax", libconfig::Setting::TypeFloat) = 100.0;
		//root.add("nXGrid", libconfig::Setting::TypeInt) = 20;
		//root.add("nRFCycles", libconfig::Setting::TypeInt) = 10;
		//root.add("nStepsPerCycle", libconfig::Setting::TypeInt) = 100;
		//root.add("nJpCycles", libconfig::Setting::TypeInt) = 6;
		//root.add("nJpPerCycle", libconfig::Setting::TypeInt) = 20;
		//root.add("eField_fName", libconfig::Setting::TypeString) = "data/kj_aorsa_1d.nc";
		//root.add("particleList_fName", libconfig::Setting::TypeString) = "data/f.nc";
		//root.add("runIdent", libconfig::Setting::TypeString) = "thisRun";
		//cfg.writeFile(cfgName.c_str());
		
		// Open the config file
		cfg.readFile(cfgName.c_str());

		string runIdent = cfg.lookup("runIdent");	

		// Read E
		string eField_fName = cfg.lookup("eField_fName");	
		cout << "Reading eField data file " << eField_fName << endl;

		// Here we are using the cxx-4 netcdf interface by Lynton Appel
		// This needs netCDF 4.1.1 or later build with
		// ./configure --enable-cxx-4 [plus other options]

		vector<float> r, b0_r, b0_p, b0_z,
				e_r_re, e_p_re, e_z_re,
				e_r_im, e_p_im, e_z_im;
		vector<C3Vec> b0_CYL, b0_XYZ;
		
		float freq;

		vector<complex<float> > e_r, e_p, e_z;	

		ifstream file(eField_fName.c_str());
		if(!file.good()) {
			cout << "ERROR: Cannot find file " << eField_fName << endl;
			exit(1);
		}


		try {
				NcFile dataFile ( eField_fName.c_str(), NcFile::read );
	
				NcDim nc_nR(dataFile.getDim("nR"));
				NcDim nc_scalar(dataFile.getDim("scalar"));
	
				int nR = nc_nR.getSize();
	
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

				nc_r.getVar(&r[0]);
				nc_freq.getVar(&freq);

				nc_b0_r.getVar(&b0_r[0]);
				nc_b0_p.getVar(&b0_p[0]);
				nc_b0_z.getVar(&b0_z[0]);

				b0_CYL.resize(nR);
				b0_XYZ.resize(nR);
				for(int i=0; i<nR; i++) {
						b0_CYL[i] = C3Vec(b0_r[i],b0_p[i],b0_z[i]);
						b0_XYZ[i] = C3Vec(cos(0.0)*b0_CYL[i].c1-sin(0.0)*b0_CYL[i].c2+0,
										sin(0.0)*b0_CYL[i].c1+cos(0.0)*b0_CYL[i].c2+0,
										0+0+1*b0_CYL[i].c3);
				}

				nc_e_r_re.getVar(&e_r_re[0]);
				nc_e_p_re.getVar(&e_p_re[0]);
				nc_e_z_re.getVar(&e_z_re[0]);
				nc_e_r_im.getVar(&e_r_im[0]);
				nc_e_p_im.getVar(&e_p_im[0]);
				nc_e_z_im.getVar(&e_z_im[0]);

				for(int i=0; i<nR; i++){
						e_r.push_back(complex<float>( e_r_re[i], e_r_im[i] ) );
						e_p.push_back(complex<float>( e_p_re[i], e_p_im[i] ) );
						e_z.push_back(complex<float>( e_z_re[i], e_z_im[i] ) );
				}

				cout << "\tR[0]: " << r[0] << ", R["<<nR<<"]: " << r[r.size()-1] << endl;
				cout << "\tfreq: " << freq << endl;
				vector<float>::iterator min = min_element(b0_p.begin(),b0_p.end());
				vector<float>::iterator max = max_element(b0_p.begin(),b0_p.end());
				cout << "\tmin(b0_p): " << *min << endl;
				cout << "\tmax(b0_p): " << *max << endl;
				cout << "\tabs(e_r[nR/2]): " << abs(e_r[nR/2]) << endl;
				cout << "\tabs(e_p[nR/2]): " << abs(e_p[nR/2]) << endl;
				cout << "\tabs(e_z[nR/2]): " << abs(e_z[nR/2]) << endl;
		}
		catch(exceptions::NcException &e) {
				cout << "NetCDF: unknown error" << endl;
				e.what();
		}

		// Rotate the e field to XYZ

		vector<C3Vec> e1Re_CYL, e1Im_CYL;
		e1Re_CYL.resize(e_r.size());
		e1Im_CYL.resize(e_r.size());

		vector<C3Vec> e1Re_XYZ, e1Im_XYZ;
		e1Re_XYZ.resize(e_r.size());
		e1Im_XYZ.resize(e_r.size());

		for(int i=0;i<e_r.size();i++) {

			e1Re_CYL[i].c1 = real(e_r[i]);
			e1Re_CYL[i].c2 = real(e_p[i]);
			e1Re_CYL[i].c3 = real(e_z[i]);
			e1Im_CYL[i].c1 = imag(e_r[i]);
			e1Im_CYL[i].c2 = imag(e_p[i]);
			e1Im_CYL[i].c3 = imag(e_z[i]);

			float _p = 0.0;

			e1Re_XYZ[i] = C3Vec( 
				cos(_p)*e1Re_CYL[i].c1-sin(_p)*e1Re_CYL[i].c2+0,
				sin(_p)*e1Re_CYL[i].c1+cos(_p)*e1Re_CYL[i].c2+0,
				0+0+1*e1Re_CYL[i].c3 );

			e1Im_XYZ[i] = C3Vec( 
				cos(_p)*e1Im_CYL[i].c1-sin(_p)*e1Im_CYL[i].c2+0,
				sin(_p)*e1Im_CYL[i].c1+cos(_p)*e1Im_CYL[i].c2+0,
				0+0+1*e1Im_CYL[i].c3 );
		}


		// Read particle list
		string particleList_fName = cfg.lookup ("particleList_fName");	
		cout << "Reading particle list " << particleList_fName << endl;

		ifstream file2(particleList_fName.c_str());
		if(!file.good()) {
			cout << "ERROR: Cannot find file " << particleList_fName << endl;
			exit(1);
		}

		vector<float> p_x, p_y, p_z, p_vx, p_vy, p_vz, p_amu, p_weight;
		vector<int> p_Z;
		float vTh;
		int nThermal;
		
		try {
				NcFile dataFile ( particleList_fName.c_str(), NcFile::read );
	
				NcDim nc_nP(dataFile.getDim("nP"));
	
				int nP = nc_nP.getSize();
	
				cout << "\tnP: " << nP << endl;

				NcVar nc_p_amu(dataFile.getVar("amu"));
				NcVar nc_p_Z(dataFile.getVar("Z"));

				NcVar nc_p_x(dataFile.getVar("x"));
				NcVar nc_p_y(dataFile.getVar("y"));
				NcVar nc_p_z(dataFile.getVar("z"));
				
				NcVar nc_p_vx(dataFile.getVar("vx"));
				NcVar nc_p_vy(dataFile.getVar("vy"));
				NcVar nc_p_vz(dataFile.getVar("vz"));

				NcVar nc_p_weight(dataFile.getVar("weight"));

				NcVar nc_nThermal(dataFile.getVar("nThermal"));
				NcVar nc_vTh(dataFile.getVar("vTh"));

				p_x.resize(nP);
				p_y.resize(nP);
				p_z.resize(nP);

				p_vx.resize(nP);
				p_vy.resize(nP);
				p_vz.resize(nP);

				p_weight.resize(nP);

				p_amu.resize(nP);
				p_Z.resize(nP);

				nc_p_x.getVar(&p_x[0]);
				nc_p_y.getVar(&p_y[0]);
				nc_p_z.getVar(&p_z[0]);

				nc_p_vx.getVar(&p_vx[0]);
				nc_p_vy.getVar(&p_vy[0]);
				nc_p_vz.getVar(&p_vz[0]);

				nc_p_weight.getVar(&p_weight[0]);

				nc_p_amu.getVar(&p_amu[0]);
				nc_p_Z.getVar(&p_Z[0]);

				nc_nThermal.getVar(&nThermal);
				nc_vTh.getVar(&vTh);

		}
		catch(exceptions::NcException &e) {
				cout << "NetCDF: unknown error" << endl;
				e.what();
				exit(1);
		}

		vector<CParticle> particles_XYZ;
		particles_XYZ.resize(p_x.size());

		for(int i=0;i<particles_XYZ.size();i++){

				CParticle thisParticle (p_x[i],p_y[i],p_z[i],
								p_vx[i],p_vy[i],p_vz[i],
								p_amu[i],p_Z[i],p_weight[i]);
				particles_XYZ[i] = thisParticle;
				particles_XYZ[i].number = i;
		}

	// Langmuir wave dispersion relation

	double wrf = freq * 2 * _pi;
	double wpe = sqrt ( 1e14*pow(_e,2)/(_me*_e0) );
	double kParSq = (pow(wrf,2)-pow(wpe,2))/(2*pow(vTh,2));
	double kPar = sqrt ( kParSq );
	double lambdaPar = 2*_pi/kPar;
	
	cout << "kParSq [m^-2]: " << kParSq << endl;
	cout << "kPar [m^-1]: " << kPar << endl;
	cout << "lambdaPar [m]: " << lambdaPar << endl;

	float xGridMin = cfg.lookup("xGridMin");
	float xGridMax = cfg.lookup("xGridMax");
	int nXGrid = cfg.lookup("nXGrid");
	vector<float> xGrid(nXGrid);
	float xGridRng = 0;
	float xGridStep = 0;
	
	if(nXGrid>1) {
		xGridRng = xGridMax-xGridMin;
		xGridStep = xGridRng/(nXGrid-1);
	}

	for(int iX=0;iX<nXGrid;iX++) {
		xGrid[iX] = xGridMin+iX*xGridStep;
	}

	vector<CParticle> particles_XYZ_0(particles_XYZ);

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

	int nRFCycles 		= cfg.lookup("nRFCycles");
	int nStepsPerCycle 	= cfg.lookup("nStepsPerCycle"); 
	float tRF 			= (2*_pi)/wrf;
	float dtMin 		= -tRF/nStepsPerCycle;
	int nSteps 			= nRFCycles*nStepsPerCycle+1;
	int nJpCycles 		= cfg.lookup("nJpCycles");
	int nJpPerCycle 	= cfg.lookup("nJpPerCycle");
	int nJp 			= nJpCycles * nJpPerCycle + 1;
	float dtJp 			= tRF / nJpPerCycle;
	int istat = 0;

	vector<float> tJp(nJp,0);
	for(int jt=0;jt<nJp;jt++) {
		tJp[jt] = jt*dtJp;
	}

	vector<float> thisT(nSteps);
	for(int i=0;i<nSteps;i++) {	
		thisT[i]=i*dtMin;
	}
	
	for(int iX=0;iX<nXGrid;iX++) {

		vector<float> j1x(nJp,0), j1y(nJp,0), j1z(nJp,0);//, tJ(nJp,0);

		cout << "xGrid " << iX << endl;

#if USEPAPI >= 1
		cpuTime0=cpuTime;realTime0=realTime;flpIns0=flpIns;
		papiReturn = PAPI_flops ( &realTime, &cpuTime, &flpIns, &mFlops );
#endif	

#if LOWMEM >= 1 // START OF THE LOWMEM CODING vvv

		#pragma omp parallel for private(istat)
		for(int iP=0;iP<particles_XYZ.size();iP++) {

			vector<C3Vec> thisOrbitE_re_XYZ(nSteps,C3Vec(0,0,0));
			vector<C3Vec> thisOrbitE_im_XYZ(nSteps,C3Vec(0,0,0));

			CParticle thisParticle_XYZ(particles_XYZ[iP]);
			thisParticle_XYZ.c1 = xGrid[iX];

			// generate orbit and get time-harmonic e along it

			vector<C3Vec> thisOrbit_XYZ(nSteps);

	 		for(int i=0;i<nSteps;i++) {	

				//thisT[i]=i*dtMin;
				if(thisParticle_XYZ.status==0) {

					thisOrbit_XYZ[i] = C3Vec(thisParticle_XYZ.c1,thisParticle_XYZ.c2,thisParticle_XYZ.c3);
					rk4_move ( thisParticle_XYZ, dtMin, thisT[i], b0_CYL, r );

					if(thisParticle_XYZ.status==0) {
						istat = 0;
						C3Vec e1ReTmp_XYZ = kj_interp1D ( thisOrbit_XYZ[i].c1, r, e1Re_XYZ, istat );
						istat = 0;
						C3Vec e1ImTmp_XYZ = kj_interp1D ( thisOrbit_XYZ[i].c1, r, e1Im_XYZ, istat );
						thisOrbitE_re_XYZ[i] = e1ReTmp_XYZ;
						thisOrbitE_im_XYZ[i] = e1ImTmp_XYZ;
					}
				}
			}

			// get Jp(t) for this spatial point
		
			for(int jt=0;jt<nJp;jt++) {

				vector<C3Vec> thisE(nSteps,C3Vec(0,0,0));

				// get Jp(t=jt*dtJp)
				for(int i=0;i<nSteps;i++) {	

					float tTmp = tJp[jt]+thisT[i];
					thisE[i] = thisOrbitE_re_XYZ[i]*cos(wrf*tTmp)-thisOrbitE_im_XYZ[i]*sin(wrf*tTmp);

					//if(tTmp>=-tRF*(nRFCycles-nJpCycles)) { 
					//	thisE[i] = thisOrbitE_re_XYZ[i]*cos(wrf*tTmp)-thisOrbitE_im_XYZ[i]*sin(wrf*tTmp);
					//}
				}

				// This is the hot piece and is done numerically. If there are no kinetic
				// effects here, then this piece should integrate to zero.
				//float trapInt1=0, trapInt2=0, trapInt3=0;
				double qOverm =  thisParticle_XYZ.q/thisParticle_XYZ.m;
				////for(int i=nSteps-2;i>-1;i--) {
				//for(int i=0;i<nSteps-1;i++) {
				//	trapInt1 += qOverm * dtMin/2.0 * (thisE[i].c1+thisE[i+1].c1);
				//	trapInt2 += qOverm * dtMin/2.0 * (thisE[i].c2+thisE[i+1].c2);
				//	trapInt3 += qOverm * dtMin/2.0 * (thisE[i].c3+thisE[i+1].c3);
				//}

				C3Vec thisV1 = qOverm * intC3VecArray ( thisT, thisE );
				//cout << trapInt1 << "   " << thisV1.c1 << endl;

				// This is the cold piece and is done analytically
				//trapInt1 += -qOverm/wrf*(thisOrbitE_re_XYZ[0].c1*cos(wrf*tJp[jt]-_pi/2)-thisOrbitE_im_XYZ[0].c1*sin(wrf*tJp[jt]-_pi/2));

				//C3Vec thisV1(trapInt1,trapInt2,trapInt3);
				float qe = thisParticle_XYZ.q;

				#pragma omp atomic
				j1x[jt] += (particles_XYZ[iP].v_c1+thisV1.c1)*particles_XYZ[iP].weight*qe;
			
			}
		}

#else // END OF LOWMEM CODING ^^^

		vector<CParticle> this_particles_XYZ(particles_XYZ);
		for(int iP=0;iP<particles_XYZ.size();iP++) {
				this_particles_XYZ[iP].c1 = xGrid[iX];
		}
		// Generate linear orbits
		vector<vector<C3Vec> > orbits_XYZ(this_particles_XYZ.size());
		vector<vector<int> > status(this_particles_XYZ.size());

		vector<int> nStepsTaken(this_particles_XYZ.size(),0);
		vector<float> t;

		t.resize(nSteps);

		for(int iP=0;iP<this_particles_XYZ.size();iP++) {

			orbits_XYZ[iP].resize(nSteps);
			status[iP].resize(nSteps);

	 		for(int i=0;i<nSteps;i++) {	

#if DEBUGLEVEL >= 3
					cout << "\tE: " << 
							0.5 * this_particles_XYZ[iP].m * 
							sqrt (pow(this_particles_XYZ[iP].v_c1,2)
											+pow(this_particles_XYZ[iP].v_c2,2)
											+pow(this_particles_XYZ[iP].v_c3,2))/_e << endl;
#endif	
					t[i]=i*dtMin;
					if(this_particles_XYZ[iP].status==0) {
						orbits_XYZ[iP][i] = C3Vec(this_particles_XYZ[iP].c1,this_particles_XYZ[iP].c2,this_particles_XYZ[iP].c3);
						rk4_move ( this_particles_XYZ[iP], dtMin, t[i], b0_CYL, r );
						if(this_particles_XYZ[iP].status==0) {
							status[iP][i] = 0;
							nStepsTaken[iP]++;
						}
						else {
							status[iP][i] = 1;
						}
					}
					else {
						status[iP][i] = 1;
					}
			}
		}

		cout << "\tnSteps: " << nSteps << endl;
		cout << "DONE" << endl;
#if USEPAPI >= 1
		cpuTime0=cpuTime;realTime0=realTime;flpIns0=flpIns;
		papiReturn = PAPI_flops ( &realTime, &cpuTime, &flpIns, &mFlops );
		printf("\nOribit calculation performance ...\n");
		printf("Real_time:\t%f\nProc_time:\t%f\nTotal flpins:\t%lld\nMFLOPS:\t\t%f\n",
					   realTime-realTime0, cpuTime-cpuTime0, flpIns-flpIns0, mFlops);
#endif
		cout << "Interpolating complex E field along trajectories for xGrid " << iX << endl;

		vector<C3Vec> dv(nSteps);	
		vector<vector<C3Vec> >e1(this_particles_XYZ.size());
		vector<vector<vector<C3Vec> > >v1(this_particles_XYZ.size());

		vector<vector<C3Vec> >e1ReHere_XYZ(this_particles_XYZ.size());
		vector<vector<C3Vec> >e1ImHere_XYZ(this_particles_XYZ.size());


		for(int iP=0;iP<this_particles_XYZ.size();iP++) {

			e1ReHere_XYZ[iP].resize(nSteps);
			e1ImHere_XYZ[iP].resize(nSteps);

			for(int i=0;i<nSteps;i++) {

					if(i<=nStepsTaken[iP]&&status[iP][i]==0) {

						stat = 0;
						C3Vec e1ReTmp_XYZ = kj_interp1D ( orbits_XYZ[iP][i].c1, r, e1Re_XYZ, stat );
						stat = 0;
						C3Vec e1ImTmp_XYZ = kj_interp1D ( orbits_XYZ[iP][i].c1, r, e1Im_XYZ, stat );

						e1ReHere_XYZ[iP][i] = e1ReTmp_XYZ;
						e1ImHere_XYZ[iP][i] = e1ImTmp_XYZ;
					}
					else {
						//printf("\t%s line: %i\n",__FILE__,__LINE__);
						//cout<<"i < nStepsTaken[iP]"<<endl;
						e1ReHere_XYZ[iP][i] = C3Vec(0,0,0);
						e1ImHere_XYZ[iP][i] = C3Vec(0,0,0);
				}

			}

		}
		
		cout << "DONE" << endl;
#if USEPAPI >= 1
		cpuTime0=cpuTime;realTime0=realTime;flpIns0=flpIns;
		papiReturn = PAPI_flops ( &realTime, &cpuTime, &flpIns, &mFlops );
		printf("\nOribit E interpolation performance ...\n");
		printf("Real_time:\t%f\nProc_time:\t%f\nTotal flpins:\t%lld\nMFLOPS:\t\t%f\n",
					   realTime-realTime0, cpuTime-cpuTime0, flpIns-flpIns0, mFlops);
#endif
		// Calculate jP1 for each time at the spatial point


#if DEBUGLEVEL >= 1
		cout << "\tnThermal: " << nThermal << endl;
		cout << "\tvTh: " << vTh << endl;
#endif

		// Check density using non-grid method
		float densityCheck = 0;
		for(int iP=0;iP<particles_XYZ_0.size();iP++){
				densityCheck += particles_XYZ_0[iP].weight;
		}

		cout << "Density on f0 using non-grid method: " << densityCheck << endl;

		cout << "DONE" << endl;


		for(int iP=0;iP<this_particles_XYZ.size();iP++) {

				e1[iP].resize(nSteps);
				v1[iP].resize(nJp);

				for(int iJ=0;iJ<nJp;iJ++) {

					v1[iP][iJ].resize(nSteps);

				}
		}
#if USEPAPI >= 1
		float eT_cpuTime=0, eT_realTime=0, eT_flpIns=0, eT_mFlops=0;
		float vT_cpuTime=0, vT_realTime=0, vT_flpIns=0, vT_mFlops=0;
#endif
		for(int jt=0;jt<nJp;jt++) {

			//tJp[jt] = jt*dtJp;
			//cout << "Create f1 for this tJp: " << tJp[jt] << endl;
			

			// Get e1 magnitude along orbit
			for(int iP=0;iP<this_particles_XYZ.size();iP++) {


				for(int i=0;i<nSteps;i++) {
					v1[iP][jt][i] = C3Vec(0,0,0);
				}
#if USEPAPI >= 1
				cpuTime0=cpuTime;realTime0=realTime;flpIns0=flpIns;
				papiReturn = PAPI_flops ( &realTime, &cpuTime, &flpIns, &mFlops );
#endif	
				for(int i=0;i<nSteps;i++) {	

					float tTmp = tJp[jt]+t[i];
					if(tTmp>=-tRF*(nRFCycles-nJpCycles)) { //i<=nStepsTaken[iP]) { 

						// Get E(t) along orbit 
						e1[iP][i] = e1ReHere_XYZ[iP][i]*cos(wrf*tTmp)
								+e1ImHere_XYZ[iP][i]*sin(wrf*tTmp);
						
					}
					else {
						//printf("\t%s line: %i\n",__FILE__,__LINE__);
						//cout << "i > nStepsTaken" << endl;
						//exit (1);
						e1[iP][i] = C3Vec(0,0,0);
					}	
				}
#if USEPAPI >= 1
				cpuTime0=cpuTime;realTime0=realTime;flpIns0=flpIns;
				papiReturn = PAPI_flops ( &realTime, &cpuTime, &flpIns, &mFlops );

				eT_realTime += (realTime-realTime0);
				eT_cpuTime += (cpuTime-cpuTime0);
				eT_flpIns += (flpIns-flpIns0);
				eT_mFlops += mFlops;
#endif
				// Intergrate e1 from t=-inf to 0 to get v1
				v1[iP][jt][nSteps-1].c1=0;v1[iP][jt][nSteps-1].c2=0;v1[iP][jt][nSteps-1].c3=0;

				float trapInt1=0, trapInt2=0, trapInt3=0;
				double qOverm =  this_particles_XYZ[iP].q/this_particles_XYZ[iP].m;
				for(int i=nSteps-2;i>-1;i--) {

					trapInt1 += qOverm * dtMin/2.0 * (e1[iP][i].c1+e1[iP][i+1].c1);
					trapInt2 += qOverm * dtMin/2.0 * (e1[iP][i].c2+e1[iP][i+1].c2);
					trapInt3 += qOverm * dtMin/2.0 * (e1[iP][i].c3+e1[iP][i+1].c3);
					//cout << "dtMin: " << dtMin << "  t[i+1]-t[i]: " << (t[i+1]-t[i]) << endl;

					v1[iP][jt][i].c1 = trapInt1;
					v1[iP][jt][i].c2 = trapInt2;
					v1[iP][jt][i].c3 = trapInt3;

				}
#if USEPAPI >= 1
				cpuTime0=cpuTime;realTime0=realTime;flpIns0=flpIns;
				papiReturn = PAPI_flops ( &realTime, &cpuTime, &flpIns, &mFlops );

				vT_realTime += (realTime-realTime0);
				vT_cpuTime += (cpuTime-cpuTime0);
				vT_flpIns += (flpIns-flpIns0);
				vT_mFlops += mFlops;
#endif
			}
#if USEPAPI >= 1
			cpuTime0=cpuTime;realTime0=realTime;flpIns0=flpIns;
			papiReturn = PAPI_flops ( &realTime, &cpuTime, &flpIns, &mFlops );
#endif
			float qe = this_particles_XYZ[0].q;

			j1x[jt] = 0;
			for(int iP=0;iP<this_particles_XYZ.size();iP++) {
					j1x[jt] += (particles_XYZ_0[iP].v_c1+v1[iP][jt][0].c1)*particles_XYZ_0[iP].weight;
					//j1x[jt] -= (particles_XYZ_0[iP].v_c1)*particles_XYZ_0[iP].weight;
			}
			j1x[jt] = j1x[jt] * qe;
			//cout << "j1x["<<jt<<"]: "<< j1x[jt]<<endl;
			
		}

#endif // END OF HIGHMEM CODING ^^^


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
		
				NcVar nc_e1_x = ncOrbitsFile.addVar("e1_x",ncFloat,nc_nPxnSteps);
				NcVar nc_e1_y = ncOrbitsFile.addVar("e1_y",ncFloat,nc_nPxnSteps);
				NcVar nc_e1_z = ncOrbitsFile.addVar("e1_z",ncFloat,nc_nPxnSteps);
		
				NcVar nc_v1_x = ncOrbitsFile.addVar("v1x",ncFloat,nc_nPxnJpxnSteps);
				NcVar nc_v1_y = ncOrbitsFile.addVar("v1y",ncFloat,nc_nPxnJpxnSteps);
				NcVar nc_v1_z = ncOrbitsFile.addVar("v1z",ncFloat,nc_nPxnJpxnSteps);
		
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

						for(int iS=0;iS<nSteps;iS++){tmpData[iS] = e1[iP][iS].c1;}
						nc_e1_x.putVar(startpA,countpA,&tmpData[0]);
						for(int iS=0;iS<nSteps;iS++){tmpData[iS] = e1[iP][iS].c2;}
						nc_e1_y.putVar(startpA,countpA,&tmpData[0]);
						for(int iS=0;iS<nSteps;iS++){tmpData[iS] = e1[iP][iS].c3;}
						nc_e1_z.putVar(startpA,countpA,&tmpData[0]);

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

						}
				}
		
				vector<size_t> startp (1,0);
				vector<size_t> countp (1,nSteps);
		
				nc_t.putVar(startp,countp,&t[0]);
		
		
		}
				catch(exceptions::NcException &e) {
						cout << "NetCDF: unknown error" << endl;
						e.what();
		}

		cout << "DONE" << endl;
#endif

		// Write current to file
	
		//cout << "Writing jP to file ... ";

		stringstream ncjPFileName;
		ncjPFileName << "output/";
		ncjPFileName << runIdent.c_str();
		// check directory exists
		struct stat st;
		if(stat(ncjPFileName.str().c_str(),&st) != 1) {
			int mkDirStat = mkdir(ncjPFileName.str().c_str(),S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
		}
		ncjPFileName << "/jP_";
		ncjPFileName << setw(3) << setfill('0') << iX;
	   	ncjPFileName << ".nc"; 	
		NcFile ncjPFile (ncjPFileName.str().c_str(), NcFile::replace);

		NcDim nc_nJp = ncjPFile.addDim("nJp", nJp);
		NcDim nc_scalar = ncjPFile.addDim("scalar", 1);

		NcVar nc_t = ncjPFile.addVar("t",ncFloat,nc_nJp);

		NcVar nc_x = ncjPFile.addVar("x",ncFloat,nc_scalar);
		NcVar nc_freq = ncjPFile.addVar("freq",ncFloat,nc_scalar);

		NcVar nc_j1x = ncjPFile.addVar("j1x",ncFloat,nc_nJp);
		NcVar nc_j1y = ncjPFile.addVar("j1y",ncFloat,nc_nJp);
		NcVar nc_j1z = ncjPFile.addVar("j1z",ncFloat,nc_nJp);

		nc_x.putVar(&xGrid[iX]);
		nc_freq.putVar(&freq);

		vector<size_t> startp (1,0);
		vector<size_t> countp (1,nJp);

		nc_t.putVar(startp,countp,&tJp[0]);

		nc_j1x.putVar(startp,countp,&j1x[0]);
		nc_j1y.putVar(startp,countp,&j1y[0]);
		nc_j1z.putVar(startp,countp,&j1z[0]);


	} // End of xGrid loop

	//ProfilerStop();

	cout << "DONE" << endl;

	return EXIT_SUCCESS;
}
