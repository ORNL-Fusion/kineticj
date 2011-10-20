#include <string>
#include <iostream>
#include <iomanip>
#include <cstdlib>
#include <netcdf>
#include <vector>
#include <algorithm>
#include <complex>
#include "constants.hpp"
//#include <google/profiler.h>

#ifdef __CUDA_ARCH__
#define PRINT cuPrintf 
#else
#define PRINT printf
#endif

class CSpecies {
		public:
				double m, q;
				double amu;
				int Z;
				std::string name;

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
		name = std::string(_s);
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
				std::complex<float> c1, c2, c3;

				C3VecI () {c1=std::complex<float>(0.0f,0.0f);c2=std::complex<float>(0.0f,0.0f);c3=std::complex<float>(0.0f,0.0f);};
				C3VecI ( std::complex<float> _c1, std::complex<float> _c2, std::complex<float> _c3 ) {c1=_c1;c2=_c2;c3=_c3;};
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
		return C3Vec(*this)+=other;
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

// First-order orbits
C3Vec rk4_evalf ( CParticle &p, const float &t, const C3Vec &v, const C3Vec &x,
				const std::vector<C3Vec> &b0Vec, const std::vector<C3VecI> &e1, const float wrf ) {

	C3Vec b0(0,0,0), F;

	C3Vec v_x_b0 ( v.c2*b0.c3-v.c3*b0.c2, -1.0*(v.c1*b0.c3-v.c3*b0.c1), v.c1*b0.c2-v.c2*b0.c1); 
	//C3Vec F ( std::real(e1.c1) * cos ( wrf * t ) + std::imag(e1.c1) * sin ( wrf * t ) + v_x_b0.c1,
	//	  	std::real(e1.c2) * cos ( wrf * t ) + std::imag(e1.c2) * sin ( wrf * t ) + v_x_b0.c2,
	//	  	std::real(e1.c3) * cos ( wrf * t ) + std::imag(e1.c3) * sin ( wrf * t ) + v_x_b0.c3 );

	return F*(p.q/p.m);	
}

C3Vec kj_interp1D ( const float &x, const std::vector<float> &xVec, const std::vector<C3Vec> &yVec ) {

	if(x<xVec.front()||x>xVec.back()) {
			printf("\t%s line: %i\n",__FILE__,__LINE__);
			std::cout << "\tERROR: Interpolation pt off grid." << std::endl;
			std::cout << "\tx: " << x << " x.front(): " << xVec.front() << " x.back(): " << xVec.back() << std::endl;
			exit (1);
	}

	float _x = (x-xVec.front())/(xVec.back()-xVec.front())*(xVec.size()-1);
	float x0 = floor(_x);
	float x1 = ceil(_x);

	// Catch for particle at point
	if(x0==x1) {
		//std::cout << "x0: " << x0 << " x1: " <<x1<< " _x: "<<_x << std::endl;
		//std::cout << "Particle at point catch: " << x0/x1 << "  "  << abs(1.0-x0/x1) << std::endl;
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
				const C3Vec &v_XYZ, const C3Vec &x, const std::vector<C3Vec> &b0Vec_CYL,
			  	const std::vector<float> &rVec ) {

	// Interpolate b0 at location in CYL
	
	float _r = sqrt ( pow(x.c1,2) + pow(x.c2,2) );
	float _p = atan2 ( x.c2, x.c1 );

#if DEBUGLEVEL >= 3
	std::cout << "\t\t\tx: " << x.c1 << " y: " << x.c2 << " z: " << x.c3 << std::endl;
	std::cout << "\t\t\tr: " << _r << " p: " << _p << std::endl;
	std::cout << "\t\t\trVec.front(): " << rVec.front() << std::endl;
	std::cout << "\t\t\tv_XYZ: " << v_XYZ.c1 << "  " << v_XYZ.c2 << "  " << v_XYZ.c3 << std::endl;
#endif

	C3Vec b0_CYL, b0_XYZ;

	b0_CYL = kj_interp1D ( _r, rVec, b0Vec_CYL );

	b0_XYZ = C3Vec( cos(_p)*b0_CYL.c1-sin(_p)*b0_CYL.c2+0,
					sin(_p)*b0_CYL.c1+cos(_p)*b0_CYL.c2+0,
					0+0+1*b0_CYL.c3 );

	C3Vec v_x_b0 ( v_XYZ.c2*b0_XYZ.c3-v_XYZ.c3*b0_XYZ.c2, 
					-1.0*(v_XYZ.c1*b0_XYZ.c3-v_XYZ.c3*b0_XYZ.c1), 
					v_XYZ.c1*b0_XYZ.c2-v_XYZ.c2*b0_XYZ.c1);

#if DEBUGLEVEL >= 3
	std::cout << "\tvxb0: " << v_x_b0.c1 << "  " << v_x_b0.c2 << "  " << v_x_b0.c3 << std::endl;
	std::cout << "\tp.q/p.m: " << p.q/p.m << std::endl;
#endif

	return v_x_b0*(p.q/p.m);	
}

// Zero-order orbits
void rk4_move ( CParticle &p, const float &dt, const float &t0, 
				const std::vector<C3Vec> &b0, const std::vector<float> &r ) {

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
		std::cout << "\tx0_XYZ: " << xn0.c1 << "  " << xn0.c2 << "  " << xn0.c3 << std::endl;
		std::cout << "\tv0_XYZ: " << yn0.c1 << "  " << yn0.c2 << "  " << yn0.c3 << std::endl;
		std::cout << "\tx1_XYZ: " << xn1.c1 << "  " << xn1.c2 << "  " << xn1.c3 << std::endl;
		std::cout << "\tv1_XYZ: " << yn1.c1 << "  " << yn1.c2 << "  " << yn1.c3 << std::endl;
		std::cout << "\tE: " << 0.5 * p.m * sqrt (pow(p.v_c1,2)+pow(p.v_c2,2)+pow(p.v_c3,2))/_e << std::endl;
#endif

}

// First-order orbits
void rk4_move ( CParticle &p, float dt, float t0, 
				const std::vector<C3Vec> &b0, const std::vector<C3VecI> &e1, const float wrf ) {

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


float maxC3VecAbs ( const std::vector<C3Vec> &input ) {

	std::vector<float> inputAbs(input.size());
	for(int i=0;i<input.size();i++) {
		inputAbs[i] = sqrt(pow(input[i].c1,2)+pow(input[i].c2,2)+pow(input[i].c3,2));
	}
	return *max_element(inputAbs.begin(),inputAbs.end());
}


// Calculate the jP given some know E and f(v)

int main ( int argc, char **argv )
{
		// Read E
	
		std::string rsfwc_fName ( "data/rsfwc_1d.nc" );	
		std::cout << "Reading rsfwc data file" << rsfwc_fName << std::endl;

		// Here we are using the cxx-4 netcdf interface by Lynton Appel
		// This needs netCDF 4.1.1 or later build with
		// ./configure --enable-cxx-4 [plus other options]

		std::vector<float> r, b0_r, b0_p, b0_z,
				e_r_re, e_p_re, e_z_re,
				e_r_im, e_p_im, e_z_im;
		std::vector<C3Vec> b0_CYL, b0_XYZ;
		
		float wrf;

		std::vector<std::complex<float> > e_r, e_p, e_z;	

		try {
				netCDF::NcFile dataFile ( rsfwc_fName.c_str(), netCDF::NcFile::read );
	
				netCDF::NcDim nc_nR(dataFile.getDim("nR"));
				netCDF::NcDim nc_scalar(dataFile.getDim("scalar"));
	
				int nR = nc_nR.getSize();
	
				std::cout << "\tnR: " << nR << std::endl;
	
				netCDF::NcVar nc_r(dataFile.getVar("r"));
				netCDF::NcVar nc_wrf(dataFile.getVar("wrf"));

				netCDF::NcVar nc_b0_r(dataFile.getVar("B0_r"));
				netCDF::NcVar nc_b0_p(dataFile.getVar("B0_p"));
				netCDF::NcVar nc_b0_z(dataFile.getVar("B0_z"));

				netCDF::NcVar nc_e_r_re(dataFile.getVar("e_r_re"));
				netCDF::NcVar nc_e_p_re(dataFile.getVar("e_p_re"));
				netCDF::NcVar nc_e_z_re(dataFile.getVar("e_z_re"));
				netCDF::NcVar nc_e_r_im(dataFile.getVar("e_r_im"));
				netCDF::NcVar nc_e_p_im(dataFile.getVar("e_p_im"));
				netCDF::NcVar nc_e_z_im(dataFile.getVar("e_z_im"));

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
				nc_wrf.getVar(&wrf);

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
						e_r.push_back(std::complex<float>( e_r_re[i], e_r_im[i] ) );
						e_p.push_back(std::complex<float>( e_p_re[i], e_p_im[i] ) );
						e_z.push_back(std::complex<float>( e_z_re[i], e_z_im[i] ) );
				}

				std::cout << "\tR[0]: " << r[0] << ", R["<<nR<<"]: " << r[r.size()-1] << std::endl;
				std::cout << "\twrf: " << wrf << std::endl;
				std::vector<float>::iterator min = std::min_element(b0_p.begin(),b0_p.end());
				std::vector<float>::iterator max = std::max_element(b0_p.begin(),b0_p.end());
				std::cout << "\tmin(b0_p): " << *min << std::endl;
				std::cout << "\tmax(b0_p): " << *max << std::endl;
				std::cout << "\tabs(e_r[nR/2]): " << std::abs(e_r[nR/2]) << std::endl;
				std::cout << "\tabs(e_p[nR/2]): " << std::abs(e_p[nR/2]) << std::endl;
				std::cout << "\tabs(e_z[nR/2]): " << std::abs(e_z[nR/2]) << std::endl;
		}
		catch(netCDF::exceptions::NcException &e) {
				std::cout << "NetCDF: unknown error" << std::endl;
				e.what();
		}

		// Create f0(v)

		std::string particleList_fName ( "data/f_1keV_electrons.nc" );	
		std::cout << "Reading particle list " << particleList_fName << std::endl;

		std::vector<float> p_x, p_y, p_z, p_vx, p_vy, p_vz, p_amu, p_weight;
		std::vector<int> p_Z;
		float vTh;
		int nThermal;
		
		try {
				netCDF::NcFile dataFile ( particleList_fName.c_str(), netCDF::NcFile::read );
	
				netCDF::NcDim nc_nP(dataFile.getDim("nP"));
	
				int nP = nc_nP.getSize();
	
				std::cout << "\tnP: " << nP << std::endl;

				netCDF::NcVar nc_p_amu(dataFile.getVar("amu"));
				netCDF::NcVar nc_p_Z(dataFile.getVar("Z"));

				netCDF::NcVar nc_p_x(dataFile.getVar("x"));
				netCDF::NcVar nc_p_y(dataFile.getVar("y"));
				netCDF::NcVar nc_p_z(dataFile.getVar("z"));
				
				netCDF::NcVar nc_p_vx(dataFile.getVar("vx"));
				netCDF::NcVar nc_p_vy(dataFile.getVar("vy"));
				netCDF::NcVar nc_p_vz(dataFile.getVar("vz"));

				netCDF::NcVar nc_p_weight(dataFile.getVar("weight"));

				netCDF::NcVar nc_nThermal(dataFile.getVar("nThermal"));
				netCDF::NcVar nc_vTh(dataFile.getVar("vTh"));

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
		catch(netCDF::exceptions::NcException &e) {
				std::cout << "NetCDF: unknown error" << std::endl;
				e.what();
		}

		std::vector<CParticle> particles_XYZ;
		particles_XYZ.resize(p_x.size());

		for(int i=0;i<particles_XYZ.size();i++){

				CParticle thisParticle (p_x[i],p_y[i],p_z[i],
								p_vx[i],p_vy[i],p_vz[i],
								p_amu[i],p_Z[i],p_weight[i]);
				particles_XYZ[i] = thisParticle;
				particles_XYZ[i].number = i;

				//std::cout << "\tamu: " << p_amu[i] << std::endl;

				//std::cout << "\tParticle[" << i << "]: " 
				//		<< particles_XYZ[i].c1 << "  " 
				//		<< particles_XYZ[i].c2 << "  " 
				//		<< particles_XYZ[i].c3 << "  " 
				//		<< particles_XYZ[i].v_c1 << "  " 
				//		<< particles_XYZ[i].v_c2 << "  " 
				//		<< particles_XYZ[i].v_c3 << "  " 
				//		<< particles_XYZ[i].q << "  " 
				//		<< particles_XYZ[i].m << "  " 
				//		<< std::endl;
		}

	float xGridMin = 9.60;
	float xGridMax = 10.4;
	float xGridRng = xGridMax-xGridMin;
	int nXGrid = 40;
	float xGridStep = xGridRng/(nXGrid-1);
	std::vector<float> xGrid(nXGrid);

	for(int iX=0;iX<nXGrid;iX++) {
			xGrid[iX] = xGridMin+iX*xGridStep;
			std::cout << "\t\txGrid[iX]: " << xGrid[iX] << std::endl;
	}

	std::vector<CParticle> particles_XYZ_0(particles_XYZ);

	//std::string googlePerfFileName = "/home/dg6/code/kineticj/googlep";
	//ProfilerStart(googlePerfFileName.c_str());

	for(int iX=0;iX<nXGrid;iX++) {

		std::vector<CParticle> particles_XYZ_thisX(particles_XYZ);

		for(int iP=0;iP<particles_XYZ_thisX.size();iP++) {
			particles_XYZ_thisX[iP].c1 = xGrid[iX];
		}

		// Generate linear orbits

		std::cout << "Calculating orbit for xGrid " << iX << std::endl;

		int nRFCycles = 4;
		int nStepsPerCycle = 100.0;

		float tRF = (2*_pi)/wrf;
		float dtMin = -tRF/nStepsPerCycle;
		//float tEnd = dtMin*nStepsPerCycle*nRFCycles;

		std::vector<std::vector<C3Vec> > orbit(particles_XYZ_thisX.size());
		std::vector<int> nStepsTaken(particles_XYZ_thisX.size(),0);
		std::vector<float> t;

		int nSteps = nRFCycles*nStepsPerCycle+1;
		//std::cout.precision(15);
		//std::cout << std::fixed << nSteps*dtMin/tRF << std::endl;
		//exit(1);
		t.resize(nSteps);

		for(int iP=0;iP<particles_XYZ_thisX.size();iP++) {
		//for(int iP=0;iP<1;iP++) {

			orbit[iP].resize(nSteps);

	 		for(int i=0;i<nSteps;i++) {	

#if DEBUGLEVEL >= 3
					std::cout << "\tE: " << 
							0.5 * particles_XYZ_thisX[iP].m * 
							sqrt (pow(particles_XYZ_thisX[iP].v_c1,2)
											+pow(particles_XYZ_thisX[iP].v_c2,2)
											+pow(particles_XYZ_thisX[iP].v_c3,2))/_e << std::endl;
#endif	
					t[i]=i*dtMin;
					if(particles_XYZ_thisX[iP].status==0) {
						orbit[iP][i] = C3Vec(particles_XYZ_thisX[iP].c1,particles_XYZ_thisX[iP].c2,particles_XYZ_thisX[iP].c3);
						nStepsTaken[iP]++;
						rk4_move ( particles_XYZ_thisX[iP], dtMin, t[i], b0_CYL, r );
					}
					else {
						printf("\t%s line: %i\n",__FILE__,__LINE__);
						std::cout << "Status != 0" << std::endl;
						exit(1);
					}
			}
		}

		//std::cout.precision(15);
		//std::cout << std::fixed << t[nSteps-1]/tRF << std::endl;
		//exit(1);

		std::cout << "\tnSteps: " << nSteps << std::endl;
		std::cout << "DONE" << std::endl;

		std::cout << "Interpolating complex E field along trajectories for xGrid " << iX << std::endl;

		std::vector<C3Vec> dv(nSteps);	
		std::vector<std::vector<C3Vec> >e1(particles_XYZ_thisX.size());
		std::vector<std::vector<std::vector<C3Vec> > >v1(particles_XYZ_thisX.size());

		std::vector<std::vector<C3Vec> >e1ReHere_CYL(particles_XYZ_thisX.size());
		std::vector<std::vector<C3Vec> >e1ImHere_CYL(particles_XYZ_thisX.size());

		std::vector<C3Vec> e1Re_CYL, e1Im_CYL;
		e1Re_CYL.resize(e_r.size());
		e1Im_CYL.resize(e_r.size());

		for(int i=0;i<e_r.size();i++) {
			e1Re_CYL[i].c1 = std::real(e_r[i]);
			e1Re_CYL[i].c2 = std::real(e_p[i]);
			e1Re_CYL[i].c3 = std::real(e_z[i]);
			e1Im_CYL[i].c1 = std::imag(e_r[i]);
			e1Im_CYL[i].c2 = std::imag(e_p[i]);
			e1Im_CYL[i].c3 = std::imag(e_z[i]);
		}

		for(int iP=0;iP<particles_XYZ_thisX.size();iP++) {

			e1ReHere_CYL[iP].resize(nSteps);
			e1ImHere_CYL[iP].resize(nSteps);

			for(int i=0;i<nSteps;i++) {

					if(i<=nStepsTaken[iP]) {

						// Interpolate e1Now to here, done in CYL
					
						float _r = sqrt ( pow(orbit[iP][i].c1,2) + pow(orbit[iP][i].c2,2) );
						float _p = atan2 ( orbit[iP][i].c2, orbit[iP][i].c1 );

						C3Vec e1ReTmp_CYL = kj_interp1D ( _r, r, e1Re_CYL );
						C3Vec e1ImTmp_CYL = kj_interp1D ( _r, r, e1Im_CYL );

						e1ReHere_CYL[iP][i] = e1ReTmp_CYL;
						e1ImHere_CYL[iP][i] = e1ImTmp_CYL;

						//std::cout 	<< "e1ReTmp.c1: "<<e1ReTmp_CYL.c1
						//			<<" e1ReTmp.c2: "<<e1ReTmp_CYL.c2
						//			<<" e1ReTmp.c3: "<<e1ReTmp_CYL.c3<< std::endl;
	
					}
					else {
						printf("\t%s line: %i\n",__FILE__,__LINE__);
						std::cout<<"i < nStepsTaken[iP]"<<std::endl;
						exit (1);
						e1ReHere_CYL[iP][i] = C3Vec(0,0,0);
						e1ImHere_CYL[iP][i] = C3Vec(0,0,0);
					}

			}

		}
		//exit(1);
		std::cout << "DONE" << std::endl;

		// Calculate jP1 for each time at the spatial point

		int nJpCycles = 2;
		int nJpPerCycle = 20;
		int nJp = nJpCycles * nJpPerCycle + 1;
		float dtJp = tRF / nJpPerCycle;
		std::vector<float> tJp(nJp,0);

		std::vector<std::vector<std::vector<float> > > f_XYZ, f_XYZ_0;
		std::vector<float> vxGrid, vyGrid, vzGrid;

#if DEBUGLEVEL >= 1
		std::cout << "\tnThermal: " << nThermal << std::endl;
		std::cout << "\tvTh: " << vTh << std::endl;
#endif

		//int nx=1000, ny=20, nz=20;
		//vxGrid.resize(nx);vyGrid.resize(ny);vzGrid.resize(nz);

		//float vxMin = -nThermal*vTh*7000;
		//float vxMax = -vxMin;
		//float vxRange = (vxMax-vxMin);
		//float dVx = vxRange / (vxGrid.size()-1);

		//float vyMin = -nThermal*vTh*5;
		//float vyMax = -vyMin;
		//float vyRange = (vyMax-vyMin);
		//float dVy = vyRange / (vyGrid.size()-1);

		//float vzMin = -nThermal*vTh*5;
		//float vzMax = -vzMin;
		//float vzRange = (vzMax-vzMin);
		//float dVz = vzRange / (vzGrid.size()-1);
		//
		//float dV = dVx * dVy * dVz;

		//std::cout << "\tdVx: " << dVx << std::endl;
		//for(int i=0;i<vxGrid.size();i++) {
		//		vxGrid[i] = vxMin + i*dVx;
		//}
		//for(int j=0;j<vyGrid.size();j++) {
		//		vyGrid[j] = vyMin + j*dVy;
		//}
		//for(int k=0;k<vzGrid.size();k++) {
		//		vzGrid[k] = vzMin + k*dVz;
		//}

		//f_XYZ.resize(nx);
		//f_XYZ_0.resize(nx);
		//for(int i=0;i<nx;i++) {
		//		f_XYZ[i].resize(ny);
		//		f_XYZ_0[i].resize(ny);
		//		for(int j=0;j<ny;j++) {
		//				f_XYZ[i][j].resize(nz);
		//				f_XYZ_0[i][j].resize(nz);
		//		}
		//}

		//std::cout << "Create f0 ..." << std::endl;

		//// Create the initial f
		//for(int iP=0;iP<particles_XYZ_0.size();iP++) {
		//		float iix = (particles_XYZ_0[iP].v_c1-vxMin)/vxRange*(vxGrid.size()-1);
		//		if(iix<0 || iix>(nx-1)){
		//				std::cout<<"Outside v grid: "<<particles_XYZ_0[iP].v_c1<<std::endl;
		//				std::cout<<"max v: "<<vxMax<<std::endl;
		//		}
		//		float iiy = (particles_XYZ_0[iP].v_c2-vyMin)/vyRange*(vyGrid.size()-1);
		//		if(iiy<0 || iiy>(ny-1)){
		//				std::cout<<"Outside v grid: "<<particles_XYZ_0[iP].v_c2<<std::endl;
		//				std::cout<<"max v: "<<vyMax<<std::endl;
		//		}
		//		float iiz = (particles_XYZ_0[iP].v_c3-vzMin)/vzRange*(vzGrid.size()-1);
		//		f_XYZ_0[iix][iiy][iiz] += particles_XYZ_0[iP].weight/dV;
		//}	

		//// Check density against expected
		//float densityCheck = 0;
		//for(int i=0;i<nx;i++){
		//		for(int j=0;j<ny;j++){
		//				for(int k=0;k<nz;k++){
		//					densityCheck += f_XYZ_0[i][j][k]*dV;
		//				}
		//		}
		//}
		//std::cout << "Density on f0: " << densityCheck << std::endl;

		// Check density using non-grid method
		float densityCheck = 0;
		for(int iP=0;iP<particles_XYZ_0.size();iP++){
				densityCheck += particles_XYZ_0[iP].weight;
		}

		std::cout << "Density on f0 using non-grid method: " << densityCheck << std::endl;

		std::cout << "DONE" << std::endl;

		std::vector<float> j1x(nJp,0), j1y(nJp,0), j1z(nJp,0), tJ(nJp,0);

		for(int iP=0;iP<particles_XYZ_thisX.size();iP++) {

				e1[iP].resize(nSteps);
				v1[iP].resize(nJp);

				for(int iJ=0;iJ<nJp;iJ++) {

					v1[iP][iJ].resize(nSteps);

				}
		}

		for(int jt=0;jt<nJp;jt++) {

			tJp[jt] = jt*dtJp;
			std::cout << "Create f1 for this tJp: " << tJp[jt] << std::endl;

			// Get e1 magnitude along orbit
			for(int iP=0;iP<particles_XYZ_thisX.size();iP++) {

				for(int i=0;i<nSteps;i++) {
					v1[iP][jt][i] = C3Vec(0,0,0);
				}

				for(int i=0;i<nSteps;i++) {	

					float tTmp = tJp[jt]+t[i];
					if(tTmp>=-tRF*(nRFCycles-nJpCycles)) { //i<=nStepsTaken[iP]) { 

						// Get E(t) along orbit 
						C3Vec e1NowAndHere_CYL;

						float _r = sqrt ( pow(orbit[iP][i].c1,2) + pow(orbit[iP][i].c2,2) );
						float _p = atan2 ( orbit[iP][i].c2, orbit[iP][i].c1 );

						e1NowAndHere_CYL = e1ReHere_CYL[iP][i]*cos(wrf*tTmp-_pi/2)+e1ImHere_CYL[iP][i]*sin(wrf*tTmp-_pi/2); 

						C3Vec e1NowAndHere_XYZ;

						e1NowAndHere_XYZ = C3Vec( 
										cos(_p)*e1NowAndHere_CYL.c1-sin(_p)*e1NowAndHere_CYL.c2+0,
										sin(_p)*e1NowAndHere_CYL.c1+cos(_p)*e1NowAndHere_CYL.c2+0,
										0+0+1*e1NowAndHere_CYL.c3 );

						e1[iP][i] = e1NowAndHere_XYZ;

					}
					else {
						//printf("\t%s line: %i\n",__FILE__,__LINE__);
						//std::cout << "i > nStepsTaken" << std::endl;
						//exit (1);
						e1[iP][i] = C3Vec(0,0,0);
					}
				}
	
				// Intergrate e1 from t=-inf to 0 to get v1
				v1[iP][jt][nSteps-1].c1=0;v1[iP][jt][nSteps-1].c2=0;v1[iP][jt][nSteps-1].c3=0;

				float trapInt1=0, trapInt2=0, trapInt3=0;
				double qOverm =  particles_XYZ_thisX[iP].q/particles_XYZ_thisX[iP].m;
				for(int i=nSteps-2;i>-1;i--) {

					trapInt1 += qOverm * dtMin/2.0 * (e1[iP][i].c1+e1[iP][i+1].c1);
					trapInt2 += qOverm * dtMin/2.0 * (e1[iP][i].c2+e1[iP][i+1].c2);
					trapInt3 += qOverm * dtMin/2.0 * (e1[iP][i].c3+e1[iP][i+1].c3);
					//std::cout << "dtMin: " << dtMin << "  t[i+1]-t[i]: " << (t[i+1]-t[i]) << std::endl;

					v1[iP][jt][i].c1 = trapInt1;
					v1[iP][jt][i].c2 = trapInt2;
					v1[iP][jt][i].c3 = trapInt3;

				}
			}


			//for(int i=0;i<nx;i++){
			//		for(int j=0;j<ny;j++){
			//				for(int k=0;k<nz;k++){
			//					f_XYZ[i][j][k]=0;
			//				}
			//		}
			//}

			//for(int iP=0;iP<particles_XYZ_thisX.size();iP++) {
			//		float iix = (particles_XYZ_0[iP].v_c1+v1[iP][jt][0].c1-vxMin)/vxRange*(vxGrid.size()-1);
			//		if(iix<0 || iix>(nx-1)){
			//				std::cout<<"\t\tError - v: "<<particles_XYZ_thisX[iP].v_c1<<std::endl;
			//				std::cout<<"\t\tError - max v: "<<vxMax<<std::endl;
			//				std::cout<<"\t\tError - v+v1: "<<particles_XYZ_0[iP].v_c1+v1[iP][jt][0].c1-vxMin<<std::endl;
			//				std::cout<<"\t\tError - iix: "<<iix<<std::endl;
			//				exit(1);
			//		}
			//		float iiy = (particles_XYZ_0[iP].v_c2+v1[iP][jt][0].c2-vyMin)/vyRange*(vyGrid.size()-1);
			//		if(iiy<0 || iiy>(ny-1)){
			//				std::cout<<"Outside v grid: "<<particles_XYZ_thisX[iP].v_c2<<std::endl;
			//				std::cout<<"max v: "<<vyMax<<std::endl;
			//				exit(1);
			//		}
			//		float iiz = (particles_XYZ_0[iP].v_c3+v1[iP][jt][0].c3-vzMin)/vzRange*(vzGrid.size()-1);
			//		if(iiz<0 || iiy>(nz-1)){
			//				std::cout<<"Outside v grid: "<<particles_XYZ_thisX[iP].v_c2<<std::endl;
			//				std::cout<<"max v: "<<vzMax<<std::endl;
			//				exit(1);
			//		}
			//		f_XYZ[iix][iiy][iiz] += particles_XYZ_thisX[iP].weight/dV;
			//}	

			//j1x[jt] = 0;
			//j1y[jt] = 0;
			//j1z[jt] = 0;

			float qe =  particles_XYZ_thisX[0].q;
			//densityCheck = 0;
			//for(int i=0;i<nx;i++){
			//		for(int j=0;j<ny;j++){
			//				for(int k=0;k<nz;k++){
			//					j1x[jt] += qe*vxGrid[i]*(f_XYZ[i][j][k]-f_XYZ_0[i][j][k])*dV;
			//					j1y[jt] += qe*vyGrid[j]*(f_XYZ[i][j][k]-f_XYZ_0[i][j][k])*dV;
			//					j1z[jt] += qe*vzGrid[k]*(f_XYZ[i][j][k]-f_XYZ_0[i][j][k])*dV;
			//					densityCheck += f_XYZ[i][j][k]*dV;
			//				}
			//		}
			//}

			//std::cout << "Density on f1: " << densityCheck << std::endl;
			//std::cout << "j1x["<<jt<<"]: "<< j1x[jt]<<std::endl;

			j1x[jt] = 0;
			for(int iP=0;iP<particles_XYZ_thisX.size();iP++) {
					j1x[jt] += (particles_XYZ_0[iP].v_c1+v1[iP][jt][0].c1)*particles_XYZ_0[iP].weight;
					//j1x[jt] -= (particles_XYZ_0[iP].v_c1)*particles_XYZ_0[iP].weight;
			}
			j1x[jt] = j1x[jt] * qe;
			std::cout << "j1x["<<jt<<"]: "<< j1x[jt]<<std::endl;

		}

#if __SAVE_ORBITS__>=1
		// Write orbits to file
	
		std::cout << "Writing orbits to file ... " << std::endl;

		std::stringstream ncOrbitsFileName;
		ncOrbitsFileName << "output/orbits_";
		ncOrbitsFileName << std::setw(3) << std::setfill('0') << iX;
	   	ncOrbitsFileName << ".nc"; 	

		try {
				// Really need to fix this but I don't know how to 
				// write a vector of structures using netCDF yet.

				netCDF::NcFile ncOrbitsFile (ncOrbitsFileName.str().c_str(), netCDF::NcFile::replace);
		
				netCDF::NcDim nc_nP = ncOrbitsFile.addDim("nP", particles_XYZ_thisX.size());
				netCDF::NcDim nc_nSteps = ncOrbitsFile.addDim("nSteps", nSteps);
				netCDF::NcDim nc_nJp = ncOrbitsFile.addDim("nJp", nJp);
	
				std::vector<netCDF::NcDim> nc_nPxnSteps(2);
				nc_nPxnSteps[0]=nc_nP;
				nc_nPxnSteps[1]=nc_nSteps;

				std::vector<netCDF::NcDim> nc_nPxnJpxnSteps(3);
				nc_nPxnJpxnSteps[0]=nc_nP;
				nc_nPxnJpxnSteps[1]=nc_nJp;
				nc_nPxnJpxnSteps[2]=nc_nSteps;
		
				netCDF::NcVar nc_t = ncOrbitsFile.addVar("t",netCDF::ncFloat,nc_nSteps);
		
				netCDF::NcVar nc_x = ncOrbitsFile.addVar("x",netCDF::ncFloat,nc_nPxnSteps);
				netCDF::NcVar nc_y = ncOrbitsFile.addVar("y",netCDF::ncFloat,nc_nPxnSteps);
				netCDF::NcVar nc_z = ncOrbitsFile.addVar("z",netCDF::ncFloat,nc_nPxnSteps);
		
				netCDF::NcVar nc_e1_x = ncOrbitsFile.addVar("e1_x",netCDF::ncFloat,nc_nPxnSteps);
				netCDF::NcVar nc_e1_y = ncOrbitsFile.addVar("e1_y",netCDF::ncFloat,nc_nPxnSteps);
				netCDF::NcVar nc_e1_z = ncOrbitsFile.addVar("e1_z",netCDF::ncFloat,nc_nPxnSteps);
		
				netCDF::NcVar nc_v1_x = ncOrbitsFile.addVar("v1x",netCDF::ncFloat,nc_nPxnJpxnSteps);
				netCDF::NcVar nc_v1_y = ncOrbitsFile.addVar("v1y",netCDF::ncFloat,nc_nPxnJpxnSteps);
				netCDF::NcVar nc_v1_z = ncOrbitsFile.addVar("v1z",netCDF::ncFloat,nc_nPxnJpxnSteps);
		
				std::vector<size_t> startpA(2);
				std::vector<size_t> countpA(2);
				for(int iP=0;iP<particles_XYZ_thisX.size();iP++) {
		
						startpA[0]=iP;
						startpA[1]=0;
						countpA[0]=1;
						countpA[1]=nSteps;
		
						std::vector<float> tmpData (nSteps,0);
						for(int iS=0;iS<nSteps;iS++){tmpData[iS] = orbit[iP][iS].c1;}
						nc_x.putVar(startpA,countpA,&tmpData[0]);
						for(int iS=0;iS<nSteps;iS++){tmpData[iS] = orbit[iP][iS].c2;}
						nc_y.putVar(startpA,countpA,&tmpData[0]);
						for(int iS=0;iS<nSteps;iS++){tmpData[iS] = orbit[iP][iS].c3;}
						nc_z.putVar(startpA,countpA,&tmpData[0]);

						for(int iS=0;iS<nSteps;iS++){tmpData[iS] = e1[iP][iS].c1;}
						nc_e1_x.putVar(startpA,countpA,&tmpData[0]);
						for(int iS=0;iS<nSteps;iS++){tmpData[iS] = e1[iP][iS].c2;}
						nc_e1_y.putVar(startpA,countpA,&tmpData[0]);
						for(int iS=0;iS<nSteps;iS++){tmpData[iS] = e1[iP][iS].c3;}
						nc_e1_z.putVar(startpA,countpA,&tmpData[0]);

				}

				std::vector<size_t> startpB(3);
				std::vector<size_t> countpB(3);
				for(int iP=0;iP<particles_XYZ_thisX.size();iP++) {
						for(int iJ=0;iJ<nJp;iJ++) {
		
							startpB[0]=iP;
							startpB[1]=iJ;
							startpB[2]=0;
							countpB[0]=1;
							countpB[1]=1;
							countpB[2]=nSteps;
		
							std::vector<float> tmpData (nSteps,0);

							for(int iS=0;iS<nSteps;iS++){tmpData[iS] = v1[iP][iJ][iS].c1;}
							nc_v1_x.putVar(startpB,countpB,&tmpData[0]);
							for(int iS=0;iS<nSteps;iS++){tmpData[iS] = v1[iP][iJ][iS].c2;}
							nc_v1_y.putVar(startpB,countpB,&tmpData[0]);
							for(int iS=0;iS<nSteps;iS++){tmpData[iS] = v1[iP][iJ][iS].c3;}
							nc_v1_z.putVar(startpB,countpB,&tmpData[0]);

						}
				}
		
				std::vector<size_t> startp (1,0);
				std::vector<size_t> countp (1,nSteps);
		
				nc_t.putVar(startp,countp,&t[0]);
		
		
		}
				catch(netCDF::exceptions::NcException &e) {
						std::cout << "NetCDF: unknown error" << std::endl;
						e.what();
		}

		std::cout << "DONE" << std::endl;
#endif

		// Write current to file
	
		std::cout << "Writing jP to file ... ";

		std::stringstream ncjPFileName;
		ncjPFileName << "output/jP_";
		ncjPFileName << std::setw(3) << std::setfill('0') << iX;
	   	ncjPFileName << ".nc"; 	
		netCDF::NcFile ncjPFile (ncjPFileName.str().c_str(), netCDF::NcFile::replace);

		netCDF::NcDim nc_nJp = ncjPFile.addDim("nJp", nJp);
		netCDF::NcDim nc_scalar = ncjPFile.addDim("scalar", 1);

		netCDF::NcVar nc_t = ncjPFile.addVar("t",netCDF::ncFloat,nc_nJp);

		netCDF::NcVar nc_x = ncjPFile.addVar("x",netCDF::ncFloat,nc_scalar);

		netCDF::NcVar nc_j1x = ncjPFile.addVar("j1x",netCDF::ncFloat,nc_nJp);
		netCDF::NcVar nc_j1y = ncjPFile.addVar("j1y",netCDF::ncFloat,nc_nJp);
		netCDF::NcVar nc_j1z = ncjPFile.addVar("j1z",netCDF::ncFloat,nc_nJp);

		nc_x.putVar(&xGrid[iX]);

		std::vector<size_t> startp (1,0);
		std::vector<size_t> countp (1,nJp);

		nc_t.putVar(startp,countp,&tJp[0]);

		nc_j1x.putVar(startp,countp,&j1x[0]);
		nc_j1y.putVar(startp,countp,&j1y[0]);
		nc_j1z.putVar(startp,countp,&j1z[0]);


	} // End of xGrid loop

	//ProfilerStop();

	std::cout << "DONE" << std::endl;

	return EXIT_SUCCESS;
}
