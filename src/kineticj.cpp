#include <string>
#include <iostream>
#include <cstdlib>
#include <netcdf>
#include <vector>
#include <algorithm>
#include <complex>
#include "constants.hpp"

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

				CParticle ();
				CParticle ( double _amu, int _Z);
				CParticle (float c1, float c2, float c3, float v_c1, float v_c2, float v_c3, double _amu, int _Z );
				CParticle (CSpecies _species);
};

CParticle::CParticle () {
}

CParticle::CParticle ( double _amu, int _Z ): CSpecies(_amu,_Z) {
}

CParticle::CParticle 
	(float _c1, float _c2, float _c3, float _v_c1, float _v_c2, float _v_c3, double _amu, int _Z ) :
	CSpecies (_amu, _Z) {

	c1 = _c1;
	c2 = _c2;
	c3 = _c3;

	v_c1 = _v_c1;
	v_c2 = _v_c2;
	v_c3 = _v_c3;
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

// Zero-order orbits
C3Vec rk4_evalf ( CParticle &p, const float &t, 
				const C3Vec &v_XYZ, const C3Vec &x, const std::vector<C3Vec> &b0Vec_CYL,
			  	const std::vector<float> &rVec ) {

	// Interpolate b0 at location in CYL
	
	float _r = sqrt ( pow(x.c1,2) + pow(x.c2,2) );
	float _p = atan2 ( x.c2, x.c1 );

	//std::cout << "\tx: " << x.c1 << " y: " << x.c2 << " z: " << x.c3 << std::endl;
	//std::cout << "\tr: " << _r << " p: " << _p << std::endl;
	//std::cout << "\trVec.front(): " << rVec.front() << std::endl;
	//std::cout << "\tv_XYZ: " << v_XYZ.c1 << "  " << v_XYZ.c2 << "  " << v_XYZ.c3 << std::endl;

	float _x = (_r-rVec.front())/(rVec.back()-rVec.front())*(rVec.size()-1);
	float x0 = floor(_x);
	float x1 = ceil(_x);

	//std::cout << "\t_x: " << _x << " x0: " << x0 << " x1: " << x1 << std::endl;

	C3Vec b0_CYL, b0_XYZ;

	if(x0>=0 && x1<=b0Vec_CYL.size()-1) {

		C3Vec y0 = b0Vec_CYL[x0];
		C3Vec y1 = b0Vec_CYL[x1];

		// Linear interpolation
		b0_CYL = y0+(_x-x0)*(y1-y0)/(x1-x0);
		b0_XYZ = C3Vec( cos(_p)*b0_CYL.c1-sin(_p)*b0_CYL.c2+0,
						sin(_p)*b0_CYL.c1+cos(_p)*b0_CYL.c2+0,
						0+0+1*b0_CYL.c3 );

		//std::cout << "\tb0_XYZ: " << b0_XYZ.c1 << "  " << b0_XYZ.c2 << "  " << b0_XYZ.c3 << std::endl;
	}
	else {
		std::cout << "\tERROR: off grid." << std::endl;
	}

	C3Vec v_x_b0 ( v_XYZ.c2*b0_XYZ.c3-v_XYZ.c3*b0_XYZ.c2, 
					-1.0*(v_XYZ.c1*b0_XYZ.c3-v_XYZ.c3*b0_XYZ.c1), 
					v_XYZ.c1*b0_XYZ.c2-v_XYZ.c2*b0_XYZ.c1);

	//std::cout << "\tvxb0: " << v_x_b0.c1 << "  " << v_x_b0.c2 << "  " << v_x_b0.c3 << std::endl;

	//std::cout << "\tp.q/p.m: " << p.q/p.m << std::endl;

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

		//std::cout << "\tx0_XYZ: " << xn0.c1 << "  " << xn0.c2 << "  " << xn0.c3 << std::endl;
		//std::cout << "\tv0_XYZ: " << yn0.c1 << "  " << yn0.c2 << "  " << yn0.c3 << std::endl;
		//std::cout << "\tx1_XYZ: " << xn1.c1 << "  " << xn1.c2 << "  " << xn1.c3 << std::endl;
		//std::cout << "\tv1_XYZ: " << yn1.c1 << "  " << yn1.c2 << "  " << yn1.c3 << std::endl;
		//std::cout << "\tE: " << 0.5 * p.m * sqrt (pow(p.v_c1,2)+pow(p.v_c2,2)+pow(p.v_c3,2))/_e << std::endl;
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

		std::string particleList_fName ( "data/f_20keV_electrons.nc" );	
		std::cout << "Reading particle list " << particleList_fName << std::endl;

		std::vector<float> p_x, p_y, p_z, p_vx, p_vy, p_vz, p_amu;
		std::vector<int> p_Z;
		
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

				p_x.resize(nP);
				p_y.resize(nP);
				p_z.resize(nP);

				p_vx.resize(nP);
				p_vy.resize(nP);
				p_vz.resize(nP);

				p_amu.resize(nP);
				p_Z.resize(nP);

				nc_p_x.getVar(&p_x[0]);
				nc_p_y.getVar(&p_y[0]);
				nc_p_z.getVar(&p_z[0]);

				nc_p_vx.getVar(&p_vx[0]);
				nc_p_vy.getVar(&p_vy[0]);
				nc_p_vz.getVar(&p_vz[0]);

				nc_p_amu.getVar(&p_amu[0]);
				nc_p_Z.getVar(&p_Z[0]);

		}
		catch(netCDF::exceptions::NcException &e) {
				std::cout << "NetCDF: unknown error" << std::endl;
				e.what();
		}

		std::vector<CParticle> particles_XYZ;
		particles_XYZ.resize(p_x.size());

		for(int i=0;i<particles_XYZ.size();i++){

				CParticle thisParticle (p_x[i],p_y[i],p_z[i],p_vx[i],p_vy[i],p_vz[i],p_amu[i],p_Z[i]);
				particles_XYZ[i] = thisParticle;

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

		// Generate linear orbits

		std::cout << "Generating linear orbit with RK4" << std::endl;

		int nRFCycles = 10;
		int nStepsPerCycle = 100;
		float tRF = (2*_pi)/wrf;
		float dtMin = tRF/nStepsPerCycle;
		float tEnd = tRF * nRFCycles;

		std::vector<std::vector<C3Vec> > orbit(particles_XYZ.size());
		std::vector<float> t;

		int nSteps = nRFCycles*nStepsPerCycle;
		t.resize(nSteps);

		for(int iP=0;iP<10;iP++) {

			orbit[iP].resize(nSteps);

	 		for(int i=0;i<nSteps;i++) {	

					//std::cout << "\tE: " << 
					//		0.5 * particles_XYZ[iP].m * 
					//		sqrt (pow(particles_XYZ[iP].v_c1,2)
					//						+pow(particles_XYZ[iP].v_c2,2)
					//						+pow(particles_XYZ[iP].v_c3,2))/_e << std::endl;
					
					t[i]=i*dtMin;
					orbit[iP][i] = C3Vec(particles_XYZ[iP].c1,particles_XYZ[iP].c2,particles_XYZ[iP].c3);
					rk4_move ( particles_XYZ[iP], dtMin, t[i], b0_CYL, r );
			}
		}

		std::cout << "\tnSteps: " << nSteps << std::endl;

	std::vector<C3Vec> dv(nSteps);	
	std::vector<C3Vec> e1(nSteps);
	std::vector<C3Vec> v1(nSteps);

	for(int iP=0;iP<10;iP++) {

		// Create f1(v) by integrating F to give dv

		for(int i=0;i<e1.size();i++) {

				std::vector<C3Vec> e1Now_CYL;
			   	e1Now_CYL.resize(e_r.size());
				for(int j=0;j<e1Now_CYL.size();j++) {
					e1Now_CYL[j] = C3Vec ( 
								std::real(e_r[j])*cos(wrf*t[i])+std::imag(e_r[j])*sin(wrf*t[i]),
								std::real(e_p[j])*cos(wrf*t[i])+std::imag(e_p[j])*sin(wrf*t[i]),
								std::real(e_z[j])*cos(wrf*t[i])+std::imag(e_z[j])*sin(wrf*t[i]) );
				}

				// Interpolate e1Now to here, done in CYL
				
				float _r = sqrt ( pow(orbit[iP][i].c1,2) + pow(orbit[iP][i].c2,2) );
				float _p = atan2 ( orbit[iP][i].c2, orbit[iP][i].c1 );

				float _x = (_r-r.front())/(r.back()-r.front())*(r.size()-1);
				float x0 = floor(_x);
				float x1 = ceil(_x);

				C3Vec e1NowAndHere_CYL, e1NowAndHere_XYZ;

				if(x0>=0 && x1<=e1Now_CYL.size()-1) {

				 	// Linear interpolation
					C3Vec y0 = e1Now_CYL[x0];
					C3Vec y1 = e1Now_CYL[x1];
					e1NowAndHere_CYL = y0+(_x-x0)*(y1-y0)/(x1-x0);

					// Rotation CYL -> XYZ
					e1NowAndHere_XYZ = C3Vec( cos(_p)*e1NowAndHere_CYL.c1-sin(_p)*e1NowAndHere_CYL.c2+0,
									sin(_p)*e1NowAndHere_CYL.c1+cos(_p)*e1NowAndHere_CYL.c2+0,
									0+0+1*e1NowAndHere_CYL.c3 );
				}
				else {
					std::cout << "\tERROR: off grid." << std::endl;
				}

				e1[i] = e1NowAndHere_XYZ;

		}

		// Integrate acceleration along zero-order orbit to get a velocity delta

		v1[0].c1=0;v1[0].c2=0;v1[0].c3=0;

		for(int i=1;i<e1.size();i++) {

			v1[i] = v1[i-1] + particles_XYZ[iP].q/particles_XYZ[iP].m *
				(t[i]-t[i-1])/6.0	* (e1[i-1]+4*(e1[i-1]+e1[i])/2.0+e1[i]);

			//std::cout << "v1: " << v1[i].c1 << "  " << v1[i].c2 << "  " << v1[i].c3 << std::endl;
		}	

		std::cout << "\tmax dV: " << maxC3VecAbs(v1) << std::endl;

	}

		// Write orbits to file
	
		std::cout << "Writing orbits to file ... ";

		netCDF::NcFile ncOrbitsFile ("output/orbits.nc", netCDF::NcFile::replace);

		netCDF::NcDim nc_nP = ncOrbitsFile.addDim("nP", particles_XYZ.size());
		netCDF::NcDim nc_nSteps = ncOrbitsFile.addDim("nSteps", nSteps);

		netCDF::NcVar nc_x = ncOrbitsFile.addVar("x",netCDF::ncFloat,nc_nSteps);
		netCDF::NcVar nc_y = ncOrbitsFile.addVar("y",netCDF::ncFloat,nc_nSteps);
		netCDF::NcVar nc_z = ncOrbitsFile.addVar("z",netCDF::ncFloat,nc_nSteps);
		netCDF::NcVar nc_t = ncOrbitsFile.addVar("t",netCDF::ncFloat,nc_nSteps);
		netCDF::NcVar nc_e1_x = ncOrbitsFile.addVar("e1_x",netCDF::ncFloat,nc_nSteps);
		netCDF::NcVar nc_e1_y = ncOrbitsFile.addVar("e1_y",netCDF::ncFloat,nc_nSteps);
		netCDF::NcVar nc_e1_z = ncOrbitsFile.addVar("e1_z",netCDF::ncFloat,nc_nSteps);

	for(int iP=0;iP<1;iP++) {
		for(int i=0;i<orbit[iP].size();i++) {
			std::vector<size_t> index;
			index.resize(1);
			index[0]=i;
			nc_x.putVar(index,orbit[iP][i].c1);
			nc_y.putVar(index,orbit[iP][i].c2);
			nc_z.putVar(index,orbit[iP][i].c3);
			nc_t.putVar(index,t[i]);
			nc_e1_x.putVar(index,e1[i].c1);
			nc_e1_y.putVar(index,e1[i].c2);
			nc_e1_z.putVar(index,e1[i].c3);
		}
	}

		std::cout << "DONE" << std::endl;


		// Calculate jP1

		float dvGuess = 3e8*0.01/50.0;
		std::vector<C3Vec> jP1(dv.size());
		for(int i=0;i<dv.size();i++) {
				jP1[i] = dv[i]*1.0e18*dvGuess;
				//std::cout << jP1[i].c1 << "  " << jP1[i].c2 << "  " << jP1[i].c3 << std::endl ;
		}

		return EXIT_SUCCESS;
}
