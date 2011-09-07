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
				C3Vec& operator *= (const C3Vec &rhs);
				C3Vec& operator *= (const float &rhs);
				C3Vec& operator /= (const C3Vec &rhs);
				C3Vec& operator /= (const float &rhs);

				C3Vec operator + (const C3Vec &other);
				C3Vec operator + (const float &other);
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

C3Vec rk4_evalf ( CParticle &p, float t, const C3Vec &k, const C3Vec &b0, const C3VecI &e1, const float wrf ) {

	C3Vec f;

	C3Vec v_x_b0 ( p.v_c2*b0.c3-p.v_c3*b0.c2, -1.0*(p.v_c1*b0.c3-p.v_c3*b0.c1), p.v_c1*b0.c2-p.v_c2*b0.c1); 
	C3Vec ( std::real(e1.c1) * cos ( wrf * t ) + std::imag(e1.c1) * sin ( wrf * t ) + v_x_b0.c1,
		  	std::real(e1.c2) * cos ( wrf * t ) + std::imag(e1.c2) * sin ( wrf * t ) + v_x_b0.c2,
		  	std::real(e1.c3) * cos ( wrf * t ) + std::imag(e1.c3) * sin ( wrf * t ) + v_x_b0.c3 );

	return f*(p.q/p.m);	
}

void rk4_move ( CParticle &p, float dt, float t0, 
				const C3Vec &b0, const C3VecI &e1, const float wrf ) {

		C3Vec yn0(p.v_c1,p.v_c2,p.v_c3), k1, k2, k3, k4, yn1; 

		k1 = rk4_evalf ( p, t0 + 0.0*dt, yn0 + 0.*yn0, b0, e1, wrf ) * dt;	
		k2 = rk4_evalf ( p, t0 + 0.5*dt, yn0 + 0.5*k1, b0, e1, wrf ) * dt;	
		k3 = rk4_evalf ( p, t0 + 0.5*dt, yn0 + 0.5*k2, b0, e1, wrf ) * dt;	
		k4 = rk4_evalf ( p, t0 + 1.0*dt, yn0 + 1.0*k3, b0, e1, wrf ) * dt;	

		yn1 = yn0 + 1.0/6.0 * (k1+2.0*k2+2.0*k3+k4);

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

		std::string particleList_fName ( "data/f_0.02keV_electrons.nc" );	
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

		std::cout << "Generating linear orbit" << std::endl;

		int nRFCycles = 10;
		float tRF = (2*_pi)/wrf;
		float dtMin = tRF/10.0;
		float tEnd = tRF * nRFCycles;

		std::vector<CParticle> orbit;

		float t=0;
		int nSteps = 0;
		int iP = 10;
		while(t<tEnd) {
				orbit.push_back(particles_XYZ[iP]);
				t+=dtMin;
				nSteps++;
		}

		std::cout << "\tnSteps: " << nSteps << std::endl;

		// Create f1(v) by integrating F to give dv

		std::vector<C3Vec> dv(orbit.size());	
		std::vector<C3Vec> e_t(orbit.size());

		for(int i=0;i<dv.size();i++) {

				float er = std::real(e_r[i])*cos(wrf*t)+std::imag(e_r[i])*sin(wrf*t);
				float ep = std::real(e_p[i])*cos(wrf*t)+std::imag(e_p[i])*sin(wrf*t);
				float ez = std::real(e_z[i])*cos(wrf*t)+std::imag(e_z[i])*sin(wrf*t);

				e_t[i] = C3Vec(er,ep,ez);	

				if(i>0) {
					dv[i] = dv[i-1]+(e_t[i]+e_t[i-1])/2*dtMin*(particles_XYZ[iP].q/particles_XYZ[iP].m);
					std::cout << dv[i].c1 << "  " << dv[i].c2 << "  " << dv[i].c3 << std::endl ;
				}
		}

		// Calculate jP1

		float dvGuess = 3e8*0.01/50.0;
		std::vector<C3Vec> jP1(dv.size());
		for(int i=0;i<dv.size();i++) {
				jP1[i] = dv[i]*1.0e18*dvGuess;
				std::cout << jP1[i].c1 << "  " << jP1[i].c2 << "  " << jP1[i].c3 << std::endl ;
		}

		return EXIT_SUCCESS;
}
