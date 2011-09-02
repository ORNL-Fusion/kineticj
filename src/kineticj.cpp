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
				float r, p, z, v_r, v_p, v_z;

				CParticle ();
				CParticle ( double _amu, int _Z);
				CParticle (float r, float p, float z, float v_r, float v_p, float v_z, double _amu, int _Z );
				CParticle (CSpecies _species);
};

CParticle::CParticle () {
}

CParticle::CParticle ( double _amu, int _Z ): CSpecies(_amu,_Z) {
}

CParticle::CParticle 
	(float _r, float _p, float _z, float _v_r, float _v_p, float _v_z, double _amu, int _Z ) :
	CSpecies (_amu, _Z) {

	r = _r;
	p = _p;
	z = _z;

	v_r = _v_r;
	v_p = _v_p;
	v_z = _v_z;
}

CParticle::CParticle ( CSpecies _species ) : CSpecies(_species) {
}

class C3Vec {
		public:
				float r, p, z;

				C3Vec () {r=0;p=0;z=0;};
				C3Vec ( float _r, float _p, float _z ) {r=_r;p=_p;z=_z;};
				C3Vec operator + (C3Vec);
				C3Vec operator + (float addMe);
				C3Vec operator * (float factor);
				C3Vec operator / (float factor);
};

C3Vec C3Vec::operator+ (C3Vec addMe) {
		C3Vec tmp;
		tmp.r = r + addMe.r;
		tmp.p = p + addMe.p;
		tmp.z = z + addMe.z;
		return (tmp);
}

C3Vec C3Vec::operator+ (float addMe) {
		C3Vec tmp;
		tmp.r = r + addMe;
		tmp.p = p + addMe;
		tmp.z = z + addMe;
		return (tmp);
}

C3Vec C3Vec::operator* (float factor) {
		C3Vec tmp;
		tmp.r = r * factor;
		tmp.p = p * factor;
		tmp.z = z * factor;
		return (tmp);
}

C3Vec C3Vec::operator/ (float factor) {
		C3Vec tmp;
		tmp.r = r / factor;
		tmp.p = p / factor;
		tmp.z = z / factor;
		return (tmp);
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

		double amu = _me_mi;
		int Z = -1;
		CSpecies species(amu,Z);
		std::cout << species.m << "  " << species.q << "  " << std:: endl;
		CParticle particle(species);
		std::vector<CParticle> particles (r.size(),particle);

		for(int i=0;i<particles.size();i++){

				particles[i].r = r[i];
				particles[i].p = 0;
				particles[i].z = 0;

				particles[i].v_r = 0;
				particles[i].v_p = 0;
				particles[i].v_z = 0;

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
				orbit.push_back(particles[iP]);
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
					dv[i] = dv[i-1]+(e_t[i]+e_t[i-1])/2*dtMin*(particles[iP].q/particles[iP].m);
					std::cout << dv[i].r << "  " << dv[i].p << "  " << dv[i].z << std::endl ;
				}
		}

		// Calculate jP1

		float dvGuess = 3e8*0.01/50.0;
		std::vector<C3Vec> jP1(dv.size());
		for(int i=0;i<dv.size();i++) {
				jP1[i] = dv[i]*1.0e18*dvGuess;
				std::cout << jP1[i].r << "  " << jP1[i].p << "  " << jP1[i].z << std::endl ;
		}

		return EXIT_SUCCESS;
}
