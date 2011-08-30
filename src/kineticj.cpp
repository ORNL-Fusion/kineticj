#include <iostream>
#include <cstdlib>
#include <netcdf>
#include <vector>
#include <algorithm>
#include <complex>

// Calculate the jP given some know E and f(v)

int main ( int argc, char **argv )
{
		// Read E
	
		std::string rsfwc_fName ( "data/rsfwc_1d.nc" );	
		std::cout << "Reading rsfwc data file" << rsfwc_fName << std::endl;

		// Here we are using the cxx-4 netcdf interface by Lynton Appel
		// This needs netCDF 4.1.1 or later build with
		// ./configure --enable-cxx-4 [plus other options]

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

				std::vector<float> r(nR), b0_r(nR), b0_p(nR), b0_z(nR),
						e_r_re(nR), e_p_re(nR), e_z_re(nR),
						e_r_im(nR), e_p_im(nR), e_z_im(nR);

				float wrf;

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

				std::vector<std::complex<float> > e_r, e_p, e_z;	
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

		// Generate linear orbits

		// Create f1(v)

		// Calculate jP1

		return EXIT_SUCCESS;
}
