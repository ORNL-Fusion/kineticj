#include <string>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <cstdlib>
#include <netcdf>
#include <vector>
#include <algorithm>
#include <complex>
#include <libconfig.h++>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <assert.h>

using namespace std;
using namespace netCDF;
using namespace exceptions;

class C3Vec {
		public:
				float c1, c2, c3;

                C3Vec (int _const) {c1=_const;c2=_const;c3=_const;};
				C3Vec () {c1=0;c2=0;c3=0;};
				C3Vec ( float _c1, float _c2, float _c3 ) {c1=_c1;c2=_c2;c3=_c3;};
				C3Vec& operator = (const C3Vec &rhs);
};

C3Vec& C3Vec::operator= (const C3Vec &rhs ) {
		if (this != &rhs) {
				c1 = rhs.c1;
				c2 = rhs.c2;
				c3 = rhs.c3;
		}
		return *this;
}

int main ( int argc, char **argv )
{
		libconfig::Config cfg;
		string cfgName = "kj.cfg";
		
		// Open the config file
		cfg.readFile(cfgName.c_str());
	    int species_number = cfg.lookup("species_number");

		// Read E
		string eField_fName = cfg.lookup("eField_fName");	
		cout << "Reading eField data file " << eField_fName << endl;

		// Here we are using the cxx-4 netcdf interface by Lynton Appel
		// This needs netCDF 4.1.1 or later build with
		// ./configure --enable-cxx-4 [plus other options]
		ifstream file(eField_fName.c_str());
		if(!file.good()) {
			cout << "ERROR: Cannot find file " << eField_fName << endl;
			exit(1);
		}

#if DIM == 1
            vector<float>  r, b0_r, b0_p, b0_z;
		try {
				NcFile dataFile ( eField_fName.c_str(), NcFile::read );
	
				NcDim nc_nR(dataFile.getDim("nR"));
				NcDim nc_scalar(dataFile.getDim("scalar"));
            
                int nR = nc_nR.getSize();
                int nZ = 1;

                cout << "\tnR: " << nR << endl;

                NcVar nc_r(dataFile.getVar("r"));
            
                NcVar nc_b0_r(dataFile.getVar("B0_r"));
                NcVar nc_b0_p(dataFile.getVar("B0_p"));
                NcVar nc_b0_z(dataFile.getVar("B0_z"));

                r.resize(nR);
            
                b0_r.resize(nR);
                b0_p.resize(nR);
                b0_z.resize(nR);

                nc_r.getVar(&r[0]);

                nc_b0_r.getVar(&b0_r[0]);
                nc_b0_p.getVar(&b0_p[0]);
                nc_b0_z.getVar(&b0_z[0]);
            }
		catch(exceptions::NcException &e) {
				cout << "NetCDF: unknown error" << endl;
				e.what();
				exit(1);
		}
#endif

#if DIM == 2
            vector<float> r, z;
			vector<vector <float> > b0_r, b0_p, b0_z;
 
		try {
				NcFile dataFile ( eField_fName.c_str(), NcFile::read );
	
				NcDim nc_nR(dataFile.getDim("nR"));
				NcDim nc_nSpec(dataFile.getDim("nSpec"));
				NcDim nc_scalar(dataFile.getDim("scalar"));

                NcDim nc_nZ(dataFile.getDim("nZ"));
    
                int nR = nc_nR.getSize();
                int nZ = nc_nZ.getSize();
                int nSpec = nc_nSpec.getSize();

                cout << "\tnR: " << nR << endl;
                cout << "\tnZ: " << nZ << endl;

                NcVar nc_r(dataFile.getVar("r"));
                NcVar nc_z(dataFile.getVar("z"));
            
                NcVar nc_b0_r = dataFile.getVar("B0_r");
                NcVar nc_b0_p = dataFile.getVar("B0_p");
                NcVar nc_b0_z = dataFile.getVar("B0_z");

                r.resize(nR);
                z.resize(nZ);

				for(int i=0; i<nR; i++)
				{
					vector<float> row(nZ);
					b0_r.push_back(row);
					b0_p.push_back(row);
					b0_z.push_back(row);
				}
       
                nc_r.getVar(&r[0]);
            
				//float* b0_r_2D = new float[nR*nZ];
				//float* b0_p_2D = new float[nR*nZ];
				//float* b0_z_2D = new float[nR*nZ];

				float b0_r_2D[nR][nZ];
				float b0_p_2D[nR][nZ];
				float b0_z_2D[nR][nZ];
            
				// Index [i][j] as [i*nZ+j]

				nc_b0_r.getVar(b0_r_2D);
				nc_b0_p.getVar(b0_p_2D);
				nc_b0_z.getVar(b0_z_2D);

				for (int i=0; i<nR; ++i) {
					for (int j=0; j<nZ; ++j) {

//						b0_r[i][j] = b0_r_2D[i*nZ+j];
//						b0_p[i][j] = b0_p_2D[i*nZ+j];
//						b0_z[i][j] = b0_z_2D[i*nZ+j];

						b0_r[i][j] = b0_r_2D[i][j];
						b0_p[i][j] = b0_p_2D[i][j];
						b0_z[i][j] = b0_z_2D[i][j];

						cout << b0_r[i][j] << "  " << b0_p[i][j] << "  " << b0_z[i][j] << endl;

					}
				}

//				delete [] b0_r_2D;
//				delete [] b0_p_2D;
//				delete [] b0_z_2D;

            }

		catch(exceptions::NcException &e) {
				cout << "NetCDF: unknown error" << endl;
				e.what();
				exit(1);
		}
#endif

    return 0;
}
