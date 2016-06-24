#include "c3vec.hpp"
#include <algorithm>
#include <fstream>
#include <iostream>
#include <netcdf>
#include <new> // std::bad_alloc
#include <string>
#include <vector>

using namespace netCDF;
using namespace exceptions;

int read_e_field(std::string eField_fName, int& species_number, float& freq,
    std::vector<float>& r, std::vector<float>& n_m3,
    std::vector<C3<std::complex<float> > >& e1_CYL, std::vector<C3<std::complex<float> > >& b1_CYL,
    std::vector<C3<float> >& e1Re_CYL, std::vector<C3<float> >& e1Im_CYL,
    std::vector<C3<float> >& b1Re_CYL, std::vector<C3<float> >& b1Im_CYL,
    std::vector<C3<float> >& b0_CYL)
{

    std::cout << "Reading eField data file " << eField_fName << std::endl;

    // Here we are using the cxx-4 netcdf interface by Lynton Appel
    // This needs netCDF 4.1.1 or later build with
    // ./configure --enable-cxx-4 [plus other options]

    std::vector<float> b0_r, b0_p, b0_z,
        e_r_re, e_p_re, e_z_re,
        e_r_im, e_p_im, e_z_im,
        b_r_re, b_p_re, b_z_re,
        b_r_im, b_p_im, b_z_im;

    std::vector<std::complex<float> > e_r, e_p, e_z;
    std::vector<std::complex<float> > b_r, b_p, b_z;

    std::ifstream file(eField_fName.c_str());
    if (!file.good()) {
        std::cout << "ERROR: Cannot find file " << eField_fName << std::endl;
        exit(1);
    }

    try {
        std::cout << "Reading E field data file ... " << eField_fName << std::endl;
        NcFile dataFile(eField_fName.c_str(), NcFile::read);

        NcDim nc_nR(dataFile.getDim("nR"));
        NcDim nc_nSpec(dataFile.getDim("nSpec"));
        NcDim nc_scalar(dataFile.getDim("scalar"));

        int nR = nc_nR.getSize();
        int nSpec = nc_nSpec.getSize();

        if (species_number > nSpec - 1) {
            std::cout << "ERROR: Asking for species that does not exist in density data" << std::endl;
            exit(1);
        }

        std::cout << "\tnR: " << nR << std::endl;

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
        std::vector<size_t> start, count;
        start.resize(2);
        count.resize(2);
        start[1] = 0;
        start[0] = species_number;
        count[1] = nR;
        count[0] = 1;

        nc_density.getVar(start, count, &n_m3[0]);

        try {
            std::cout << "nR : " << nR << std::endl;
            b0_CYL.resize(nR);
        } catch (const std::bad_alloc& error) {
            std::cout << "Allocation error at " << __FILE__ << __LINE__ << std::endl;
            std::cout << error.what();
        }

        for (int i = 0; i < nR; i++) {
            b0_CYL[i] = C3<float>(b0_r[i], b0_p[i], b0_z[i]);
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

        for (int i = 0; i < nR; i++) {
            e_r.push_back(std::complex<float>(e_r_re[i], e_r_im[i]));
            e_p.push_back(std::complex<float>(e_p_re[i], e_p_im[i]));
            e_z.push_back(std::complex<float>(e_z_re[i], e_z_im[i]));
        }

        for (int i = 0; i < nR; i++) {
            b_r.push_back(std::complex<float>(b_r_re[i], b_r_im[i]));
            b_p.push_back(std::complex<float>(b_p_re[i], b_p_im[i]));
            b_z.push_back(std::complex<float>(b_z_re[i], b_z_im[i]));
        }

        std::vector<float>::iterator min = std::min_element(b0_p.begin(), b0_p.end());
        std::vector<float>::iterator max = std::max_element(b0_p.begin(), b0_p.end());
#if DEBUGLEVEL >= 1
        std::cout << "\tR[0]: " << r[0] << ", R[" << nR << "]: " << r[r.size() - 1] << std::endl;
        std::cout << "\tfreq: " << freq << std::endl;
        std::cout << "\tmin(b0_p): " << *min << std::endl;
        std::cout << "\tmax(b0_p): " << *max << std::endl;
        std::cout << "\tabs(e_r[nR/2]): " << abs(e_r[nR / 2]) << std::endl;
        std::cout << "\tabs(e_p[nR/2]): " << abs(e_p[nR / 2]) << std::endl;
        std::cout << "\tabs(e_z[nR/2]): " << abs(e_z[nR / 2]) << std::endl;
#endif
    } catch (exceptions::NcException& e) {
        std::cout << "NetCDF: unknown error" << std::endl;
        e.what();
        exit(1);
    }

    e1Re_CYL.resize(e_r.size());
    e1Im_CYL.resize(e_r.size());
    b1Re_CYL.resize(e_r.size());
    b1Im_CYL.resize(e_r.size());

    e1_CYL.resize(e_r.size());
    b1_CYL.resize(b_r.size());

    for (int i = 0; i < e_r.size(); i++) {

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

        e1_CYL[i].c1 = e_r[i];
        e1_CYL[i].c2 = e_p[i];
        e1_CYL[i].c3 = e_z[i];

        b1_CYL[i].c1 = b_r[i];
        b1_CYL[i].c2 = b_p[i];
        b1_CYL[i].c3 = b_z[i];
    }

    std::cout << "End of " << __FILE__ << std::endl;
    return (0);
}
