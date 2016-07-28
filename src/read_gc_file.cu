#include "c3vec.hpp"
#include "cspecies.hpp"
#include <fstream>
#include <iostream>
#include <netcdf>
#include <string>
#include <vector>

using namespace netCDF;
using namespace exceptions;

int read_gc_file(std::string fName,
    std::vector<float>& r_gc, std::vector<C3<float> >& curv_CYL, std::vector<C3<float> >& grad_CYL,
    std::vector<float>& bDotGradB)
{

    // Read the guiding center terms from file
    std::cout << "Reading GC terms data file " << fName << std::endl;

    std::vector<float> curv_r, curv_p, curv_z,
        grad_r, grad_p, grad_z;

    std::ifstream gc_file(fName.c_str());
    if (!gc_file.good()) {
        std::cout << "ERROR: Cannot find file " << fName << std::endl;
        return (1);
    }

    NcFile dataFile(fName.c_str(), NcFile::read);
    NcDim gc_nc_nR(dataFile.getDim("nR"));
    NcDim gc_nc_scalar(dataFile.getDim("scalar"));
    if (!dataFile.getVar("z").isNull())
        throw NcException("NcException", "This is a 2D GC terms file", __FILE__, __LINE__);

    int nR_gc = gc_nc_nR.getSize();
    std::cout << "nR_gc: " << nR_gc << std::endl;
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
    for (int i = 0; i < nR_gc; i++) {
        curv_CYL[i] = C3<float>(curv_r[i], curv_p[i], curv_z[i]);
        grad_CYL[i] = C3<float>(grad_r[i], grad_p[i], grad_z[i]);
    }
    std::cout << "Finished reading gc_terms file" << std::endl;

    return (0);
}
