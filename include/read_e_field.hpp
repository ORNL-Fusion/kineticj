#ifndef READ_E_FIELD_HPP
#define READ_E_FIELD_HPP

#include <string>
#include <netcdf>
#include <vector>
#include "c3vec.hpp"
#include <iostream>
#include <fstream>

using namespace netCDF;
using namespace exceptions;

int read_e_field( std::string eField_fName, int &species_number, float &freq, 
                std::vector<float> &r, std::vector<float> &n_m3, 
                std::vector<C3VecI> &e1_CYL, std::vector<C3VecI> &b1_CYL, 
                std::vector<C3Vec> &e1Re_CYL,std::vector<C3Vec> &e1Im_CYL,
                std::vector<C3Vec> &b1Re_CYL,std::vector<C3Vec> &b1Im_CYL, 
                std::vector<C3Vec> &b0_CYL);

#endif
