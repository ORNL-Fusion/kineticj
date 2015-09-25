#ifndef READ_GC_FILE_HPP
#define READ_GC_FILE_HPP

#include <string>
#include <netcdf>
#include <vector>
#include "c3vec.hpp"

int read_gc_file( std::string fName, std::vector<float> &r_gc, std::vector<C3Vec> &curv_CYL, std::vector<C3Vec> &grad_CYL, 
               std::vector<float> &bDotGradB );

#endif

