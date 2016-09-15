#ifndef CSPECIES_HPP
#define CSPECIES_HPP
#include <string>
#include "constants.hpp"

#ifdef __CUDACC__
#define HOST __host__ 
#define DEVICE __device__
#else
#define HOST 
#define DEVICE
#endif

#ifdef __CUDACC__
#define PRAGMA #pragma hd_warning_disable 
#else
#define PRAGMA
#endif

class CSpecies {
		public:
				double m, q;
				double amu;
				int Z;
                //std::string name;

                PRAGMA
                HOST DEVICE
				CSpecies () {};

                PRAGMA
                HOST DEVICE
				CSpecies ( const double amu, const int Z);
};

#endif
