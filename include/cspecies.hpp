#ifndef CSPECIES_HPP
#define CSPECIES_HPP
#include <string>
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

#endif
