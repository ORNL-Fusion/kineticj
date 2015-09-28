#include "cspecies.hpp"

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
