#include "cspecies.hpp"

CSpecies::CSpecies(const double _amu, const int _Z)
{

    amu = _amu;
    Z = _Z;
    m = amu * physConstants::amu;
    q = Z * physConstants::e;
}

