#include "cparticle.hpp"

CParticle::CParticle()
{
    status = 0;
}

CParticle::CParticle(double _amu, int _Z)
    : CSpecies(_amu, _Z)
{
    status = 0;
}

CParticle::CParticle(float _c1, float _c2, float _c3,
    float _v_c1, float _v_c2, float _v_c3,
    double _amu, int _Z, float _weight)
    : CSpecies(_amu, _Z)
{

    c1 = _c1;
    c2 = _c2;
    c3 = _c3;

    v_c1 = _v_c1;
    v_c2 = _v_c2;
    v_c3 = _v_c3;

    weight = _weight;

    status = 0;
}

CParticle::CParticle(CSpecies _species)
    : CSpecies(_species)
{
}
