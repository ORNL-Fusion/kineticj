#ifndef CREATEPARTICLES_HPP
#define CREATEPARTICLES_HPP

#include "c3vec.hpp"
#include "cparticle.hpp"
#include "interp.hpp"
#include "rotation.hpp"

using namespace std;

float GetGyroPhase ( const C3Vec v_abp );
float maxwellian ( float vx, float vy, float vz, float vTh );
float get_vTh ( const float _amu, const float _Z, const float _T_keV );
C3Vec maxwellian_df0_dv (const C3Vec _v, const float _T_keV, const float _n_m3, const float _amu, const float _Z );
vector<CParticle> create_particles ( float x, float amu, float Z, float T_keV, float n_m3, int nPx, int nPy, int nPz, int nThermal, float &dv, vector<float> &r, vector<C3Vec> &b0_CYL);

#endif
