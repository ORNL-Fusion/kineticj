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
vector<CParticle> create_particles ( float x, float amu, float Z, float T_keV, float n_m3, int nPx, int nPy, int nPz, int nThermal, float &dv, float *r, C3Vec *b0_CYL, int nX);

// Functor to wrap df0_dv 

struct get_df0_dv
{
    C3Vec operator() (CParticle &p) {
        C3Vec thisVel_XYZ(p.v_c1, p.v_c2, p.v_c3);
        C3Vec gradv_f0_XYZ = maxwellian_df0_dv(thisVel_XYZ, p.T, p.n, p.amu, p.Z);
        return gradv_f0_XYZ;
    }
};

struct set_vx
{
#ifdef __CUDACC__
    __host__ __device__
#endif
    float operator() (float &x, CParticle &p) {
        return p.v_c1;
    }
};

struct set_vy
{
#ifdef __CUDACC__
    __host__ __device__
#endif
    float operator() (float &x, CParticle &p) {
        return p.v_c2;
    }
};

struct set_vz
{
#ifdef __CUDACC__
    __host__ __device__
#endif
   float operator() (float &x, CParticle &p) {
        return p.v_c3;
    }
};

#endif
