#ifndef CREATEPARTICLES_HPP
#define CREATEPARTICLES_HPP

#include "c3vec.hpp"
#include "cparticle.hpp"
#include "interp.hpp"
#include "rotation.hpp"

#if defined(__CUDACC__) 
#define HOST __host__ 
#define DEVICE __device__
#else
#define HOST 
#define DEVICE
#endif

#if defined(__CUDACC__) 
#define PRAGMA #pragma hd_warning_disable 
#else
#define PRAGMA
#endif

using namespace std;

float GetGyroPhase ( const C3<float> v_abp );

PRAGMA
HOST DEVICE
float maxwellian ( float vx, float vy, float vz, float vTh );

PRAGMA
HOST DEVICE
float get_vTh ( const float _amu, const float _Z, const float _T_keV );

PRAGMA
HOST DEVICE
C3<float> maxwellian_df0_dv (const C3<float> _v, const float _T_keV, const float _n_m3, const float _amu, const float _Z );

vector<CParticle> create_particles ( float x, float amu, float Z, float T_keV, float n_m3, int nPx, int nPy, int nPz, int nThermal, float &dv, float *r, C3<float> *b0_CYL, int nX);

// Functor to wrap df0_dv 

struct get_df0_dv
{
    PRAGMA
    HOST DEVICE
    C3<float> operator() (CParticle &p) {
        C3<float> thisVel_XYZ(p.v_c1, p.v_c2, p.v_c3);
        C3<float> gradv_f0_XYZ = maxwellian_df0_dv(thisVel_XYZ, p.T, p.n, p.amu, p.Z);
        return gradv_f0_XYZ;
    }
};

struct set_vx
{
    PRAGMA
    HOST DEVICE
    float operator() (float &x, CParticle &p) {
        return p.v_c1;
    }
};

struct set_vy
{
    PRAGMA
    HOST DEVICE
    float operator() (float &x, CParticle &p) {
        return p.v_c2;
    }
};

struct set_vz
{
    PRAGMA
    HOST DEVICE                
    float operator() (float &x, CParticle &p) {
        return p.v_c3;
    }
};

#endif
