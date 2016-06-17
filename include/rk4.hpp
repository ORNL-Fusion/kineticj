#ifndef RK4_HPP
#define RK4_HPP

#include "cparticle.hpp"
#include "c3vec.hpp"
#include "getFields.hpp"
#include "gcTerms.hpp"

#ifdef __CUDACC__
#include <thrust/device_ptr.h>
#endif

using namespace std;

C3Vec rk4_evalf ( CParticle &p, const float &t, const C3Vec &v_XYZ, const C3Vec &x, float *rVec, C3Vec *b0Vec_CYL, int nR ); 
int rk4_move_gc ( CParticle &p, const float &dt, const float &t0, const vector<float> &r_b0, const vector<C3Vec> &b0_CYL, const vector<float> &r_GC, const vector<C3Vec> &curv_CYL, const vector<C3Vec> &grad_CYL, const vector<float> &bDotGradB, const float wrf );
int rk4_move ( CParticle &p, float dt, float t0, const vector<C3Vec> &b0, const vector<C3VecI> &e1, const float wrf );
int rk4_move ( CParticle &p, const float &dt, float *r, C3Vec *b0, int n );

#include "rk4.tpp"

// Functor to wrap particle move 
//
struct moveParticle
{

    float dt;
    float *r;
    C3Vec *b;
    int n;

    moveParticle( float &_dt, float *_r, C3Vec *_b, int _n) : dt(_dt), r(_r), b(_b), n(_n) {}
#ifdef __CUDACC__
    moveParticle( float &_dt, thrust::device_ptr<float> _r, thrust::device_ptr<C3Vec> _b, int _n) : 
            dt(_dt), r(thrust::raw_pointer_cast(_r)), b(thrust::raw_pointer_cast(_b)), n(_n) {}
    __host__ __device__
#endif
    void operator() (CParticle &p) {
        rk4_move(p,dt,r,b,n);
    }

};

#endif

