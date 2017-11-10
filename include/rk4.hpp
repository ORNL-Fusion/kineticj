#ifndef RK4_HPP
#define RK4_HPP

#include "cparticle.hpp"
#include "c3vec.hpp"
#include "getFields.hpp"
#include "gcTerms.hpp"

#if defined(__CUDACC__) 
#define HOST __host__ 
#define DEVICE __device__
#else
#define HOST 
#define DEVICE
#endif

using namespace std;

HOST DEVICE
C3<float> rk4_evalf ( CParticle &p, const float &t, const C3<float> &v_XYZ, const C3<float> &x, float *rVec, C3<float> *b0Vec_CYL, int nR ); 

HOST DEVICE
int rk4_move ( CParticle &p, const float &dt, float *r, C3<float> *b0, int n );

HOST DEVICE
int rk4_move_gc ( CParticle &p, const float dt, const float t0, const float *r_b0, const C3<float> *b0_CYL, int nR, 
                const float *r_GC, const C3<float> *curv_CYL, const C3<float> *grad_CYL, const float *bDotGradB, int nR_GC );

#include "rk4.tpp"

// Functor to wrap particle move 

struct moveParticle
{

    float dt;
    float *r;
    C3<float> *b;
    int n;

    moveParticle( float &_dt, float *_r, C3<float> *_b, int _n) : dt(_dt), r(_r), b(_b), n(_n) {}

    HOST DEVICE 
    void operator() (CParticle &p) {
        rk4_move(p,dt,r,b,n);
    }

};

struct moveParticle_gc
{

    float dt, t;
    float *r, *r_GC;
    C3<float> *b, *curv, *grad;
    float *bDotGrad;
    int nR, nR_GC;

    moveParticle_gc( float _dt, float _t, float *_r, C3<float> *_b, int _nR, 
                    float *_r_GC, C3<float> *_curv, C3<float> *_grad, float *_bDotGrad, int _nR_GC) 
            : dt(_dt), t(_t), r(_r), b(_b), nR(_nR), 
                r_GC(_r_GC), curv(_curv), grad(_grad), bDotGrad(_bDotGrad), nR_GC(_nR_GC) {}

    HOST DEVICE 
    void operator() (CParticle &p) {
        int status = rk4_move_gc(p,dt,t,r,b,nR,r_GC,curv,grad,bDotGrad,nR_GC);
    }

};

#endif

