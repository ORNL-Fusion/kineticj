#ifndef RK4_HPP
#define RK4_HPP

#include "cparticle.hpp"
#include "c3vec.hpp"
#include "getFields.hpp"
#include "gcTerms.hpp"

#ifdef __CUDACC__
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
int rk4_move_gc ( CParticle &p, const float &dt, const float &t0, const vector<float> &r_b0, const vector<C3<float> > &b0_CYL, const vector<float> &r_GC, const vector<C3<float> > &curv_CYL, const vector<C3<float> > &grad_CYL, const vector<float> &bDotGradB, const float wrf );

//int rk4_move ( CParticle &p, float dt, float t0, const vector<C3<float>> &b0, const vector<C3<float>I> &e1, const float wrf );

HOST DEVICE
int rk4_move ( CParticle &p, const float &dt, float *r, C3<float> *b0, int n );

#include "rk4.tpp"

// Functor to wrap particle move 
//
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

#endif

