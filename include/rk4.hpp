#ifndef RK4_HPP
#define RK4_HPP

#include "cparticle.hpp"
#include "c3vec.hpp"
#include "getFields.hpp"
#include "gcTerms.hpp"

using namespace std;

C3Vec rk4_evalf ( CParticle &p, const float &t, const C3Vec &v, const C3Vec &x, const vector<C3Vec> &b0Vec, const vector<C3VecI> &e1, const float wrf );
C3Vec rk4_evalf ( CParticle &p, const float &t, const C3Vec &v_XYZ, const C3Vec &x, const vector<C3Vec> &b0Vec_CYL, const vector<float> &rVec ); 
int rk4_move ( CParticle &p, const float &dt, const float &t0, const vector<float> &r, const vector<C3Vec> &b0 );
int rk4_move_gc ( CParticle &p, const float &dt, const float &t0, const vector<float> &r_b0, const vector<C3Vec> &b0_CYL, const vector<float> &r_GC, const vector<C3Vec> &curv_CYL, const vector<C3Vec> &grad_CYL, const vector<float> &bDotGradB, const float wrf );
int rk4_move ( CParticle &p, float dt, float t0, const vector<C3Vec> &b0, const vector<C3VecI> &e1, const float wrf );

#endif

