#ifndef GCTERMS_HPP
#define GCTERMS_HPP

#include "c3vec.hpp"
#include "cparticle.hpp"
#include "interp.hpp"

using namespace std;

float eval_aPar ( CParticle &p, const C3Vec r, const vector<float> &r_GC, const vector<float> &bDotGradB );
float eval_vPer ( CParticle &p, const C3Vec r, const vector<float> &r_b0, const vector<C3Vec> &b0_CYL );
C3Vec eval_vGC ( CParticle &p, const C3Vec r, const float vPer, const float vPar, const vector<float> &r_b0, const vector<C3Vec> &b0_CYL, const vector<float> &r_GC, const vector<C3Vec> &curv_CYL, const vector<C3Vec> &grad_CYL );
float GetAlpComp ( const float vPer, const float phs );
float GetBetComp ( const float vPer, const float phs );

#endif
