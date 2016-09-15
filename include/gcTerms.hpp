#ifndef GCTERMS_HPP
#define GCTERMS_HPP

#include "c3vec.hpp"
#include "cparticle.hpp"
#include "interp.hpp"

using namespace std;

float eval_aPar ( CParticle &p, const C3<float> r, const float *r_GC, const float *bDotGradB, int nGC );
float eval_vPer ( CParticle &p, const C3<float> r, const float *r_b0, const C3<float> *b0_CYL, int n );
C3<float> eval_vGC ( CParticle &p, const C3<float> r, const float vPer, const float vPar, const float *r_b0, const C3<float> *b0_CYL, int n, const float *r_GC, const C3<float> *curv_CYL, const C3<float> *grad_CYL, int nGC );
float GetAlpComp ( const float vPer, const float phs );
float GetBetComp ( const float vPer, const float phs );

#endif
