#ifndef GCTERMS_HPP
#define GCTERMS_HPP

#include "c3vec.hpp"
#include "cparticle.hpp"
#include "interp.hpp"

using namespace std;

float eval_aPar ( CParticle &p, const C3<float> r, float *r_GC, float *bDotGradB, int nGC );
float eval_vPer ( CParticle &p, const C3<float> r, float *r_b0, C3<float> *b0_CYL, int n );
C3<float> eval_vGC ( CParticle &p, const C3<float> r, const float vPer, const float vPar, float *r_b0, C3<float> *b0_CYL, int n, float *r_GC, C3<float> *curv_CYL, C3<float> *grad_CYL, int nGC );
float GetAlpComp ( const float vPer, const float phs );
float GetBetComp ( const float vPer, const float phs );

#endif
