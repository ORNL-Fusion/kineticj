#ifndef GETFIELDS_HPP
#define GETFIELDS_HPP

#include "cparticle.hpp"
#include "c3vec.hpp"
#include "interp.hpp"
#include "rotation.hpp"

using namespace std;

C3Vec getB_XYZ ( CParticle &p_XYZ, float *rVec, C3Vec *b0Vec_CYL, int nR );
C3VecI getE1orB1_XYZ ( CParticle &p_XYZ, float *rVec, C3VecI *E1Vec_CYL, int nR, int nPhi );

// Functor to wrap these 

struct getPerturbedField 
{

    float *r;
    C3VecI *field_CYL;
    int nR;
    int nPhi;
    float weight;

    getPerturbedField( float *_r, C3VecI *_field_CYL, int _nR, int _nPhi, float _weight) : 
            r(_r), field_CYL(_field_CYL), nR(_nR), nPhi(_nPhi), weight(_weight) {}
    C3VecI operator() (CParticle &p) {
        // --- CHECK THIS COMMENTED PIECE ---
        //complex<float> _i(0, 1);
        //E1_XYZ = weight * exp(-_i * wrf * thisT[i]) * getE1orB1_XYZ(p, r, f, nPhi);
        C3VecI E1_XYZ = weight * getE1orB1_XYZ(p, r, field_CYL, nR, nPhi);
        C3VecI field_XYZ = E1_XYZ * (1 - p.status);
        return field_XYZ;
    }
};



#endif
