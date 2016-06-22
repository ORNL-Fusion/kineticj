#ifndef GETFIELDS_HPP
#define GETFIELDS_HPP

#include "cparticle.hpp"
#include "c3vec.hpp"
#include "interp.hpp"
#include "rotation.hpp"

#ifdef __CUDACC__
#define HOST __host__ 
#define DEVICE __device__
#else
#define HOST 
#define DEVICE
#endif

#ifdef __CUDACC__
#define PRAGMA #pragma hd_warning_disable 
#else
#define PRAGMA
#endif

using namespace std;

#include "getFields.tpp"

// Functor to wrap these 

struct getPerturbedField 
{

    float *r;
    C3<std::complex<float> > *field_CYL;
    int nR;
    int nPhi;
    float weight;

    getPerturbedField( float *_r, C3<std::complex<float> > *_field_CYL, int _nR, int _nPhi, float _weight) : 
            r(_r), field_CYL(_field_CYL), nR(_nR), nPhi(_nPhi), weight(_weight) {}
    C3<std::complex<float> > operator() (CParticle &p) {
        // --- CHECK THIS COMMENTED PIECE ---
        //complex<float> _i(0, 1);
        //E1_XYZ = weight * exp(-_i * wrf * thisT[i]) * getE1orB1_XYZ(p, r, f, nPhi);
        C3<std::complex<float> > E1_XYZ = weight * getE1orB1_XYZ(p, r, field_CYL, nR, nPhi);
        C3<std::complex<float> > field_XYZ = E1_XYZ * (1 - p.status);
        return field_XYZ;
    }
};

#endif
