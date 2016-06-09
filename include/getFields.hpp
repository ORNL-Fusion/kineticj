#ifndef GETFIELDS_HPP
#define GETFIELDS_HPP

#include "cparticle.hpp"
#include "c3vec.hpp"
#include "interp.hpp"
#include "rotation.hpp"

using namespace std;

C3Vec getB_XYZ ( CParticle &p_XYZ, const vector<float> &rVec, const vector<C3Vec> &b0Vec_CYL );
C3VecI getE1orB1_XYZ ( CParticle &p_XYZ, const vector<float> &rVec, const vector<C3VecI> &E1Vec_CYL, int nPhi );

#endif
