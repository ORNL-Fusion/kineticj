#include "getFields.hpp"

C3Vec getB_XYZ ( CParticle &p_XYZ, const vector<float> &rVec, const vector<C3Vec> &b0Vec_CYL ) {

	float _r = sqrt ( pow(p_XYZ.c1,2) + pow(p_XYZ.c2,2) );
	float _p = atan2 ( p_XYZ.c2, p_XYZ.c1 );

	C3Vec b0_CYL, b0_XYZ;

    int status = 0;
	b0_CYL = kj_interp1D ( _r, rVec, b0Vec_CYL, status );
    p_XYZ.status = max(p_XYZ.status,status);

    b0_XYZ = rot_CYL_to_XYZ ( _p, b0_CYL, 1 );

    return b0_XYZ;
}

C3VecI getE1orB1_XYZ ( CParticle &p_XYZ, const vector<float> &rVec, const vector<C3VecI> &E1Vec_CYL, int nPhi ) {

	float _r = sqrt ( pow(p_XYZ.c1,2) + pow(p_XYZ.c2,2) );
	float _p = atan2 ( p_XYZ.c2, p_XYZ.c1 );

	C3VecI E1_CYL, E1_XYZ;

    int status = 0;
	E1_CYL = kj_interp1D ( _r, rVec, E1Vec_CYL, status );
    p_XYZ.status = max(p_XYZ.status,status);

    complex<float> ii(0,1);

    E1_XYZ = exp(ii*float(nPhi*_p)) * rot_CYL_to_XYZ ( _p, E1_CYL, 1 ) ;

    return E1_XYZ;
}
