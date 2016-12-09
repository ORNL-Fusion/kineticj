#include "getFields.hpp"

#ifdef __CUDACC__
HOST DEVICE
C3<thrust::complex<float> > getE1orB1_XYZ_fromCYL(CParticle& p_XYZ, float *rVec, C3<thrust::complex<float> > *E1Vec_CYL, int nR, int nPhi)
{

    float _r = sqrt(pow(p_XYZ.c1, 2) + pow(p_XYZ.c2, 2));
    float _p = atan2(p_XYZ.c2, p_XYZ.c1);

    C3<thrust::complex<float> > E1_CYL, E1_XYZ;

    int status = 0;
    E1_CYL = kj_interp1D(_r, rVec, E1Vec_CYL, nR, status);
    p_XYZ.status = max(p_XYZ.status, status);

    thrust::complex<float> ii(0, 1);

    E1_XYZ = thrust::exp(ii * float(nPhi * _p)) * rot_CYL_to_XYZ(_p, E1_CYL, 1);

    return E1_XYZ;
}
#endif

HOST
C3<std::complex<float> > getE1orB1_XYZ_fromCYL(CParticle& p_XYZ, float *rVec, C3<std::complex<float> > *E1Vec_CYL, int nR, int nPhi)
{

    float _r = sqrt(pow(p_XYZ.c1, 2) + pow(p_XYZ.c2, 2));
    float _p = atan2(p_XYZ.c2, p_XYZ.c1);

    C3<std::complex<float> > E1_CYL, E1_XYZ;

    int status = 0;
    E1_CYL = kj_interp1D(_r, rVec, E1Vec_CYL, nR, status);
    p_XYZ.status = max(p_XYZ.status, status);

    std::complex<float> ii(0, 1);

    E1_XYZ = std::exp(ii * float(nPhi * _p)) * rot_CYL_to_XYZ(_p, E1_CYL, 1);

    return E1_XYZ;
}

#ifdef __CUDACC__
HOST DEVICE
C3<thrust::complex<float> > getE1orB1_XYZ_fromXYZ(CParticle& p_XYZ, float *rVec, C3<thrust::complex<float> > *E1Vec_CYL, int nR, float ky)
{

    float _x = p_XYZ.c1;
    float _y = p_XYZ.c2;
    float _z = p_XYZ.c3;

    C3<thrust::complex<float> > E1_XYZ;

    int status = 0;
    thrust::complex<float> ii(0, 1);

    E1_XYZ = thrust::exp(ii * float(ky*_y)) * kj_interp1D(_x, rVec, E1Vec_CYL, nR, status);

    p_XYZ.status = max(p_XYZ.status, status);

    return E1_XYZ;
}
#endif

HOST
C3<std::complex<float> > getE1orB1_XYZ_fromXYZ(CParticle& p_XYZ, float *rVec, C3<std::complex<float> > *E1Vec_CYL, int nR, float ky)
{

    float _x = p_XYZ.c1;
    float _y = p_XYZ.c2;
    float _z = p_XYZ.c3;

    C3<std::complex<float> > E1_XYZ;

    int status = 0;
    std::complex<float> ii(0, 1);

    E1_XYZ = std::exp(ii * float(ky * _y)) * kj_interp1D(_x, rVec, E1Vec_CYL, nR, status);

    p_XYZ.status = max(p_XYZ.status, status);

    return E1_XYZ;
}


