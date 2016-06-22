PRAGMA
template <typename T>
HOST DEVICE
C3<T> getB_XYZ(CParticle& p_XYZ, float *rVec, C3<T> *b0Vec_CYL, int nB)
{

    float _r = sqrt(pow(p_XYZ.c1, 2) + pow(p_XYZ.c2, 2));
    float _p = atan2(p_XYZ.c2, p_XYZ.c1);

    C3<T> b0_CYL, b0_XYZ;

    int status = 0;
    b0_CYL = kj_interp1D(_r, rVec, b0Vec_CYL, nB, status);
    p_XYZ.status = max(p_XYZ.status, status);

    b0_XYZ = rot_CYL_to_XYZ(_p, b0_CYL, 1);

    return b0_XYZ;
}

PRAGMA
template <typename T>
HOST DEVICE
C3<T> getE1orB1_XYZ(CParticle& p_XYZ, float *rVec, C3<T> *E1Vec_CYL, int nR, int nPhi)
{

    float _r = sqrt(pow(p_XYZ.c1, 2) + pow(p_XYZ.c2, 2));
    float _p = atan2(p_XYZ.c2, p_XYZ.c1);

    C3<T> E1_CYL, E1_XYZ;

    int status = 0;
    E1_CYL = kj_interp1D(_r, rVec, E1Vec_CYL, nR, status);
    p_XYZ.status = max(p_XYZ.status, status);

    complex<float> ii(0, 1);

    E1_XYZ = exp(ii * float(nPhi * _p)) * rot_CYL_to_XYZ(_p, E1_CYL, 1);

    return E1_XYZ;
}
