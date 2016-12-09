PRAGMA
template <typename T>
HOST DEVICE
C3<T> getB_XYZ(CParticle& p_XYZ, float *rVec, C3<T> *b0Vec_CYL, int nB)
{

#if CYLINDRICAL_INPUT_FIELDS >= 1

    float _r = sqrt(pow(p_XYZ.c1, 2) + pow(p_XYZ.c2, 2));
    float _p = atan2(p_XYZ.c2, p_XYZ.c1);

    C3<T> b0_CYL, b0_XYZ;

    int status = 0;
    b0_CYL = kj_interp1D(_r, rVec, b0Vec_CYL, nB, status);
    p_XYZ.status = max(p_XYZ.status, status);

    b0_XYZ = rot_CYL_to_XYZ(_p, b0_CYL, 1);

#else

    float _x = p_XYZ.c1;

    C3<T> b0_XYZ;

    int status = 0;
    b0_XYZ = kj_interp1D(_x, rVec, b0Vec_CYL, nB, status);
    p_XYZ.status = max(p_XYZ.status, status);

#endif

    return b0_XYZ;
}


