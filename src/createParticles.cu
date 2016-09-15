#include "createParticles.hpp"
#include "rk4.hpp"
#include <cmath>

float GetGyroPhase(const C3<float> v_abp)
{

    // alp is mostly in the x / r direction
    // bet is mostly z direction

    float alp = v_abp.c1;
    float bet = v_abp.c2;

    return atan2(alp, bet);
}

PRAGMA
HOST DEVICE
float maxwellian(float vx, float vy, float vz, float vTh)
{

    float weight_x = 1.0 / (vTh * std::sqrt(physConstants::pi)) * std::exp(-std::pow(vx, 2) / std::pow(vTh, 2));
    float weight_y = 1.0 / (vTh * std::sqrt(physConstants::pi)) * std::exp(-std::pow(vy, 2) / std::pow(vTh, 2));
    float weight_z = 1.0 / (vTh * std::sqrt(physConstants::pi)) * std::exp(-std::pow(vz, 2) / std::pow(vTh, 2));

    return weight_x * weight_y * weight_z;
}

PRAGMA
HOST DEVICE
float get_vTh(const float _amu, const float _Z, const float _T_keV)
{

    float m = _amu * physConstants::amu;
    float kT_joule = _T_keV * 1e3 * physConstants::e; // This may actually be E_keV so may need a 3/2 somewhere
    float vTh = std::sqrt(2.0 * kT_joule / m);

    return vTh;
}

PRAGMA
HOST DEVICE
C3<float> maxwellian_df0_dv(const C3<float> _v, const float _T_keV, const float _n_m3, const float _amu, const float _Z)
{

    C3<float> df0_dv;

    float vTh = get_vTh(_amu, _Z, _T_keV);

    float _vx = _v.c1;
    float _vy = _v.c2;
    float _vz = _v.c3;

    // Get the 3 components of df0_dv at this point in velocity space

    float h = vTh / 1000.0;
    float vxL = _vx - h;
    float vxR = _vx + h;
    float fL = maxwellian(vxL, _vy, _vz, vTh);
    float fR = maxwellian(vxR, _vy, _vz, vTh);
    float _df0_dv = (-fL + fR) / (2 * h);

    df0_dv.c1 = _df0_dv * _n_m3;

    float vyL = _vy - h;
    float vyR = _vy + h;
    fL = maxwellian(_vx, vyL, _vz, vTh);
    fR = maxwellian(_vx, vyR, _vz, vTh);
    _df0_dv = (-fL + fR) / (2 * h);

    df0_dv.c2 = _df0_dv * _n_m3;

    float vzL = _vz - h;
    float vzR = _vz + h;
    fL = maxwellian(_vx, _vy, vzL, vTh);
    fR = maxwellian(_vx, _vy, vzR, vTh);
    _df0_dv = (-fL + fR) / (2 * h);

    df0_dv.c3 = _df0_dv * _n_m3;

    return df0_dv;
}

vector<CParticle> create_particles(float x, float amu, float Z, float T_keV, float n_m3,
    int nPx, int nPy, int nPz, int nThermal, float& dv, float *r, C3<float> *b0_CYL, int nR)
{

    vector<CParticle> pList;

    int nP = nPx * nPy * nPz;
    pList.resize(nP);

    float vTh = get_vTh(amu, Z, T_keV);

#if DEBUG_MAXWELLIAN >= 1
    cout << "amu: " << amu << endl;
    cout << "Z: " << Z << endl;
    cout << "vTh: " << vTh << endl;
#endif

    float vxRange = vTh * nThermal * 2;
    float vxMin = -vxRange / 2.0;
    float dvx = vxRange / (nPx - 1);

    float vyRange = vTh * nThermal * 2;
    float vyMin = -vyRange / 2.0;
    float dvy = vyRange / (nPy - 1);

    float vzRange = vTh * nThermal * 2;
    float vzMin = -vzRange / 2.0;
    float dvz = vzRange / (nPz - 1);

    dv = dvx * dvy * dvz; // Return the Jacobian (volume element for integration later)

    float TestIntegratedValue = 0;

    int cnt = 0;
    for (int i = 0; i < nPx; i++) {
        for (int j = 0; j < nPy; j++) {
            for (int k = 0; k < nPz; k++) {

                float thisvx = vxMin + i * dvx;
                float thisvy = vyMin + j * dvy;
                float thisvz = vzMin + k * dvz;

                float weight = maxwellian(thisvx, thisvy, thisvz, vTh) * n_m3;

                TestIntegratedValue += weight * dv;

                CParticle p(x, 0.0, 0.0, thisvx, thisvy, thisvz, amu, Z, weight, T_keV, n_m3);
                pList[cnt] = p;
                pList[cnt].number = cnt;
                pList[cnt].vTh = vTh;

                pList[cnt].d3v = dv;

                // Get vPar, vPer and mu for guiding center integration

                C3<float> thisV_XYZ(thisvx, thisvy, thisvz);
                int iStat = 0;
                C3<float> this_b0_CYL = kj_interp1D(x, r, b0_CYL, nR, iStat);
                if(iStat>0) {
                    cout << "ERROR : Interpolation failure on b0_CYL" << endl;
                    exit(1);
                }
                C3<float> this_b0_XYZ = rot_CYL_to_XYZ(0, this_b0_CYL, 1);
                float bMag = mag(this_b0_XYZ);
                float vMag = mag(thisV_XYZ);

                C3<float> thisV_abp = rot_XYZ_to_abp(thisV_XYZ, this_b0_XYZ, 0);

                pList[cnt].vPar = thisV_abp.c3;
                //std::cout << "vPar 1: " << pList[cnt].vPar << std::endl;
                //std::cout << "vPar 2: " << dot(thisV_XYZ,this_b0_XYZ/mag(this_b0_XYZ)) << std::endl;
                pList[cnt].vPer = std::sqrt(std::pow(thisV_abp.c1, 2) + std::pow(thisV_abp.c2, 2));
                //std::cout << "vPer 1: " << pList[cnt].vPer << std::endl;
                //std::cout << "vPer 2: " << std::sqrt(std::pow(mag(thisV_XYZ),2)-std::pow(pList[cnt].vPar,2)) << std::endl;
 
                pList[cnt].gyroPhase = GetGyroPhase(thisV_abp);
                pList[cnt].u = pList[cnt].m * std::pow(pList[cnt].vPer, 2) / (2.0 * bMag);

#if GC_ORBITS > 0
                // Update the starting point to be at the guiding center

                int nTGC = 40;
                float wc = std::abs(pList[cnt].q * bMag / pList[cnt].m);
                float dTGC = 2*physConstants::pi/wc/nTGC;
                CParticle pGC(pList[cnt]);
                float averageX=0;
                float averageY=0;
                float averageZ=0;
                int MoveStatus = 0;
                for(int iGC=0; iGC<nTGC; iGC++) {
                    MoveStatus = rk4_move(pGC, dTGC, r, b0_CYL, nR);

                    averageX += pGC.c1;
                    averageY += pGC.c2;
                    averageZ += pGC.c3;
                }
                averageX = averageX/nTGC;
                averageY = averageY/nTGC;
                averageZ = averageZ/nTGC;

                pList[cnt].c1 = averageX;
                pList[cnt].c2 = averageY;
                pList[cnt].c3 = averageZ;
#endif

#if DEBUG_MAXWELLIAN >= 2
                cout << "ThisVx: " << thisvx << endl;
                cout << "ThisVy: " << thisvy << endl;
                cout << "ThisVz: " << thisvz << endl;
                cout << "b0_XYZ: " << this_b0_XYZ.c1 << ", " << this_b0_XYZ.c2 << ", " << this_b0_XYZ.c3 << endl;
                cout << "vMag: " << vMag << endl;
                cout << "vPer: " << pList[cnt].vPer << endl;
                cout << "vPar: " << pList[cnt].vPar << endl;
                cout << "u: " << pList[cnt].u << endl
                     << endl;
                if (isnan(pList[cnt].u))
                    exit(1);
                if (vMag > 3e8)
                    exit(1);
#endif
                cnt++;
            }
        }
    }

#if DEBUG_MAXWELLIAN >= 1
    cout << "TestIntegratedValue: " << TestIntegratedValue << endl;
#endif
    return pList;
}
