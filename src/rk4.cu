#include "rk4.hpp"

// Zero-order orbits
HOST DEVICE
C3<float> rk4_evalf(CParticle& p, const float& t,
    const C3<float>& v_XYZ, const C3<float>& x, float *rVec, C3<float> *b0Vec_CYL, int nR)
{

    C3<float> b0_XYZ;
    b0_XYZ = getB_XYZ(p, rVec, b0Vec_CYL, nR);

    C3<float> v_x_b0 = cross(v_XYZ, b0_XYZ);

    return v_x_b0 * (p.q / p.m);
}

// Zero-order orbits
HOST DEVICE
int rk4_move(CParticle& p, const float& dt, float *r, C3<float> *b0, int nR)
{

    float t0 = 0;            

    C3<float> yn0(p.v_c1, p.v_c2, p.v_c3), xn0(p.c1, p.c2, p.c3);
    C3<float> k1, k2, k3, k4, yn1, x1, x2, x3, x4, xn1;

    k1 = rk4_evalf(p, t0 + 0.0 * dt, yn0, xn0, r, b0, nR) * dt;
    x1 = yn0 * dt;

    k2 = rk4_evalf(p, t0 + 0.5 * dt, yn0 + 0.5 * k1, xn0 + 0.5 * x1, r, b0, nR) * dt;
    x2 = (yn0 + 0.5 * k1) * dt;

    k3 = rk4_evalf(p, t0 + 0.5 * dt, yn0 + 0.5 * k2, xn0 + 0.5 * x2, r, b0, nR) * dt;
    x3 = (yn0 + 0.5 * k2) * dt;

    k4 = rk4_evalf(p, t0 + 1.0 * dt, yn0 + 1.0 * k3, xn0 + 1.0 * x3, r, b0, nR) * dt;
    x4 = (yn0 + 1.0 * k3) * dt;

    yn1 = yn0 + 1.0 / 6.0 * (k1 + 2.0 * k2 + 2.0 * k3 + k4) * (1 - p.status); // the * (1-p.status) sets the move to zero for dead particles;
    xn1 = xn0 + 1.0 / 6.0 * (x1 + 2.0 * x2 + 2.0 * x3 + x4) * (1 - p.status);

    p.c1 = xn1.c1;
    p.c2 = xn1.c2;
    p.c3 = xn1.c3;
    p.v_c1 = yn1.c1;
    p.v_c2 = yn1.c2;
    p.v_c3 = yn1.c3;

#if _PARTICLE_BOUNDARY == 1
// Particle absorbing walls
#elif _PARTICLE_BOUNDARY == 2
    // Periodic
    if (p.c1 < r[0]) {
#if DEBUGLEVEL >= 1
        cout << "Particle went left" << endl;
#endif
        p.c1 = r[nR-1] - (r[0] - p.c1);
    }
    if (p.c1 > r[nR-1]) {
#if DEBUGLEVEL >= 1
        cout << "Particle went right" << endl;
#endif
        p.c1 = r[0] + (p.c1 - r[nR-1]);
    }
#elif _PARTICLE_BOUNDARY == 3
    // Particle reflecting walls
    if (p.c1 < r[0]) {
        cout << "Particle hit the left wall" << endl;
        cout << "r.front(): " << r[0] << endl;
        p.c1 = r[0] + (r[0] - p.c1);
        p.v_c1 = -p.v_c1;
    }
    if (p.c1 > r[nR-1]) {
        cout << "Particle hit the right wall" << endl;
        cout << "r[nR-1]: " << r[nR-1] << endl;
        p.c1 = r[nR-1] - (p.c1 - r[nR-1]);
        p.v_c1 = -p.v_c1;
    }
#endif

#if DEBUGLEVEL >= 3
    cout << "\tdt: " << dt << endl;
    cout << "\tr.front(): " << r.front() << endl;
    cout << "\tr.back(): " << r.back() << endl;
    cout << "\tx0_XYZ: " << xn0.c1 << "  " << xn0.c2 << "  " << xn0.c3 << endl;
    cout << "\tv0_XYZ: " << yn0.c1 << "  " << yn0.c2 << "  " << yn0.c3 << endl;
    cout << "\tx1_XYZ: " << xn1.c1 << "  " << xn1.c2 << "  " << xn1.c3 << endl;
    cout << "\tv1_XYZ: " << yn1.c1 << "  " << yn1.c2 << "  " << yn1.c3 << endl;
    cout << "\tE: " << 0.5 * p.m * sqrt(pow(p.v_c1, 2) + pow(p.v_c2, 2) + pow(p.v_c3, 2)) / _e << endl;
#endif
    return p.status;
}

// Guiding center orbit
HOST DEVICE
int rk4_move_gc(CParticle& p, const float dt, const float t0,
    const float *r_b0, const C3<float> *b0_CYL, int nB, const float *r_GC,
    const C3<float> *curv_CYL, const C3<float> *grad_CYL,
    const float *bDotGradB, int nGC)
{

    C3<float> xn0_XYZ(p.c1, p.c2, p.c3);
    C3<float> xn0 = XYZ_to_CYL(xn0_XYZ);

    float This_vPer = eval_vPer(p, xn0, r_b0, b0_CYL, nB);
#if DEBUG_GC >= 2
    cout << "p.vPer: " << p.vPer << endl;
    cout << "p.vPar: " << p.vPar << endl;
    cout << "This_vPer: " << This_vPer << endl;
    if (isnan(p.vPer))
        exit(1);
#endif
    C3<float> This_vGC = eval_vGC(p, xn0, This_vPer, p.vPar + 0, r_b0, b0_CYL, nB, r_GC, curv_CYL, grad_CYL, nGC);
    float k1_vPar = dt * eval_aPar(p, xn0, r_GC, bDotGradB, nGC);
    C3<float> k1_vgc = dt * This_vGC;
#if DEBUG_GC >= 2
    std::cout<<"k1_vgc: "<<k1_vgc<<std::endl;
    std::cout<<"xn0: " <<xn0<<std::endl;
    cout << "Status: " << p.status << endl;
    if (isnan(k1_vgc) || isinf(k1_vgc) || isnan(xn0) || isinf(xn0) || p.status > 0) {
        p.status = 1;
        return p.status;
    }
#endif
    This_vPer = eval_vPer(p, xn0 + k1_vgc / 2.0, r_b0, b0_CYL, nB);
    This_vGC = eval_vGC(p, xn0 + k1_vgc / 2.0, This_vPer, p.vPar + k1_vPar / 2.0, r_b0, b0_CYL, nB, r_GC, curv_CYL, grad_CYL, nGC);
    float k2_vPar = dt * eval_aPar(p, xn0 + k1_vgc / 2.0, r_GC, bDotGradB, nGC);
    C3<float> k2_vgc = dt * This_vGC;
#if DEBUG_GC >= 2
    std::cout<<"k2_vgc: "<<k2_vgc<<std::endl;
    if (isnan(k2_vgc) || isinf(k2_vgc) || isnan(xn0) || isinf(xn0) || p.status > 0) {
        p.status = 1;
        return p.status;
    }
#endif
    This_vPer = eval_vPer(p, xn0 + k2_vgc / 2.0, r_b0, b0_CYL, nB);
    This_vGC = eval_vGC(p, xn0 + k2_vgc / 2.0, This_vPer, p.vPar + k2_vPar / 2.0, r_b0, b0_CYL, nB, r_GC, curv_CYL, grad_CYL, nGC);
    float k3_vPar = dt * eval_aPar(p, xn0 + k2_vgc / 2.0, r_GC, bDotGradB, nGC);
    C3<float> k3_vgc = dt * This_vGC;
#if DEBUG_GC >= 2
    std::cout<<"k3_vgc: "<<k3_vgc<<std::endl;
    if (isnan(k3_vgc) || isinf(k3_vgc) || isnan(xn0) || isinf(xn0) || p.status > 0) {
        p.status = 1;
        return p.status;
    }
#endif
    This_vPer = eval_vPer(p, xn0 + k3_vgc, r_b0, b0_CYL, nB);
    This_vGC = eval_vGC(p, xn0 + k3_vgc, This_vPer, p.vPar + k3_vPar, r_b0, b0_CYL, nB, r_GC, curv_CYL, grad_CYL, nGC);
    float k4_vPar = dt * eval_aPar(p, xn0 + k3_vgc, r_GC, bDotGradB, nGC);
    C3<float> k4_vgc = dt * This_vGC;
#if DEBUG_GC >= 2
    std::cout<<"k4_vgc: "<<k4_vgc<<std::endl;
    if (isnan(k4_vgc) || isinf(k4_vgc) || isnan(xn0) || isinf(xn0) || p.status > 0) {
        p.status = 1;
        return p.status;
    }
#endif
    float vPar1 = p.vPar + (k1_vPar + 2.0 * k2_vPar + 2.0 * k3_vPar + k4_vPar) / 6.0 * (1 - p.status);
    C3<float> xn1 = xn0 + (k1_vgc + 2.0 * k2_vgc + 2.0 * k3_vgc + k4_vgc) / 6.0 * (1 - p.status);

#if DEBUG_GC >= 1
    if (isnan(xn1) || isinf(xn1)) {
        p.status = 1;
        return p.status;
    }
#endif

    // Update particle with moved position and new vPar & vPer

    float vPer1 = eval_vPer(p, xn1, r_b0, b0_CYL, nB);

    p.vPar = vPar1;
    p.vPer = vPer1;

    C3<float> xn1_XYZ = CYL_to_XYZ(xn1);

    p.c1 = xn1_XYZ.c1;
    p.c2 = xn1_XYZ.c2;
    p.c3 = xn1_XYZ.c3;

    // Update the XYZ velocity also

    int status = 0;
    C3<float> this_b0_CYL = kj_interp1D(xn1.c1, r_b0, b0_CYL, nB, status);
    p.status = max(p.status, status);

    C3<float> this_b0_XYZ = rot_CYL_to_XYZ(xn1.c2, this_b0_CYL, 1);

    C3<float> v_abp;

    float this_wc = p.q * mag(this_b0_CYL) / p.m;
    p.phs = this_wc * t0 + p.gyroPhase;
    v_abp.c1 = GetAlpComp(vPer1, p.phs);
    v_abp.c2 = GetBetComp(vPer1, p.phs);
    v_abp.c3 = vPar1;

    p.vAlp = v_abp.c1;
    p.vBet = v_abp.c2;

    C3<float> this_v_XYZ = rot_XYZ_to_abp(v_abp, this_b0_XYZ, -1);

    p.v_c1 = this_v_XYZ.c1;
    p.v_c2 = this_v_XYZ.c2;
    p.v_c3 = this_v_XYZ.c3;

    return p.status;
}

