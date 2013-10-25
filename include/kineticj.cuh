#ifndef KINETICJ_H
#define KINETICJ_H

#include "cusp/complex.h"
#include "vecs.cuh"

typedef struct PARAMS_
{
  int nXGrid;
  float dv;
  int nSteps;
  int nV;
  float dtMin;
  double wrf;
  int nJp;

} params;

typedef struct GPU_MEM_
{
    complex<float> *j1xc;
    float *thisT;
    float *tJp;
    float *hanningWeight;
    float *r_kjGrid;
    C3Vec *e1Re_XYZ_kjGrid;
    C3Vec *e1Im_XYZ_kjGrid;
    CParticle_PODS *particles_XYZ_PODS;
    CParticle_PODS *particles_XYZ_0_PODS;
    float *xGrid;
    C3Vec *b0_CYL_kjGrid;
    float *df0_dv;
    C3Vec *thisOrbitE_re_XYZ;
    C3Vec *thisOrbitE_im_XYZ;
    C3Vec *thisOrbit_XYZ;
} gpu_mem;

void copyToDevice(complex<float> *j1xc, float *thisT, float *tJp, float *hanningWeight, float *r_kjGrid, C3Vec *e1Re_XYZ_kjGrid,
                  C3Vec *e1Im_XYZ_kjGrid, CParticle_PODS *particles_XYZ_PODS, CParticle_PODS *particles_XYZ_0_PODS,  float *xGrid, C3Vec *b0_CYL_kjGrid,
                  float *df0_dv, params *p, gpu_mem *gmem);

void launchKernel(params *p, gpu_mem *gmem);

void copyToHost(complex<float> *j1xc, params *p, gpu_mem *gmem);

__device__ C3Vec kj_interp1D ( float x, const float xVec[], const C3Vec yVec[], int n, int &stat );

__device__ C3Vec rk4_evalf ( CParticle_PODS &p, const float &t,
                const C3Vec &v_XYZ, const C3Vec &x,
                const C3Vec b0_CYL[],
                const float r[],
                const int n );

__device__ void rk4_move  ( CParticle_PODS &p, const float &dt, const float &t0,
                                const C3Vec b0[], const float r[], const int n);

__global__ void low_mem_kernel(complex<float> *j1xc, float *thisT, float *tJp, float *hanningWeight, float *r_kjGrid, C3Vec *e1Re_XYZ_kjGrid,
                               C3Vec *e1Im_XYZ_kjGrid, CParticle_PODS *particles_XYZ_PODS, CParticle_PODS *particles_XYZ_0_PODS,  float *xGrid, C3Vec *b0_CYL_kjGrid,
                               float *df0_dv, C3Vec *thisOrbitE_re_XYZ, C3Vec *thisOrbitE_im_XYZ, C3Vec *thisOrbit_XYZ, const __restrict params p);

#endif
