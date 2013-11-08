#include "kineticj.cuh"
#include "vecs.cuh"
#include "cusp/complex.h"
#include "constants.hpp"
#include "grid_sizes.hpp"
#include "iostream"
#include "stdio.h"

// used for complex
// usinfg namespace is probably not a good idea here...lazy
using namespace cusp;

// Host wrapper functions
void copyToDevice(complex<float> *j1xc, float *thisT, float *tJp, float *hanningWeight, float *r_kjGrid, C3Vec *e1Re_XYZ_kjGrid,
                  C3Vec *e1Im_XYZ_kjGrid, CParticle_PODS *particles_XYZ_PODS, CParticle_PODS *particles_XYZ_0_PODS,  float *xGrid, C3Vec *b0_CYL_kjGrid,
                  float *df0_dv, complex<float> *all_j1xc, params *p, gpu_mem *gmem)
{
    std::cout <<"copying memory"<<std::endl;

    size_t bytes;

    bytes = sizeof(complex<float>)*p->nXGrid*p->nJp;
    cudaMalloc((void**)&(gmem->j1xc),bytes);
    cudaMemcpy(gmem->j1xc, j1xc,  bytes,  cudaMemcpyHostToDevice);

    bytes = sizeof(float)*p->nSteps;
    cudaMalloc((void**)&(gmem->thisT),bytes);
    cudaMemcpy(gmem->thisT, thisT,  bytes,  cudaMemcpyHostToDevice);

    bytes = sizeof(float)*p->nJp;
    cudaMalloc((void**)&(gmem->tJp),bytes);
    cudaMemcpy(gmem->tJp, tJp,  bytes,  cudaMemcpyHostToDevice);

    bytes = sizeof(float)*p->nSteps;
    cudaMalloc((void**)&(gmem->hanningWeight),bytes);
    cudaMemcpy(gmem->hanningWeight, hanningWeight,  bytes,  cudaMemcpyHostToDevice);

    bytes = sizeof(float)*_N_DATA;
    cudaMalloc((void**)&(gmem->r_kjGrid),bytes);
    cudaMemcpy(gmem->r_kjGrid, r_kjGrid,  bytes,  cudaMemcpyHostToDevice);

    bytes = sizeof(C3Vec)*_N_DATA;
    cudaMalloc((void**)&(gmem->e1Re_XYZ_kjGrid),bytes);
    cudaMemcpy(gmem->e1Re_XYZ_kjGrid, e1Re_XYZ_kjGrid,  bytes,  cudaMemcpyHostToDevice);

    bytes = sizeof(C3Vec)*_N_DATA;
    cudaMalloc((void**)&(gmem->e1Im_XYZ_kjGrid),bytes);
    cudaMemcpy(gmem->e1Im_XYZ_kjGrid, e1Im_XYZ_kjGrid,  bytes,  cudaMemcpyHostToDevice);

    bytes = sizeof(CParticle_PODS)*p->nV;
    cudaMalloc((void**)&(gmem->particles_XYZ_PODS),bytes);
    cudaMemcpy(gmem->particles_XYZ_PODS, particles_XYZ_PODS,  bytes,  cudaMemcpyHostToDevice);

    bytes = sizeof(CParticle_PODS)*p->nV;
    cudaMalloc((void**)&(gmem->particles_XYZ_0_PODS),bytes);
    cudaMemcpy(gmem->particles_XYZ_0_PODS, particles_XYZ_0_PODS,  bytes,  cudaMemcpyHostToDevice);

    bytes = sizeof(float)*p->nXGrid;
    cudaMalloc((void**)&(gmem->xGrid),bytes);
    cudaMemcpy(gmem->xGrid, xGrid,  bytes,  cudaMemcpyHostToDevice);

    bytes = sizeof(C3Vec)*_N_DATA;
    cudaMalloc((void**)&(gmem->b0_CYL_kjGrid),bytes);
    cudaMemcpy(gmem->b0_CYL_kjGrid, b0_CYL_kjGrid,  bytes,  cudaMemcpyHostToDevice);

    bytes = sizeof(float)*p->nV;
    cudaMalloc((void**)&(gmem->df0_dv),bytes);
    cudaMemcpy(gmem->df0_dv, df0_dv,  bytes,  cudaMemcpyHostToDevice);

    // Create following local arrays
    bytes = sizeof(complex<float>)*p->nV*p->nXGrid;
    cudaMalloc((void**)&(gmem->all_j1xc), bytes);
    cudaMemcpy(gmem->all_j1xc, all_j1xc,  bytes,  cudaMemcpyHostToDevice);

    //bytes = sizeof(C3Vec)*p->nSteps*p->nV*p->nXGrid;
    //cudaMalloc((void**)&(gmem->thisOrbitE_re_XYZ), bytes);

    //bytes = sizeof(C3Vec)*p->nSteps*p->nV*p->nXGrid;
    //cudaMalloc((void**)&(gmem->thisOrbitE_im_XYZ), bytes);

    //bytes = sizeof(C3Vec)*p->nSteps*p->nV*p->nXGrid;
    //cudaMalloc((void**)&(gmem->thisOrbit_XYZ), bytes);

    // check for error
    cudaError_t error = cudaGetLastError();
    if(error != cudaSuccess)
    {
      // print the CUDA error message and exit
      printf("CUDA error: %s\n", cudaGetErrorString(error));
      exit(-1);
    }

}

void launchKernel(params *p, gpu_mem *gmem)
{

  std::cout <<"kernel"<<std::endl;

  // Blarg
  int block_size = 32;
  int num_blocks = ceilf((p->nXGrid*p->nV)/(float)block_size);
  
  std::cout<<"nxGrid: "<<p->nXGrid<<std::endl;
  std::cout<<"num blocks: "<<num_blocks<<std::endl;

  cudaFuncSetCacheConfig(low_mem_kernel,cudaFuncCachePreferL1);

  low_mem_kernel<<<num_blocks, block_size>>>(gmem->j1xc, gmem->thisT, gmem->tJp, gmem->hanningWeight, gmem->r_kjGrid, gmem->e1Re_XYZ_kjGrid,
                               gmem->e1Im_XYZ_kjGrid, gmem->particles_XYZ_PODS, gmem->particles_XYZ_0_PODS,  gmem->xGrid, gmem->b0_CYL_kjGrid,
                               gmem->df0_dv, gmem->all_j1xc, *p);

    // check for error
    cudaDeviceSynchronize();
    cudaError_t error = cudaGetLastError();
    if(error != cudaSuccess)
    {
      // print the CUDA error message and exit
      printf("CUDA error: %s\n", cudaGetErrorString(error));
      exit(-1);
    }

}

void copyToHost(complex<float> *j1xc, complex<float> *all_j1xc, params *p, gpu_mem *gmem)
{

    std::cout <<"copying memory"<<std::endl;

    size_t bytes;
    bytes = sizeof(complex<float>)*p->nXGrid*p->nJp;
    cudaMemcpy(j1xc, gmem->j1xc,  bytes,  cudaMemcpyDeviceToHost);

    bytes = sizeof(complex<float>)*p->nXGrid*p->nV;
    cudaMemcpy(all_j1xc, gmem->all_j1xc,  bytes,  cudaMemcpyDeviceToHost);

    // check for error
    cudaError_t error = cudaGetLastError();
    if(error != cudaSuccess)
    {
      // print the CUDA error message and exit
      printf("CUDA error: %s\n", cudaGetErrorString(error));
      exit(-1);
    }

}

__device__ C3Vec kj_interp1D( float x, const float xVec[], const C3Vec yVec[], int n, int &stat )
{
	float _x, x0, x1;

	if(x<xVec[0]||x>xVec[n-1]||stat>0) {
			// Particle absorbing walls
			++stat;
			return C3Vec(0.0f,0.0f,0.0f);
	}

	_x = (x-xVec[0])/(xVec[n-1]-xVec[0])*(n-1);
	int xF = floorf(_x);
	int xC = ceilf(_x);

	x0 = floorf(_x);
	x1 = ceilf(_x);
	
	// Catch for particle at point
	if(xF==xC) {
		return yVec[xF];
	}
	else {

		C3Vec y0 = yVec[xF];
		C3Vec y1 = yVec[xC];
		C3Vec tmpA = y0+(_x-x0)*(y1-y0)/(x1-x0);
		return y0+(_x-x0)*(y1-y0)/(x1-x0);
	}
}

// Zero-order orbits
__device__ C3Vec rk4_evalf ( CParticle_PODS &p, const float &t,
                const C3Vec &v_XYZ, const C3Vec &x,
                const C3Vec b0_CYL[],
                const float r[],
                const int n ) 
 {

	// Interpolate b0 at location in CYL
	
	float _r = sqrtf ( powf(x.c1,2.0f) + powf(x.c2,2.0f) );
	float _p = atan2f ( x.c2, x.c1 );

	C3Vec thisb0_CYL, thisb0_XYZ;

	thisb0_CYL = kj_interp1D ( _r, r, b0_CYL, n, p.status );

	thisb0_XYZ = C3Vec( cosf(_p)*thisb0_CYL.c1-sinf(_p)*thisb0_CYL.c2,
					sinf(_p)*thisb0_CYL.c1+cosf(_p)*thisb0_CYL.c2,
					thisb0_CYL.c3 );

	C3Vec thisv_x_b0 ( v_XYZ.c2*thisb0_XYZ.c3-v_XYZ.c3*thisb0_XYZ.c2, 
					-1.0f*(v_XYZ.c1*thisb0_XYZ.c3-v_XYZ.c3*thisb0_XYZ.c1), 
					v_XYZ.c1*thisb0_XYZ.c2-v_XYZ.c2*thisb0_XYZ.c1);

	return thisv_x_b0*(p.q/p.m);	
}

// Zero-order orbits
__device__ void rk4_move ( CParticle_PODS &p, const float &dt, const float &t0,
                                const C3Vec b0[], const float r[], const int n) {

		C3Vec k1, k2, k3, k4, yn1, x1, x2, x3, x4, xn1; 

		C3Vec yn0(p.v_c1,p.v_c2,p.v_c3), xn0(p.c1, p.c2, p.c3);
		k1 = rk4_evalf ( p, t0 + 0.0f*dt, yn0         , xn0         , b0, r, n ) * dt;	
		x1 = yn0 * dt;
		k2 = rk4_evalf ( p, t0 + 0.5f*dt, yn0 + 0.5f*k1, xn0 + 0.5f*x1, b0, r, n ) * dt;	
		x2 = (yn0 + 0.5f*k1) * dt;
		k3 = rk4_evalf ( p, t0 + 0.5f*dt, yn0 + 0.5f*k2, xn0 + 0.5f*x2, b0, r, n ) * dt;	
		x3 = (yn0 + 0.5f*k2) * dt;
		k4 = rk4_evalf ( p, t0 + 1.0f*dt, yn0 + 1.0f*k3, xn0 + 1.0f*x3, b0, r, n ) * dt;	
		x4 = (yn0 + 1.0f*k3) * dt;

		//printf("dt: %f\n",dt);
		//printf("k1: %f, k2: %f, k3: %f, k4: %f\n", k1.c1, k2.c1, k3.c1, k4.c1);

		yn1 = yn0 + 1.0f/6.0f * (k1+2.0f*k2+2.0f*k3+k4);
		xn1 = xn0 + 1.0f/6.0f * (x1+2.0f*x2+2.0f*x3+x4);

		p.c1 = xn1.c1;
		p.c2 = xn1.c2;
		p.c3 = xn1.c3;
		p.v_c1 = yn1.c1;
		p.v_c2 = yn1.c2;
		p.v_c3 = yn1.c3;
}

// Need to find correct memory space for these. Constant cache, shared, etc...
__global__ void low_mem_kernel(complex<float> *j1xc, float *thisT, float *tJp, float *hanningWeight, float *r_kjGrid, C3Vec *e1Re_XYZ_kjGrid,
			       C3Vec *e1Im_XYZ_kjGrid, CParticle_PODS *particles_XYZ_PODS, CParticle_PODS *particles_XYZ_0_PODS,  float *xGrid, C3Vec *b0_CYL_kjGrid, 
			       float *df0_dv, complex<float> *all_j1xc, const __restrict params p)
{

    int tid = blockIdx.x * blockDim.x + threadIdx.x; 
    //int nXGrid = p.nXGrid;

    float dv = p.dv;
    int nSteps = p.nSteps;
    int nV = p.nV;
    double dtMin = p.dtMin;
    double wrf = p.wrf;
    //int nJp = p.nJp;

	int iX = tid / nV;
	int iP = tid % nV;

    //if(iX < nXGrid ) {

	//printf("tid: %i, iX: %i, iP: %i\n", tid, iX, iP);
	//printf("dv: %f, nSteps: %i, nV: %i, dtMin: %e, wrf: %f, nJp: %i\n", dv, nSteps, nV, dtMin, wrf, nJp);

    int i, istat;
    complex<float> this_j1xc;
    complex<float> f1c;

    CParticle_PODS thisParticle_XYZ;

    //for(iP=0;iP<nV;iP++) {

        thisParticle_XYZ = particles_XYZ_PODS[iP];
        thisParticle_XYZ.c1 = xGrid[iX];

        double qOverm =  thisParticle_XYZ.q/thisParticle_XYZ.m;
        float qe = thisParticle_XYZ.q;
        float h = dv * qe;
        
        // generate orbit and get time-harmonic e along it
        C3Vec e1ReTmp_XYZ, e1ImTmp_XYZ;
       
        // get Jp(t) for this spatial point
        C3VecI thisEc;
        C3VecI thisV1c;
       
        for(i=0;i<nSteps;i++) {

			e1ReTmp_XYZ = C3Vec(0,0,0);
 			e1ImTmp_XYZ = C3Vec(0,0,0);

            if(thisParticle_XYZ.status==0) {
                
                rk4_move ( thisParticle_XYZ, dtMin, thisT[i], b0_CYL_kjGrid, r_kjGrid, _N_DATA );
   
                if(thisParticle_XYZ.status==0) {
		    		istat = 0;
                    e1ReTmp_XYZ = kj_interp1D ( thisParticle_XYZ.c1, r_kjGrid, e1Re_XYZ_kjGrid, _N_DATA, istat);
		    		istat = 0;
                    e1ImTmp_XYZ = kj_interp1D ( thisParticle_XYZ.c1, r_kjGrid, e1Im_XYZ_kjGrid, _N_DATA, istat);
                }
            }           

            float tTmp = thisT[i];
            float weight = hanningWeight[i];
            float phs = -(wrf*tTmp);
            
            thisEc = C3VecI(
                            weight*complex<float>(
                                                  e1ReTmp_XYZ.c1*cosf(phs)-e1ImTmp_XYZ.c1*sinf(phs),
                                                  e1ImTmp_XYZ.c1*cosf(phs)+e1ReTmp_XYZ.c1*sinf(phs)),
                            weight*complex<float>(
                                                  e1ReTmp_XYZ.c2*cosf(phs)-e1ImTmp_XYZ.c2*sinf(phs),
                                                  e1ImTmp_XYZ.c2*cosf(phs)+e1ReTmp_XYZ.c2*sinf(phs)),
                            weight*complex<float>( 
                                                  e1ReTmp_XYZ.c3*cosf(phs)-e1ImTmp_XYZ.c3*sinf(phs),
                                                  e1ImTmp_XYZ.c3*cosf(phs)+e1ReTmp_XYZ.c3*sinf(phs))
                            );
            
            int N = nSteps - 1;
            float A = i % N;
            float B = A / N;
            int factor = ceilf(B)+1;
            
            thisV1c += -qOverm * dtMin/2 * ( factor * thisEc);
        }
      
        f1c = -thisV1c.c1*df0_dv[iP];
        
        float v0_i = particles_XYZ_0_PODS[iP].v_c1;
        
        int N = nV - 1;
        float A = iP % N;
        float B = A / N;
        int factor = ceilf(B)+1;
        
        //this_j1xc += h/2 * ( factor * v0_i*f1c);
		all_j1xc[iX*nV+iP] = h/2 * ( factor * v0_i*f1c);
    //}
    
    //[iX][0]
    //j1xc[iX*nJp] = this_j1xc;

    //}
}
