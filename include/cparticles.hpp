#ifndef CPARTICLES_HPP
#define CPARTICLES_HPP

#include "cspecies.hpp"

#include "array.h"
#include "managed_allocation.h"

#ifdef __CUDACC__
#define CUDA_CALLABLE_MEMBER __host__ __device__
#else
#define CUDA_CALLABLE_MEMBER
#endif

#ifdef __CUDACC__
#include <thrust/host_vector.h>
#include <thrust/device_vector.h>
#include <thrust/copy.h>
#include <thrust/random.h>
#include <curand_kernel.h>
#endif

class CParticles: public CSpecies, public ManagedAllocation {
		public:
                std::size_t nParticles;
				sim::Array<float> c1, c2, c3, v_c1, v_c2, v_c3;
                sim::Array<int> number;
                sim::Array<float> weight;
                sim::Array<int> status;
                sim::Array<float> dvx, dvy, dvz, d3v;
                sim::Array<float> vPar, vPer, gyroPhase, u, vTh;
                sim::Array<float> vAlp, vBet, phs;
                sim::Array<float> T, n; // Temp and density of the Maxwellian at this point

                CUDA_CALLABLE_MEMBER
                CParticles(std::size_t nP) :
		            c1{nP},
                    c2{nP},
                    c3{nP},
                    v_c1{nP},
                    v_c2{nP},
                    v_c3{nP},
                    number{nP},
                    weight{nP},
                    status{nP},
                    dvx{nP},
                    dvy{nP},
                    dvz{nP},
                    d3v{nP},
                    vPar{nP},
                    vPer{nP},
                    gyroPhase{nP},
                    u{nP},
                    vTh{nP},
                    vAlp{nP},
                    vBet{nP},
                    phs{nP},
                    T{nP},
                    n{nP} 
                {};

                CUDA_CALLABLE_MEMBER
                void setParticle( 
                        int indx, float c1, float c2, float c3,
                        float v_c1, float v_c2, float v_c3,
                        int number, float weight, int status,
                        float dvx, float dvy, float dvz, float d3v,
                        float vPar, float vPer, float gyroPhase,
                        float u,
                        float vTh, float vAlp, float vBet,
                        float phs, 
                        float T,
                        float n) {

		                    this->c1[indx] = c1;
                            this->c2[indx] = c2;
                            this->c3[indx] = c3;
                            this->v_c1[indx] = v_c1;
                            this->v_c2[indx] = v_c2;
                            this->v_c3[indx] = v_c3;
                            this->number[indx] = number;
                            this->weight[indx] = weight;
                            this->status[indx] = status;
                            this->dvx[indx] = dvx;
                            this->dvy[indx] = dvy;
                            this->dvz[indx] = dvz;
                            this->d3v[indx] = d3v;
                            this->vPar[indx] = vPar;
                            this->vPer[indx] = vPer;
                            this->gyroPhase[indx] = gyroPhase;
                            this->u[indx] = u;
                            this->vTh[indx] = vTh;
                            this->vAlp[indx] = vAlp;
                            this->vBet[indx] = vBet;
                            this->phs[indx] = phs;
                            this->T[indx] = T;
                            this->n[indx] = n;

                        };



};

#endif
