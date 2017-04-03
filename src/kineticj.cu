#include "c3vec.hpp"
#include "constants.hpp"
#include "cparticle.hpp"
#include "createParticles.hpp"
#include "cspecies.hpp"
#include "interp.hpp"
#include "read_e_field.hpp"
#include "read_gc_file.hpp"
#include "rk4.hpp"
#include "rotation.hpp"
#include <algorithm>
#include <assert.h>
#include <complex>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <libconfig.h++>
#include <netcdf>
#include <new> // for stl::bad_alloc
#include <omp.h>
#include <string>
#include <sys/stat.h>
#include <sys/types.h>
#include <unistd.h>
#include <vector>
#include <numeric>

#ifdef __CUDACC__
#include <thrust/device_vector.h>
#include <thrust/host_vector.h>
#include <thrust/complex.h>
#endif

#if CLOCK >= 1
#include <ctime>
#endif

#if USEPAPI >= 1
#include <papi.h>
#endif

#if LOWMEM_USEPAPI >= 1
#include <papi.h>
#endif

//#include <google/profiler.h>

#ifdef __CUDA_ARCH__
#define PRINT cuPrintf
#else
#define PRINT printf
#endif

using namespace netCDF;
using namespace exceptions;

// Calculate the jP given some know E and f(v)

int main(int argc, char** argv)
{

#ifdef __CUDACC__

    int num_gpus = 0;   // number of CUDA GPUs

    printf("%s Starting...\n\n", argv[0]);

    // determine the number of CUDA capable GPUs
    cudaGetDeviceCount(&num_gpus);

    if (num_gpus < 1)
    {
        printf("no CUDA capable devices were detected\n");
        return 1;
    }

    // display CPU and GPU configuration
    printf("number of host CPUs:\t%d\n", omp_get_num_procs());
    printf("number of CUDA devices:\t%d\n", num_gpus);

    for (int i = 0; i < num_gpus; i++)
    {
        cudaDeviceProp dprop;
        cudaGetDeviceProperties(&dprop, i);
        printf("   %d: %s\n", i, dprop.name);
    }

#endif

    // Make sure the "output/" directory exists

    stringstream outputDirName;
    outputDirName << "output/";

    // check directory exists
    struct stat st;
    int dirTest = stat(outputDirName.str().c_str(), &st);
    if (dirTest != 0) {
        std::cout << "Had to create output/ directory" << std::endl;
        int mkDirStat = mkdir(outputDirName.str().c_str(),
            S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
    }

#if CLOCK >= 1
    clock_t ProgramTime = clock();
#endif

#if (USEPAPI >= 1 || LOWMEM_USEPAPI >= 1)
    float realTime0, cpuTime0, realTime = 0, cpuTime = 0, mFlops = 0;
    long long flpIns0, flpIns = 0;
    int papiReturn;

    cpuTime0 = cpuTime;
    realTime0 = realTime;
    flpIns0 = flpIns;
    papiReturn = PAPI_flops(&realTime, &cpuTime, &flpIns, &mFlops);
    if (papiReturn < 0) {
        std::cout << "ERROR: PAPI Failed to initialize with error code: " << papiReturn
             << std::endl;
        std::cout << "ERROR: See papi.h for error code explanations " << std::endl;
        exit(1);
    }
    printf("Real_time:\t%f\nProc_time:\t%f\nTotal flpins:\t%lld\nMFLOPS:\t\t%f\n",
        realTime, cpuTime, flpIns, mFlops);

    papiReturn = PAPI_flops(&realTime, &cpuTime, &flpIns, &mFlops);
    if (papiReturn < 0) {
        std::cout << "ERROR: PAPI Failed to initialize with error code: " << papiReturn
             << std::endl;
        std::cout << "ERROR: See papi.h for error code explanations " << std::endl;
        exit(1);
    } else {
        std::cout << "PAPI called successfully with return code: " << papiReturn << std::endl;
    }
    printf("Real_time:\t%f\nProc_time:\t%f\nTotal flpins:\t%lld\nMFLOPS:\t\t%f\n",
        realTime, cpuTime, flpIns, mFlops);
#endif

    // Read config file

    libconfig::Config cfg;
    string cfgName = "kj.cfg";

    try {
        cfg.readFile(cfgName.c_str());
    } catch (const libconfig::FileIOException& fioex) {
        std::cerr << "I/O error while reading file." << std::endl;
        return (EXIT_FAILURE);
    } catch (const libconfig::ParseException& pex) {
        std::cerr << "Parse error at " << pex.getFile() << ":" << pex.getLine()
                  << " - " << pex.getError() << std::endl;
        return (EXIT_FAILURE);
    }

    int species_number = cfg.lookup("species_number");
    float T_keV_cfg = cfg.lookup("T_keV");

    // Read E
    string input_fName = cfg.lookup("input_fName");
    vector<C3<std::complex<float> > > e1_CYL, b1_CYL;
    vector<C3<float> > b0_CYL, b0_XYZ;
    vector<float> r, n_m3;
    float freq;
    int eReadStat = read_e_field(input_fName, species_number, freq, r, n_m3, 
                    e1_CYL, b1_CYL, b0_CYL);
#if GC_ORBITS >= 1

    // Read GC terms
    std::string gc_fName;
    if(cfg.lookupValue("gc_fName",gc_fName)) {
    } else {
        gc_fName = "gc_terms.nc";
    }
    //string gc_fName = cfg.lookup("gc_fName");

    vector<C3<float> > curv_CYL, grad_CYL;
    std::vector<float> r_gc, bDotGradB;
    int gcReadStat = read_gc_file(gc_fName, r_gc, curv_CYL, grad_CYL, bDotGradB);

#endif

    float wrf = freq * 2 * physConstants::pi;
    float xGridMin = cfg.lookup("xGridMin");
    float xGridMax = cfg.lookup("xGridMax");
    int nXGrid = cfg.lookup("nXGrid");
    std::cout << "nXGrid: " << nXGrid << std::endl;

    vector<float> xGrid(nXGrid);
    vector<float> density_m3(nXGrid);
    vector<float> T_keV(nXGrid);
    vector<float> wrf_wc(nXGrid);
    vector<float> bMag_kjGrid(nXGrid);

    float xGridRng = 0;
    float xGridStep = 0;

    if (nXGrid > 1) {
        xGridRng = xGridMax - xGridMin;
        xGridStep = xGridRng / (nXGrid - 1);
    }

    for (int iX = 0; iX < nXGrid; iX++) {
        xGrid[iX] = xGridMin + iX * xGridStep;
        int iStat=0;
        density_m3[iX] = kj_interp1D(xGrid[iX], &r[0], &n_m3[0], r.size(), iStat);
        if(iStat>0) {
            std::cout << "INTERPOLATION ERROR for Density" << std::endl;
            exit(1);
        }
        iStat=0;
        C3<float> this_b0 = kj_interp1D(xGrid[iX], &r[0], &b0_CYL[0], r.size(), iStat);
        if(iStat>0) {
            std::cout << "INTERPOLATION ERROR for B0" << std::endl;
            exit(1);
        }
        bMag_kjGrid[iX] = mag(this_b0);
        T_keV[iX] = T_keV_cfg; // kj_interp1D(xGrid[iX],r,n_m3);
    }

    float MaxB0 = *max_element(bMag_kjGrid.begin(), bMag_kjGrid.end());

#if USEPAPI >= 1
    cpuTime0 = cpuTime;
    realTime0 = realTime;
    flpIns0 = flpIns;
    papiReturn = PAPI_flops(&realTime, &cpuTime, &flpIns, &mFlops);
    printf("\nStartup performance:\n");
    printf("Real_time:\t%f\nProc_time:\t%f\nTotal flpins:\t%lld\nMFLOPS:\t\t%f\n",
        realTime - realTime0, cpuTime - cpuTime0, flpIns - flpIns0, mFlops);
#endif

    float nRFCycles = cfg.lookup("nRFCycles");
    float nStepsPerCyclotronPeriod = cfg.lookup("nStepsPerCyclotronPeriod");
    float tRF = (2 * physConstants::pi) / wrf;
    int nPhi = cfg.lookup("nPhi");
    float ky = cfg.lookup("ky"); // Only used for -DCYLINDRICAL_INPUT_FIELDS=0
    float kz = cfg.lookup("kz"); // Only used for -DCYLINDRICAL_INPUT_FIELDS=0
    int istat = 0;
    int nPx = cfg.lookup("nP_Vx");
    int nPy = cfg.lookup("nP_Vy");
    int nPz = cfg.lookup("nP_Vz");
    float amu = cfg.lookup("species_amu");
    float Z = cfg.lookup("species_Z");
    int nThermal = cfg.lookup("nThermal");
    long int nP = nPx * nPy * nPz;
    float wc = std::abs ( Z * physConstants::e * MaxB0 / (amu * physConstants::mi) );
    float cyclotronPeriod = 2 * physConstants::pi / wc;
    float dtMin = -cyclotronPeriod / nStepsPerCyclotronPeriod;

    int SanityCheck = 0;

    if (std::isinf(cyclotronPeriod)) ++SanityCheck;

    if (SanityCheck > 0) {
        std::cout<<"SanityCheck Failure"<<std::endl;
        exit(SanityCheck);
    }

    int nSteps = nRFCycles * tRF / std::abs(dtMin) + 1;

    for (int iX = 0; iX < nXGrid; iX++) {
        float this_wc = Z * physConstants::e * bMag_kjGrid[iX] / (amu * physConstants::amu);
        wrf_wc[iX] = wrf / this_wc;
        std::cout<<"mass: "<<amu*physConstants::mi<<std::endl;
        std::cout<<"wrf_wc[iX]: "<<wrf_wc[iX]<<std::endl;
    }

#if PRINT_INFO >= 1
    std::cout << "dtMin [s]: " << dtMin << std::endl;
    std::cout << "Cyclotron Period: " << cyclotronPeriod << std::endl;
    std::cout << "RF Period: " << tRF << std::endl;
    std::cout << "nSteps: " << nSteps << std::endl;
    std::cout << "nStepsPerCyclotronPeriod: " << nStepsPerCyclotronPeriod << std::endl;
    std::cout << "freq: " << freq << std::endl;
    std::cout << "Max B0: " << MaxB0 << std::endl;
#endif

    vector<float> thisT;
    try {
        thisT.resize(nSteps);
    } catch (const std::bad_alloc& error) {
        std::cout << "Allocation error at " << __FILE__ << __LINE__ << std::endl;
        std::cout << error.what();
    }

    for (int i = 0; i < nSteps; i++) {
        thisT[i] = i * dtMin; //+1.5*dtMin;
    }

    vector<float> hanningWeight(nSteps);
    vector<float> expWeight(nSteps);
    vector<float> linearWeight(nSteps);
    for (int i = 0; i < nSteps; i++) {
        // linearWeight[i]=thisT[i]*1.0/(tRF*nRFCycles)+1.0;
        hanningWeight[i] = 0.5 * (1 - cos(2 * physConstants::pi * i / (nSteps - 1))); // Regular
        // hanningWeight[i]=0.5*(1-cos(2*physConstants::pi*i/(nSteps*0.25-1))); //Sharper
        // hanningWeight[i] = linearWeight[i];
        if (i < nSteps / 2)
            hanningWeight[i] = 1; // Regular
        // if(i<nSteps*7.0/8.0) hanningWeight[i]=1; //Sharper
        // complex<float> _i (0.0,1.0);
        // complex<float> wrf_c (wrf,wrf*0.0025);
        // expWeight[i] = 1.0;//std::abs(exp(-_i*wrf_c*thisT[i]));
        // hanningWeight[i] = hanningWeight[i] * expWeight[i];
    }

    vector<vector<float> > j1x(nXGrid), j1y(nXGrid), j1z(nXGrid);
    vector<complex<float> > j1xc(nXGrid), j1yc(nXGrid), j1zc(nXGrid);

#if defined(_OPENMP)
    int nThreads, tid, spoken = 0;
#endif

#if CLOCK >= 1
        clock_t startTimeFunctor = clock();
#endif

    float dv;

    // Create worklist of nX * nP particles

    long int nWork = nXGrid * nP;

    vector<CParticle> particleWorkList;
    for (int iX = 0; iX < nXGrid; iX++) {

        vector<CParticle> moreWork(
            create_particles(xGrid[iX], amu, Z, T_keV[iX], density_m3[iX], nPx, nPy,
                nPz, nThermal, dv, &r[0], &b0_CYL[0], r.size() ));

        particleWorkList.insert( particleWorkList.end(), moreWork.begin(), moreWork.end() );
    }

#ifdef __CUDACC__
    std::cout<<"Copying particle worklist to device ..."<<std::endl;
    thrust::device_vector<CParticle> particleWorkList_device = particleWorkList;

    thrust::device_vector<float> vx_device(nWork,0);
    thrust::device_vector<float> vy_device(nWork,0);
    thrust::device_vector<float> vz_device(nWork,0);

    thrust::transform( vx_device.begin(), vx_device.end(), particleWorkList_device.begin(), vx_device.begin(), set_vx() );
    thrust::transform( vy_device.begin(), vy_device.end(), particleWorkList_device.begin(), vy_device.begin(), set_vy() );
    thrust::transform( vz_device.begin(), vz_device.end(), particleWorkList_device.begin(), vz_device.begin(), set_vz() );

    std::cout<<"DONE"<<std::endl;
#endif

    vector<float> vx(nWork,0);
    vector<float> vy(nWork,0);
    vector<float> vz(nWork,0);

    transform( vx.begin(), vx.end(), particleWorkList.begin(), vx.begin(), set_vx() );
    transform( vy.begin(), vy.end(), particleWorkList.begin(), vy.begin(), set_vy() );
    transform( vz.begin(), vz.end(), particleWorkList.begin(), vz.begin(), set_vz() );

    // Velocity space calculation

    vector<C3<float> > df0_dv_XYZ(nWork,0);
    vector<C3<std::complex<float> > > E1(nWork,0);
    vector<C3<std::complex<float> > > B1(nWork,0);
    vector<C3<std::complex<float> > > vCrossB(nWork,0);
    vector<C3<std::complex<float> > > vCrossB_E1(nWork,0);
    vector<complex<float> > forceDotGradf0(nWork,0);
    vector<complex<float> > dtIntegral(nWork,0);
    vector<complex<float> > f1(nWork,0);
    vector<complex<float> > vxf1(nWork,0);
    vector<complex<float> > vyf1(nWork,0);
    vector<complex<float> > vzf1(nWork,0);

#ifdef __CUDACC__

    thrust::device_vector<C3<float> > df0_dv_XYZ_device(nWork,0);
    thrust::device_vector<C3<thrust::complex<float> > > E1_device(nWork,0);
    thrust::device_vector<C3<thrust::complex<float> > > B1_device(nWork,0);
    thrust::device_vector<C3<thrust::complex<float> > > vCrossB_device(nWork,0);
    thrust::device_vector<C3<thrust::complex<float> > > vCrossB_E1_device(nWork,0);
    thrust::device_vector<thrust::complex<float> > forceDotGradf0_device(nWork,0);
    thrust::device_vector<thrust::complex<float> > dtIntegral_device(nWork,0);
    thrust::device_vector<thrust::complex<float> > f1_device(nWork,0);
    thrust::device_vector<thrust::complex<float> > vxf1_device(nWork,0);
    thrust::device_vector<thrust::complex<float> > vyf1_device(nWork,0);
    thrust::device_vector<thrust::complex<float> > vzf1_device(nWork,0);

    thrust::host_vector<CParticle> p_host(nWork);
    thrust::host_vector<C3<float> > df0_dv_XYZ_host(nWork,0);
    thrust::host_vector<C3<thrust::complex<float> > > E1_host(nWork,0);
    thrust::host_vector<C3<thrust::complex<float> > > B1_host(nWork,0);
    thrust::host_vector<C3<thrust::complex<float> > > vCrossB_host(nWork,0);
    thrust::host_vector<C3<thrust::complex<float> > > vCrossB_E1_host(nWork,0);
    thrust::host_vector<thrust::complex<float> > forceDotGradf0_host(nWork,0);
    thrust::host_vector<thrust::complex<float> > dtIntegral_host(nWork,0);
    thrust::host_vector<thrust::complex<float> > f1_host(nWork,0);
    thrust::host_vector<thrust::complex<float> > vxf1_host(nWork,0);
    thrust::host_vector<thrust::complex<float> > vyf1_host(nWork,0);
    thrust::host_vector<thrust::complex<float> > vzf1_host(nWork,0);
 
    // Also copy across the fields to be interpolated to the device
    // For some reason, and only for the thrust::complex<> type, it
    // seems required to go via a thrust::host_vector<> rather than
    // directly to the thrust::device_vector<> as with the floats.

    thrust::host_vector<C3<thrust::complex<float> > > e1_CYL_host(e1_CYL);
    thrust::host_vector<C3<thrust::complex<float> > > b1_CYL_host(b1_CYL);

    thrust::device_vector<float> r_device = r;
    thrust::device_vector<C3<float> > b0_CYL_device = b0_CYL;
    thrust::device_vector<C3<thrust::complex<float> > > e1_CYL_device = e1_CYL_host;
    thrust::device_vector<C3<thrust::complex<float> > > b1_CYL_device = b1_CYL_host;

    float *r_dPtr_raw = thrust::raw_pointer_cast(r_device.data());
    C3<float> *b0_dPtr_raw = thrust::raw_pointer_cast(b0_CYL_device.data());
    C3<thrust::complex<float> > *e1_dPtr_raw = thrust::raw_pointer_cast(e1_CYL_device.data());
    C3<thrust::complex<float> > *b1_dPtr_raw = thrust::raw_pointer_cast(b1_CYL_device.data());

#endif

    // Move particles
    std::cout << "Moving particles with for_each ..." << std::endl;

    for (int i = 0; i < nSteps; i++) {

        //std::cout<<"Step "<<i<<" of "<<nSteps<<std::endl;

        float dtIntFac = 1;
        if (i > 0) dtIntFac = 2;

        dtIntFac = dtMin / 2.0 * dtIntFac;

#ifdef __CUDACC__
        // Move particle
        thrust::for_each( particleWorkList_device.begin(), particleWorkList_device.end(), 
                        moveParticle(dtMin, r_dPtr_raw, b0_dPtr_raw, r.size()) ); 
        thrust::copy(particleWorkList_device.begin(),particleWorkList_device.end(),p_host.begin());

        // df0(v)/dv 
        thrust::transform( particleWorkList_device.begin(), particleWorkList_device.end(), df0_dv_XYZ_device.begin(), 
                        get_df0_dv() ); 
        thrust::copy(df0_dv_XYZ_device.begin(),df0_dv_XYZ_device.end(),df0_dv_XYZ_host.begin());

        // E1(x) 
        thrust::transform( particleWorkList_device.begin(), particleWorkList_device.end(), E1_device.begin(), 
                        getPerturbedField_device(r_dPtr_raw,e1_dPtr_raw,r.size(),nPhi,ky,kz,hanningWeight[i],wrf,thisT[i]) ); 
        thrust::copy(E1_device.begin(),E1_device.end(),E1_host.begin());

        // B1(x) 
        thrust::transform( particleWorkList_device.begin(), particleWorkList_device.end(), B1_device.begin(), 
                        getPerturbedField_device(r_dPtr_raw,b1_dPtr_raw,r.size(),nPhi,ky,kz,hanningWeight[i],wrf,thisT[i]) ); 
        thrust::copy(B1_device.begin(),B1_device.end(),B1_host.begin());

        // v x B1 
        thrust::transform( particleWorkList_device.begin(), particleWorkList_device.end(), B1_device.begin(), vCrossB_device.begin(), 
                        vCross<thrust::complex<float> >() );
        thrust::copy(vCrossB_device.begin(),vCrossB_device.end(),vCrossB_host.begin());

        // E1 + v x B1
        thrust::transform( E1_device.begin(), E1_device.end(), vCrossB_device.begin(), vCrossB_E1_device.begin(), 
                        thrust::plus<C3<thrust::complex<float> > >() );
        thrust::copy(vCrossB_E1_device.begin(),vCrossB_E1_device.end(),vCrossB_E1_host.begin());

        //  (E1 + v x B1) . grad_v(f0(v))
        thrust::transform( vCrossB_E1_device.begin(), vCrossB_E1_device.end(), df0_dv_XYZ_device.begin(), forceDotGradf0_device.begin(), 
                        doDotProduct_device() );
        thrust::copy(forceDotGradf0_device.begin(),forceDotGradf0_device.end(),forceDotGradf0_host.begin());

        // int( (E1 + v x B1) . grad_v(f0(v)), dt ) via running dt integral
        thrust::transform( dtIntegral_device.begin(), dtIntegral_device.end(), forceDotGradf0_device.begin(), dtIntegral_device.begin(), 
                        runningIntegral<thrust::complex<float> >(dtIntFac) );
        thrust::copy(dtIntegral_device.begin(),dtIntegral_device.end(),dtIntegral_host.begin());

        // f1(v) = -q/m * int( (E1 + v x B1) . grad_v(f0(v)), dt )
        thrust::transform( dtIntegral_device.begin(), dtIntegral_device.end(), particleWorkList_device.begin(), f1_device.begin(), 
                        multiplyByChargeOverMass<thrust::complex<float> >() ); 
        thrust::copy(f1_device.begin(),f1_device.end(),f1_host.begin());

        // q . f1(v) // first step in velocity momemnt for current calculation 
        thrust::transform( f1_device.begin(), f1_device.end(), particleWorkList_device.begin(), f1_device.begin(), 
                        multiplyByCharge<thrust::complex<float> >() ); 
        // q . vx . f1(v) 
        thrust::transform( f1_device.begin(), f1_device.end(), vx_device.begin(), vxf1_device.begin(), 
                        thrust::multiplies<thrust::complex<float> >() ); 

        thrust::copy(vxf1_device.begin(),vxf1_device.end(),vxf1_host.begin());

        // q . vy . f1(v) 
        thrust::transform( f1_device.begin(), f1_device.end(), vy_device.begin(), vyf1_device.begin(), 
                        thrust::multiplies<thrust::complex<float> >() ); 
        // q . vz . f1(v) 
        thrust::transform( f1_device.begin(), f1_device.end(), vz_device.begin(), vzf1_device.begin(), 
                        thrust::multiplies<thrust::complex<float> >() ); 


#endif 

#if DO_CPU_ITERATOR_APPROACH > 0

        // Move particle
#if GC_ORBITS >= 1
        for_each( particleWorkList.begin(), particleWorkList.end(), 
                        moveParticle_gc(dtMin, thisT[i], &r[0], &b0_CYL[0], r.size(), &r_gc[0], &curv_CYL[0], &grad_CYL[0], &bDotGradB[0], r_gc.size() ) ); 
#else
        for_each( particleWorkList.begin(), particleWorkList.end(), 
                        moveParticle(dtMin, &r[0], &b0_CYL[0], r.size() ) ); 
#endif
#ifdef __CUDACC__
        std::cout<<"move CPU: "<<particleWorkList[0].c1<<particleWorkList[0].c2<<particleWorkList[0].c3
                <<" GPU: "<<p_host[0].c1<<p_host[0].c2<<p_host[0].c3<<std::endl;
#endif
        // df0(v)/dv 
        transform( particleWorkList.begin(), particleWorkList.end(), df0_dv_XYZ.begin(), 
                        get_df0_dv() ); 

#ifdef __CUDACC__
       std::cout<<"df0_dv_XYZ CPU: "<<df0_dv_XYZ[0]<<" GPU: "<<df0_dv_XYZ_host[0]<<std::endl;
#endif
        // E1(x) 
        transform( particleWorkList.begin(), particleWorkList.end(), E1.begin(), 
                        getPerturbedField(&r[0],&e1_CYL[0],r.size(),nPhi,ky,kz,hanningWeight[i],wrf,thisT[i]) ); 
        
#ifdef __CUDACC__
        std::cout<<"E1 CPU: "<<E1[0]<<" GPU: "<<E1_host[0]<<std::endl;
#endif
        // B1(x) 
        transform( particleWorkList.begin(), particleWorkList.end(), B1.begin(), 
                        getPerturbedField(&r[0],&b1_CYL[0],r.size(),nPhi,ky,kz,hanningWeight[i],wrf,thisT[i]) ); 

#ifdef __CUDACC__
        std::cout<<"B1 CPU: "<<B1[0]<<" GPU: "<<B1_host[0]<<std::endl;
#endif
        // v x B1 
        transform( particleWorkList.begin(), particleWorkList.end(), B1.begin(), vCrossB.begin(), 
                        vCross<std::complex<float> >() );

#ifdef __CUDACC__
        std::cout<<"vCross CPU: "<<vCrossB[0]<<" GPU: "<<vCrossB_host[0]<<std::endl;
#endif
        // E1 + v x B1
        transform( E1.begin(), E1.end(), vCrossB.begin(), vCrossB_E1.begin(), 
                        std::plus<C3<std::complex<float> > >() );

#ifdef __CUDACC__
        std::cout<<"vCrosB_E1 CPU: "<<vCrossB_E1[0]<<" GPU: "<<vCrossB_E1_host[0]<<std::endl;
#endif
        //  (E1 + v x B1) . grad_v(f0(v))
        transform( vCrossB_E1.begin(), vCrossB_E1.end(), df0_dv_XYZ.begin(), forceDotGradf0.begin(), 
                        doDotProduct() );

#ifdef __CUDACC__
        std::cout<<"forceDotGradf0 CPU: "<<forceDotGradf0[0]<<" GPU: "<<forceDotGradf0_host[0]<<std::endl;
#endif
        // int( (E1 + v x B1) . grad_v(f0(v)), dt ) via running dt integral
        transform( dtIntegral.begin(), dtIntegral.end(), forceDotGradf0.begin(), dtIntegral.begin(), 
                        runningIntegral<std::complex<float> >(dtIntFac) );
#ifdef __CUDACC__
        std::cout<<"dtIntegral CPU: "<<dtIntegral[0]<<" GPU: "<<dtIntegral_host[0]<<std::endl;
#endif
        // f1(v) = -q/m * int( (E1 + v x B1) . grad_v(f0(v)), dt )
        transform( dtIntegral.begin(), dtIntegral.end(), particleWorkList.begin(), f1.begin(), 
                        multiplyByChargeOverMass<std::complex<float> >() ); 

        // q . f1(v) // first step in velocity momemnt for current calculation 
        transform( f1.begin(), f1.end(), particleWorkList.begin(), f1.begin(), 
                        multiplyByCharge<std::complex<float> >() ); 

        // q . vx . f1(v) 
        transform( f1.begin(), f1.end(), vx.begin(), vxf1.begin(), 
                        std::multiplies< complex<float> >() ); 

#ifdef __CUDACC__
        std::cout<<"CPU: "<<vxf1[0]<<" GPU: "<<vxf1_host[0]<<std::endl;
#endif
        // q . vy . f1(v) 
        transform( f1.begin(), f1.end(), vy.begin(), vyf1.begin(), 
                        std::multiplies< complex<float> >() ); 

        // q . vz . f1(v) 
        transform( f1.begin(), f1.end(), vz.begin(), vzf1.begin(), 
                        std::multiplies< complex<float> >() ); 
#endif

    }

    // Reduce velocity space to current via the first velocity moment

#if DO_CPU_ITERATOR_APPROACH > 0
    for (int i=0;i<nXGrid;i++) {
        j1xc[i] = dv * accumulate( vxf1.begin()+nP*i, vxf1.begin()+nP*i+nP, complex<float>(0) );
        j1yc[i] = dv * accumulate( vyf1.begin()+nP*i, vyf1.begin()+nP*i+nP, complex<float>(0) );
        j1zc[i] = dv * accumulate( vzf1.begin()+nP*i, vzf1.begin()+nP*i+nP, complex<float>(0) );
        std::cout << j1xc[i].real() << "  " << j1xc[i].imag() << std::endl;
    }
#endif

#ifdef __CUDACC__

    // Copy data back from GPU

    thrust::copy(vxf1_device.begin(),vxf1_device.end(),vxf1_host.begin());
    thrust::copy(vyf1_device.begin(),vyf1_device.end(),vyf1_host.begin());
    thrust::copy(vzf1_device.begin(),vzf1_device.end(),vzf1_host.begin());

    for (int i=0;i<nXGrid;i++) {
        j1xc[i] = dv * accumulate( vxf1_host.begin()+nP*i, vxf1_host.begin()+nP*i+nP, thrust::complex<float>(0) );
        j1yc[i] = dv * accumulate( vyf1_host.begin()+nP*i, vyf1_host.begin()+nP*i+nP, thrust::complex<float>(0) );
        j1zc[i] = dv * accumulate( vzf1_host.begin()+nP*i, vzf1_host.begin()+nP*i+nP, thrust::complex<float>(0) );
        std::cout << j1xc[i].real() << "  " << j1xc[i].imag() << std::endl;
    }
#endif

    stringstream ncjPFileName2("output/jP2.nc");

    NcFile ncjPFile(ncjPFileName2.str().c_str(), NcFile::replace);

    NcDim nc_nX = ncjPFile.addDim("nJp", nXGrid);

    NcVar nc_x = ncjPFile.addVar("x", ncFloat, nc_nX);

    NcVar nc_j1xc_re = ncjPFile.addVar("j1xc_re", ncFloat, nc_nX);
    NcVar nc_j1xc_im = ncjPFile.addVar("j1xc_im", ncFloat, nc_nX);

    NcVar nc_j1yc_re = ncjPFile.addVar("j1yc_re", ncFloat, nc_nX);
    NcVar nc_j1yc_im = ncjPFile.addVar("j1yc_im", ncFloat, nc_nX);

    NcVar nc_j1zc_re = ncjPFile.addVar("j1zc_re", ncFloat, nc_nX);
    NcVar nc_j1zc_im = ncjPFile.addVar("j1zc_im", ncFloat, nc_nX);

    vector<float> JxRe(nXGrid,0);
    vector<float> JxIm(nXGrid,0);
    vector<float> JyRe(nXGrid,0);
    vector<float> JyIm(nXGrid,0);
    vector<float> JzRe(nXGrid,0);
    vector<float> JzIm(nXGrid,0);

    for (int i=0;i<nXGrid;i++) {
       JxRe[i] = j1xc[i].real(); 
       JxIm[i] = j1xc[i].imag(); 
       JyRe[i] = j1yc[i].real(); 
       JyIm[i] = j1yc[i].imag(); 
       JzRe[i] = j1zc[i].real(); 
       JzIm[i] = j1zc[i].imag(); 
    }
    nc_x.putVar(&xGrid[0]);
    nc_j1xc_re.putVar(&JxRe[0]);
    nc_j1xc_im.putVar(&JxIm[0]);
    nc_j1yc_re.putVar(&JyRe[0]);
    nc_j1yc_im.putVar(&JyIm[0]);
    nc_j1zc_re.putVar(&JzRe[0]);
    nc_j1zc_im.putVar(&JzIm[0]);

    std::cout << "DONE" << std::endl;

#if DO_CPU_APPROACH > 0

#if CLOCK >= 1
#if not defined(_OPENMP)
        clock_t endTimeFunctor = clock();
        double timeInSecondsFunctor = (endTimeFunctor - startTimeFunctor) / (double)CLOCKS_PER_SEC;
        std::cout << "Time for this spatial point: " << timeInSecondsFunctor << std::endl;
        std::cout << "Time per particle: " << timeInSecondsFunctor / nWork << std::endl;
#endif
#endif

std::cout << "Continuing with non functor approach ..." << std::endl;

int write_iX = 0;//31;//15;
int write_iP = 180;//52;//33;

#pragma omp parallel for private(istat, tid, spoken)
    for (int iX = 0; iX < nXGrid; iX++) {

#if defined(_OPENMP)
        nThreads = omp_get_num_threads();
        tid = omp_get_thread_num();
        if (tid == 0 && spoken == 0) {
            std::cout << "tid : " << tid << std::endl;
            std::cout << "OMP_NUM_THREADS: " << nThreads << std::endl;
            spoken = 1;
        }
#endif
        vector<CParticle> ThisParticleList(
            create_particles(xGrid[iX], amu, Z, T_keV[iX], density_m3[iX], nPx, nPy,
                nPz, nThermal, dv, &r[0], &b0_CYL[0], r.size() ));

#if CLOCK >= 1
        clock_t startTime = clock();
#endif
        j1xc[iX] = complex<float>(0, 0);
        j1yc[iX] = complex<float>(0, 0);
        j1zc[iX] = complex<float>(0, 0);

#if LOWMEM_USEPAPI >= 1
        cpuTime0 = cpuTime;
        realTime0 = realTime;
        flpIns0 = flpIns;
        papiReturn = PAPI_flops(&realTime, &cpuTime, &flpIns, &mFlops);
#endif

#if F1_WRITE >= 1
        int f1_write_iX = write_iX;
        ofstream f1File;
        if (iX == f1_write_iX) {
            f1File.open("output/f1.txt", ios::out | ios::trunc);
            f1File << " vx  vy  vz  re(f1) im(f1) " << std::endl;
        }
#endif

        vector<float> f1(nP);
        vector<complex<float> > f1c(nP);

        for (int iP = 0; iP < nP; iP++) {

            vector<C3<float> > thisOrbitE1_re_XYZ(nSteps, C3<float>(0, 0, 0));
            vector<C3<float> > thisOrbitE1_im_XYZ(nSteps, C3<float>(0, 0, 0));

            vector<C3<float> > thisOrbitB1_re_XYZ(nSteps, C3<float>(0, 0, 0));
            vector<C3<float> > thisOrbitB1_im_XYZ(nSteps, C3<float>(0, 0, 0));

            CParticle thisParticle_XYZ(ThisParticleList[iP]);

            float qOverm = thisParticle_XYZ.q / thisParticle_XYZ.m;

            float Ze = thisParticle_XYZ.q;
#if LOWMEM_ORBIT_WRITE >= 1
            ofstream OrbitFile;
            ofstream v1File;
            ofstream e1_dot_grad_File;
            ofstream df0dv_File;

            if (iX == write_iX && iP == write_iP) {
                std::cout << "Write Particle Properties:" << std::endl;
                std::cout << " vTh: " << thisParticle_XYZ.vTh << std::endl;
                std::cout << " v1: " << thisParticle_XYZ.v_c1 << std::endl;
                std::cout << " v2: " << thisParticle_XYZ.v_c2 << std::endl;
                std::cout << " v3: " << thisParticle_XYZ.v_c3 << std::endl;

                OrbitFile.open("output/orbit.txt", ios::out | ios::trunc);
                OrbitFile << "wc / wrf: " << wrf_wc[iX] << std::endl;
                OrbitFile << " t  x  y  z  re(e1)  im(e1)  re(e2)  im(e2)  re(e3)  "
                             "im(e3)  re(b1)  im(b1)  re(b2)  im(b2)  re(b3)  im(b3) "
                             "status"
                          << std::endl;
                v1File.open("output/orbit_v1.txt", ios::out | ios::trunc);
                v1File << " t  re(v11)  im(v11)  re(v12)  im(v12)  re(v13)  im(v13)"
                       << std::endl;
                e1_dot_grad_File.open("output/orbit_e1_dot_grad_df0_dv.txt",
                    ios::out | ios::trunc);
                e1_dot_grad_File << " t  re(e1.Gradvf0)  im(e1.Gradvf0)  re(e1_per.Gradvf0_per)  "
                                    "im(e1_per.Gradvf0_per)  re(e1_par.Gradvf0_par)  im(e1_par.Gradvf0_par)"
                                 << std::endl;
                df0dv_File.open("output/df0dv.txt", ios::out | ios::trunc);
                df0dv_File << " t  vx  vy  vz  valp  vbet  vpar  vper  gyroAngle  "
                              "df0dv_x  df0dv_y  df0dv_z"
                           << std::endl;
            }
#endif
            // generate orbit and get time-harmonic e along it

            vector<C3<float> > thisOrbit_XYZ(nSteps);
            vector<C3<std::complex<float> > > thisE1c_XYZ(nSteps, C3<std::complex<float> >());
            vector<C3<std::complex<float> > > thisB1c_XYZ(nSteps, C3<std::complex<float> >());
            C3<std::complex<float> > thisV1c_(0, 0, 0), dVc(0, 0, 0), crossTerm(0, 0, 0);
            vector<complex<float> > this_e1_dot_gradvf0(nSteps);
            vector<complex<float> > this_e1_dot_gradvf0_parOnly(nSteps);
            vector<complex<float> > this_e1_dot_gradvf0_perOnly(nSteps);
            vector<C3<std::complex<float> > > this_vCrossB1(nSteps);

            for (int i = 0; i < nSteps; i++) {
#if DEBUG_MOVE >= 1
                std::cout << "Position Before Move: " << thisParticle_XYZ.c1 << "  "
                     << thisParticle_XYZ.c2 << "  " << thisParticle_XYZ.c3 << std::endl;
                std::cout << "p.status: " << thisParticle_XYZ.status << std::endl;
#endif
                thisOrbit_XYZ[i] = C3<float>(thisParticle_XYZ.c1, thisParticle_XYZ.c2,
                    thisParticle_XYZ.c3);
#if GC_ORBITS >= 1
                int MoveStatus = rk4_move_gc(thisParticle_XYZ, dtMin, thisT[i], &r[0], &b0_CYL[0], r.size(), 
                                &r_gc[0], &curv_CYL[0], &grad_CYL[0], &bDotGradB[0], r_gc.size());
#else
                int MoveStatus = rk4_move(thisParticle_XYZ, dtMin, &r[0], &b0_CYL[0], r.size());
#endif
                int OverallStatus = max(thisParticle_XYZ.status, MoveStatus);
#if DEBUG_MOVE >= 1
                std::cout << "Position After Move: " << thisParticle_XYZ.c1 << "  "
                         << thisParticle_XYZ.c2 << "  " << thisParticle_XYZ.c3 << std::endl;
                if (MoveStatus > 0) {
                    std::cout << "MoveStatus: " << MoveStatus << std::endl;
                }
#endif
                C3<float> thisPos_XYZ(thisParticle_XYZ.c1, thisParticle_XYZ.c2, thisParticle_XYZ.c3);
                C3<float> thisPos_CYL = XYZ_to_CYL(thisPos_XYZ);
                C3<float> thisB0 = kj_interp1D(thisPos_CYL.c1, &r[0], &b0_CYL[0], r.size(), thisParticle_XYZ.status);

#if GC_ORBITS >= 1
                C3<float> par = thisB0 / mag(thisB0);
                C3<float> per = cross(par,C3<float>(1,0,0));
                per = per / mag(per);

                // We can just pick any gyrophase and use that to compute gradv_f0, 
                // since gradvf0.par is independent of gyrophase (proven in the commented
                // section below)

                C3<float> this_v = par * thisParticle_XYZ.vPar + per * thisParticle_XYZ.vPer;
                C3<float> this_gradv_f0_GC = ( maxwellian_df0_dv(this_v, T_keV[iX], density_m3[iX], thisParticle_XYZ.amu, thisParticle_XYZ.Z) );
                float this_gradv_f0_GC_dot_par = dot(this_gradv_f0_GC,par);
 
                //// Prove that gradv_f0_par is independent of gyrophase
                //std::vector< C3<float> > gradv_f0_XYZ_GC;
                //int nTh = 12;
                //float dTh = 360.0 / nTh;
                //for(int iTh=0; iTh<nTh; iTh++){
                //        float th = iTh*dTh;
                //        C3<float> thisPer = rot_axis_angle(per,par,th); 
                //        C3<float> this_v_2 = par * thisParticle_XYZ.vPar + thisPer * thisParticle_XYZ.vPer;

                //        gradv_f0_XYZ_GC.push_back( maxwellian_df0_dv(this_v_2, T_keV[iX], density_m3[iX], thisParticle_XYZ.amu, thisParticle_XYZ.Z) );
                //        if(thisParticle_XYZ.status>0) gradv_f0_XYZ_GC[iTh] = 0; // To account for the / 0 above.

                //        std::cout<<"iTh: "<<iTh<<std::endl;
                //        std::cout<<"th: "<<th<<std::endl;
                //        std::cout<<"thisPer: "<<thisPer<<std::endl;
                //        std::cout<<"mag(thisPer): "<<mag(thisPer)<<std::endl;
                //        std::cout<<"mag(thisVel): "<<mag(this_v_2)<<std::endl;
                //        std::cout<<"mag(v): "<<std::sqrt(std::pow(thisParticle_XYZ.vPar,2)+std::pow(thisParticle_XYZ.vPer,2))<<std::endl;
                //        std::cout<<"this_v_2: "<<this_v_2<<std::endl;
                //        std::cout<<"grad_f0: "<<gradv_f0_XYZ_GC[iTh]<<std::endl;
                //        std::cout<<"grad_f0: "<<maxwellian_df0_dv(this_v_2, T_keV[iX], density_m3[iX], thisParticle_XYZ.amu, thisParticle_XYZ.Z)<<std::endl;
                //        std::cout<<"grad_f0.par: "<<dot(gradv_f0_XYZ_GC[iTh],par)<<std::endl;

                //        C3<float> this_v_2_above = par * thisParticle_XYZ.vPar + per * thisParticle_XYZ.vPer;
                //        C3<float> this_gradv_f0_GC_above = ( maxwellian_df0_dv(this_v_2_above, T_keV[iX], density_m3[iX], thisParticle_XYZ.amu, thisParticle_XYZ.Z) );
                //        float this_gradv_f0_GC_dot_par = dot(this_gradv_f0_GC_above,par);
                //        std::cout<<"grad_f0.par (actual): "<<this_gradv_f0_GC_dot_par<<std::endl;

                //}

                //exit(1);

#if DEBUG_MOVE >= 2
                std::cout << "vPar: " << thisParticle_XYZ.vPar << " vPer: " << thisParticle_XYZ.vPer << std::endl;
                std::cout << "status: " << thisParticle_XYZ.status << std::endl;
                std::cout << "c1: " << thisOrbit_XYZ[i].c1 << std::endl;
                std::cout << "thisB0: " << thisB0 << std::endl;
                std::cout << "mag(thisB0): " << mag(thisB0) << std::endl;
                std::cout << "thisParticle_XYZ.vPar: " << thisParticle_XYZ.vPar << std::endl;
                std::cout << "thisVel_XYZ: " << thisVel_XYZ << std::endl;
#endif

#else
                C3<float> thisVel_XYZ(thisParticle_XYZ.v_c1, thisParticle_XYZ.v_c2, thisParticle_XYZ.v_c3);
                C3<float> gradv_f0_XYZ = maxwellian_df0_dv(thisVel_XYZ, T_keV[iX], density_m3[iX], thisParticle_XYZ.amu, thisParticle_XYZ.Z);
                if(thisParticle_XYZ.status>0) gradv_f0_XYZ = 0; // To account for the / 0 above.
#endif

                complex<float> _i(0, 1);

                // why is this exp(-iwt) here? surely it's not required for the freq domain calc?

                C3<std::complex<float> > E1_XYZ;
#if CYLINDRICAL_INPUT_FIELDS >=1 
                E1_XYZ = hanningWeight[i] * exp(-_i * wrf * thisT[i]) * getE1orB1_XYZ_fromCYL(thisParticle_XYZ, &r[0], &e1_CYL[0], r.size(), nPhi);
#else
                E1_XYZ = hanningWeight[i] * exp(-_i * wrf * thisT[i]) * getE1orB1_XYZ_fromXYZ(thisParticle_XYZ, &r[0], &e1_CYL[0], r.size(), ky, kz);
#endif
                thisE1c_XYZ[i] = E1_XYZ * (1 - thisParticle_XYZ.status);

                C3<std::complex<float> > B1_XYZ;
#if CYLINDRICAL_INPUT_FIELDS >=1 
                B1_XYZ = hanningWeight[i] * exp(-_i * wrf * thisT[i]) * getE1orB1_XYZ_fromCYL(thisParticle_XYZ, &r[0], &b1_CYL[0], r.size(), nPhi);
#else
                B1_XYZ = hanningWeight[i] * exp(-_i * wrf * thisT[i]) * getE1orB1_XYZ_fromXYZ(thisParticle_XYZ, &r[0], &b1_CYL[0], r.size(), ky, kz);
#endif
                thisB1c_XYZ[i] = B1_XYZ * (1 - thisParticle_XYZ.status);

                //if (iX == write_iX && iP == write_iP) {
                //    std::cout<<"E1_XYZ: "<<E1_XYZ<<std::endl;
                //    std::cout<<"status: "<<thisParticle_XYZ.status<<std::endl;
                //    std::cout<<"hanningWeight: "<<hanningWeight[i]<<std::endl;
                //}

#if DEBUG_MOVE >= 2
                std::cout << "thisE1c[i].c1: " << thisE1c_XYZ[i].c1 << std::endl;
                std::cout << "thisE1c[i].c2: " << thisE1c_XYZ[i].c2 << std::endl;
                std::cout << "thisE1c[i].c3: " << thisE1c_XYZ[i].c3 << std::endl;

                std::cout << "thisB1c[i].c1: " << thisB1c_XYZ[i].c1 << std::endl;
                std::cout << "thisB1c[i].c2: " << thisB1c_XYZ[i].c2 << std::endl;
                std::cout << "thisB1c[i].c3: " << thisB1c_XYZ[i].c3 << std::endl;
#endif
#if DEBUG_FORCE_TERM >= 1
                std::cout << "thisE1c[i].c1: " << thisE1c_XYZ[i].c1 << std::endl;
                std::cout << "thisE1c[i].c2: " << thisE1c_XYZ[i].c2 << std::endl;
                std::cout << "thisE1c[i].c3: " << thisE1c_XYZ[i].c3 << std::endl;

                std::cout << "thisB1c[i].c1: " << thisB1c_XYZ[i].c1 << std::endl;
                std::cout << "thisB1c[i].c2: " << thisB1c_XYZ[i].c2 << std::endl;
                std::cout << "thisB1c[i].c3: " << thisB1c_XYZ[i].c3 << std::endl;

                std::cout << "thisVel_XYZ.c1: " << thisVel_XYZ.c1 << std::endl;
                std::cout << "thisVel_XYZ.c2: " << thisVel_XYZ.c2 << std::endl;
                std::cout << "thisVel_XYZ.c3: " << thisVel_XYZ.c3 << std::endl;

#endif

#if GC_ORBITS >= 1

                C3<std::complex<float> > this_force = thisE1c_XYZ[i];
                this_e1_dot_gradvf0[i] = dot(this_force,par) * this_gradv_f0_GC_dot_par;
                if(thisParticle_XYZ.status>0) this_e1_dot_gradvf0[i] = 0; // To account for the / 0 above.

#else
                this_vCrossB1[i] = cross(thisVel_XYZ, thisB1c_XYZ[i]);
                ////C3<std::complex<float> > this_force = this_vCrossB1[i] + thisE1c_XYZ[i];
                C3<std::complex<float> > this_force = thisE1c_XYZ[i];
                this_e1_dot_gradvf0[i] = dot(this_force, gradv_f0_XYZ);
                if(thisParticle_XYZ.status>0) this_e1_dot_gradvf0[i] = 0; // To account for the / 0 above.

                // Get the per and par contributions to the total e1DotGradvF0 for debugging 

                C3<float> par = thisB0 / mag(thisB0);
                C3<float> per = cross(par,C3<float>(1,0,0));
                per = per / mag(per);

                C3<std::complex<float> > par_force = par*dot(this_force,par);
                C3<float> par_gradf = par*dot(gradv_f0_XYZ,par);

                C3<std::complex<float> > per_force = this_force - par_force;
                C3<float> per_gradf = gradv_f0_XYZ - par_gradf;

                this_e1_dot_gradvf0_perOnly[i] = dot(per_force, per_gradf);
                this_e1_dot_gradvf0_parOnly[i] = dot(par_force, par_gradf);

#endif

#if LOWMEM_ORBIT_WRITE >= 1
#if GC_ORBITS == 0
                if (iX == write_iX && iP == write_iP) {
                    df0dv_File << scientific;
                    df0dv_File << thisT[i] << "    " << thisVel_XYZ.c1 << "    "
                               << thisVel_XYZ.c2 << "    " << thisVel_XYZ.c3 << "    "
                               << thisParticle_XYZ.vAlp << "    " << thisParticle_XYZ.vBet
                               << "    " << thisParticle_XYZ.vPar << "    "
                               << thisParticle_XYZ.vPer << "    " << thisParticle_XYZ.phs
                               << "    " << gradv_f0_XYZ.c1 << "    " << gradv_f0_XYZ.c2
                               << "    " << gradv_f0_XYZ.c3 << std::endl;
                }
#endif
                if (iX == write_iX && iP == write_iP) {
                    OrbitFile << scientific;
                    OrbitFile << thisT[i] << "    " 
                              << thisPos_XYZ.c1 << "    " 
                              << thisPos_XYZ.c2 << "    " 
                              << thisPos_XYZ.c3 << "    " 
                              << real(thisE1c_XYZ[i].c1) << "    " 
                              << imag(thisE1c_XYZ[i].c1) << "    "
                              << real(thisE1c_XYZ[i].c2) << "    "
                              << imag(thisE1c_XYZ[i].c2) << "    "
                              << real(thisE1c_XYZ[i].c3) << "    "
                              << imag(thisE1c_XYZ[i].c3) << "    "
                              << real(thisB1c_XYZ[i].c1) << "    "
                              << imag(thisB1c_XYZ[i].c1) << "    "
                              << real(thisB1c_XYZ[i].c2) << "    "
                              << imag(thisB1c_XYZ[i].c2) << "    "
                              << real(thisB1c_XYZ[i].c3) << "    "
                              << imag(thisB1c_XYZ[i].c3) << "    "
                              << real(this_vCrossB1[i].c1) << "    "
                              << imag(this_vCrossB1[i].c1) << "    "
                              << real(this_vCrossB1[i].c2) << "    "
                              << imag(this_vCrossB1[i].c2) << "    "
                              << real(this_vCrossB1[i].c3) << "    "
                              << imag(this_vCrossB1[i].c3) << "    "
                              << thisParticle_XYZ.status << std::endl;
                }
                if (iX == write_iX && iP == write_iP) {
                    e1_dot_grad_File << scientific;
                    e1_dot_grad_File << thisT[i] 
                                    << "    " << real(this_e1_dot_gradvf0[i])
                                    << "    " << imag(this_e1_dot_gradvf0[i]) 
                                    << "    " << real(this_e1_dot_gradvf0_perOnly[i])
                                    << "    " << imag(this_e1_dot_gradvf0_perOnly[i])
                                    << "    " << real(this_e1_dot_gradvf0_parOnly[i])
                                    << "    " << imag(this_e1_dot_gradvf0_parOnly[i])
                                    << std::endl;
                }
#endif
            }
#if LOWMEM_ORBIT_WRITE >= 1
            if (iX == write_iX && iP == write_iP) {
                OrbitFile.close();
            }
#endif

            complex<float> this_f1c = -qOverm * intVecArray(thisT, this_e1_dot_gradvf0);

#if GC_ORBITS >= 1

            // Add the offset to the GC time integration 

            C3<float> StartingPos_XYZ(ThisParticleList[iP].c1, ThisParticleList[iP].c2, ThisParticleList[iP].c3);
            C3<float> StartingPos_CYL = XYZ_to_CYL(StartingPos_XYZ)
 
            C3<float> thisB0 = kj_interp1D(StartingPos_CYL.c1, &r[0], &b0_CYL[0], r.size(), thisParticle_XYZ.status);
            C3<float> par = thisB0 / mag(thisB0);
            C3<float> per = cross(par,C3<float>(1,0,0));
            per = per / mag(per);

            // Calculate the initial gyrophase offset to be added to
            // the guiding center calculation ...

            C3<std::complex<float> > initial_force = thisE1c_XYZ[0];
            C3<float> initialV_XYZ(ThisParticleList[iP].v_c1, ThisParticleList[iP].v_c2, ThisParticleList[iP].v_c3);
            C3<float> initial_gradv_f0_XYZ = maxwellian_df0_dv(initialV_XYZ, T_keV[iX], density_m3[iX], thisParticle_XYZ.amu, thisParticle_XYZ.Z);

            // Angle between perp only components

            C3<std::complex<float> > perp_force = initial_force - par*dot(initial_force,par);
            C3<float> perp_gradf = initial_gradv_f0_XYZ - par*dot(initial_gradv_f0_XYZ,par);

            C3<float> perp_force_re(perp_force.c1.real(),perp_force.c2.real(),perp_force.c3.real());
            C3<float> perp_force_im(perp_force.c1.imag(),perp_force.c2.imag(),perp_force.c3.imag());

            float this_angle_perp_re = std::acos(dot(perp_force_re,perp_gradf) / (mag(perp_force_re)*mag(perp_gradf)));
            float this_angle_perp_im = std::acos(dot(perp_force_im,perp_gradf) / (mag(perp_force_im)*mag(perp_gradf)));

            // what is the total angle, not just the perp one?

            std::complex<float> this_angle_perp = std::complex<float>(this_angle_perp_re,this_angle_perp_im);

            // Normalize angle to the 0<th<360 range

            if(this_angle_perp_re<0) this_angle_perp_re += 2 * physConstants::pi;
            if(this_angle_perp_im<0) this_angle_perp_im += 2 * physConstants::pi;

            // Double check the angle
            if( (this_angle_perp_re<0) || (this_angle_perp_re>2*physConstants::pi) ) {
                    std::cout<<"ERROR: angle out of range"<<std::endl;
                    exit(1);
            }

            if( (this_angle_perp_im<0) || (this_angle_perp_im>2*physConstants::pi) ) {
                    std::cout<<"ERROR: angle out of range"<<std::endl;
                    exit(1);
            }

            float angleEnd_re = 0;
            float angleEnd_im = 0;

            if(this_angle_perp_re>=0) angleEnd_re = physConstants::pi;
            if(this_angle_perp_re>=physConstants::pi) angleEnd_re = 2*physConstants::pi;

            if(this_angle_perp_im>=0) angleEnd_im = physConstants::pi/2;
            if(this_angle_perp_im>=physConstants::pi/2) angleEnd_im = physConstants::pi/2+physConstants::pi;
            if(this_angle_perp_im>=physConstants::pi/2+physConstants::pi) angleEnd_im = physConstants::pi/2+2*physConstants::pi;

            float offsetReal = +(( std::sin(angleEnd_re) - std::sin(this_angle_perp_re) ) * mag(perp_force)*mag(perp_gradf)).real();
            float offsetImag = -(( std::cos(this_angle_perp_im) - std::cos(angleEnd_im) ) * mag(perp_force)*mag(perp_gradf)).imag();

            // Convert gyro angle integral to time integral

            float offset_wc = std::abs(ThisParticleList[iP].q*mag(thisB0)/ThisParticleList[iP].m);
            float offset_period = 2.0f*physConstants::pi/offset_wc;
            std::complex<float> offset_dt = this_angle_perp / float(2.0f*physConstants::pi)*offset_period;
 
            std::complex<float> offset = std::complex<float>(offsetReal,offsetImag) * dtMin * float(2*physConstants::pi);

            // Account for this_force==0 due to hanningWeight==0 at last point, or e1.gradvf0==0.
            if(isnan(this_angle_perp.real())) {
                offset = 0;
            }

            this_f1c += -qOverm * offset;

            if (iX == write_iX && iP == write_iP) {
                    std::cout<<"Offset: "<<offset<<std::endl;
            } 

            //average_e1_dot_gradvf0 = dot(this_force,par) * dot(initial_gradv_f0_XYZ,par);

            if (iX == write_iX && iP == write_iP) {
                std::cout<<"Offset wc: "<<offset_wc<<std::endl;
                std::cout<<"Offset T: "<<offset_period<<std::endl;
                std::cout<<"Offset_dt: "<<offset_dt<<std::endl;
                std::cout<<"dtMin: "<<dtMin<<std::endl;
                std::cout<<"par: "<<par<<std::endl;
                std::cout<<"thisE1c_XYZ: "<<thisE1c_XYZ[0]<<std::endl;
                std::cout<<"this_force: "<<initial_force<<std::endl;
                std::cout<<"this_angle_perp: "<<this_angle_perp*float(180.0f/physConstants::pi)<<std::endl;
                std::cout<<"perp_force: "<<perp_force<<std::endl;
                std::cout<<"perp_gradv: "<<perp_gradf<<std::endl;
                std::cout<<"initial_gradv_f0_XYZ:"<<initial_gradv_f0_XYZ<<std::endl;
            }
#endif

            if (iX == write_iX && iP == write_iP) {
                    //for(int i=0; i<nSteps;i++){
                    //    std::cout<<"this_e1_dot_gradvf0[i]: "<<this_e1_dot_gradvf0[i]<<std::endl;
                    //}
                    std::cout<<"this_f1c: "<<this_f1c<<std::endl;
            }

#if LOWMEM_ORBIT_WRITE >= 1
            if (iX == write_iX && iP == write_iP) {

                complex<float> tmp = 0.0;
                for (int i = 0; i < nSteps; i++) {
                    tmp += -qOverm * this_e1_dot_gradvf0[i] * dtMin;
                    v1File << thisT[i] << "    " << real(tmp) << "    " << imag(tmp)
                           << std::endl;
                }
            }
#endif
            f1c[iP] = -this_f1c;

            float v0x_i = ThisParticleList[iP].v_c1;
            float v0y_i = ThisParticleList[iP].v_c2;
            float v0z_i = ThisParticleList[iP].v_c3;

            float h = dv * Ze;

#pragma omp critical // "atomic" does not work for complex numbers
            {
                j1xc[iX] += h * (v0x_i * f1c[iP]);
                j1yc[iX] += h * (v0y_i * f1c[iP]);
                j1zc[iX] += h * (v0z_i * f1c[iP]);
#if DEBUG_MOVE >= 2
                std::cout << "v0x_i: " << v0x_i << std::endl;
                std::cout << "v0y_i: " << v0y_i << std::endl;
                std::cout << "v0z_i: " << v0z_i << std::endl;

                std::cout << "f1c[iP]: " << f1c[iP] << std::endl;
                std::cout << "qOverm: " << qOverm << std::endl;
                std::cout << "dtMin: " << dtMin << std::endl;

                std::cout << "j1xc[iX]: " << j1xc[iX] << std::endl;
                std::cout << "j1yc[iX]: " << j1yc[iX] << std::endl;
                std::cout << "j1zc[iX]: " << j1zc[iX] << std::endl;
                //exit(1);
#endif
            }

            if (iX == write_iX && iP == write_iP) {
                std::cout<<"iX: "<<iX<<std::endl;
                std::cout<<"iP: "<<iP<<std::endl;
                std::cout<<"h: "<<h<<std::endl;
                std::cout<<"dv: "<<dv<<std::endl;
                std::cout<<"f1c[iP]: "<<f1c[iP]<<std::endl;
                std::cout<<"j1xc[iX]: "<<j1xc[iX]<<std::endl;
                std::cout<<"j1yc[iX]: "<<j1yc[iX]<<std::endl;
                std::cout<<"j1zc[iX]: "<<j1zc[iX]<<std::endl;
            }

#if F1_WRITE >= 1
            if (iX == f1_write_iX) {
                f1File << scientific;
                f1File << showpos;
                f1File << v0x_i << "    " << v0y_i << "    " << v0z_i << "    "
                       << real(f1c[iP]) << "    " << imag(f1c[iP]) << std::endl;
            }
#endif
        }

#if CLOCK >= 1
#if not defined(_OPENMP)
        clock_t endTime = clock();
        double timeInSeconds = (endTime - startTime) / (double)CLOCKS_PER_SEC;
        std::cout << "Time for this spatial point: " << timeInSeconds << std::endl;
        std::cout << "Time per particle: " << timeInSeconds / nP << std::endl;
#endif
#endif

#if LOWMEM_USEPAPI >= 1
        cpuTime0 = cpuTime;
        realTime0 = realTime;
        flpIns0 = flpIns;
        papiReturn = PAPI_flops(&realTime, &cpuTime, &flpIns, &mFlops);
        printf("\nLOWMEM Oribit calculation performance ...\n");
        printf(
            "Real_time:\t%f\nProc_time:\t%f\nTotal flpins:\t%lld\nMFLOPS:\t\t%f\n",
            realTime - realTime0, cpuTime - cpuTime0, flpIns - flpIns0, mFlops);
#endif

#if USEPAPI >= 1
        printf("\nGet e(t) and integrate performance ...\n");
        printf(
            "Real_time:\t%f\nProc_time:\t%f\nTotal flpins:\t%lld\nMFLOPS:\t\t%f\n",
            eT_realTime, eT_cpuTime, eT_flpIns, eT_mFlops );
        printf("\nGet v(t) and integrate performance ...\n");
        printf(
            "Real_time:\t%f\nProc_time:\t%f\nTotal flpins:\t%lld\nMFLOPS:\t\t%f\n",
            vT_realTime, vT_cpuTime, vT_flpIns, vT_mFlops );

        cpuTime0 = cpuTime;
        realTime0 = realTime;
        flpIns0 = flpIns;
        papiReturn = PAPI_flops(&realTime, &cpuTime, &flpIns, &mFlops);
        printf("\nj(t) performance ...\n");
        printf(
            "Real_time:\t%f\nProc_time:\t%f\nTotal flpins:\t%lld\nMFLOPS:\t\t%f\n",
            realTime - realTime0, cpuTime - cpuTime0, flpIns - flpIns0, mFlops);
#endif


    } // End of xGrid loop

    // Write current(s) to file

    // std::cout << "Writing jP to file ... ";

    for (int iX = 0; iX < nXGrid; iX++) {

        stringstream ncjPFileName;
        ncjPFileName << "output/";
        // check directory exists
        struct stat st;
        if (stat(ncjPFileName.str().c_str(), &st) != 1) {
            int mkDirStat = mkdir(ncjPFileName.str().c_str(),
                S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
        }
        ncjPFileName << "/jP_";
        ncjPFileName << setw(3) << setfill('0') << iX;
        ncjPFileName << ".nc";
#if DEBUGLEVEL >= 1
        std::cout << ncjPFileName.str().c_str() << std::endl;
#endif

        NcFile ncjPFile(ncjPFileName.str().c_str(), NcFile::replace);

        NcDim nc_scalar = ncjPFile.addDim("scalar", 1);

        NcVar nc_x = ncjPFile.addVar("x", ncFloat, nc_scalar);
        NcVar nc_freq = ncjPFile.addVar("freq", ncFloat, nc_scalar);

        NcVar nc_j1xc_re = ncjPFile.addVar("j1xc_re", ncFloat, nc_scalar);
        NcVar nc_j1xc_im = ncjPFile.addVar("j1xc_im", ncFloat, nc_scalar);

        NcVar nc_j1yc_re = ncjPFile.addVar("j1yc_re", ncFloat, nc_scalar);
        NcVar nc_j1yc_im = ncjPFile.addVar("j1yc_im", ncFloat, nc_scalar);

        NcVar nc_j1zc_re = ncjPFile.addVar("j1zc_re", ncFloat, nc_scalar);
        NcVar nc_j1zc_im = ncjPFile.addVar("j1zc_im", ncFloat, nc_scalar);

        nc_x.putVar(&xGrid[iX]);
        nc_freq.putVar(&freq);

        vector<size_t> startp(1, 0);

        float tmpJxRe = real(j1xc[iX]);
        float tmpJxIm = imag(j1xc[iX]);
        nc_j1xc_re.putVar(&tmpJxRe);
        nc_j1xc_im.putVar(&tmpJxIm);

        float tmpJyRe = real(j1yc[iX]);
        float tmpJyIm = imag(j1yc[iX]);
        nc_j1yc_re.putVar(&tmpJyRe);
        nc_j1yc_im.putVar(&tmpJyIm);

        float tmpJzRe = real(j1zc[iX]);
        float tmpJzIm = imag(j1zc[iX]);
        nc_j1zc_re.putVar(&tmpJzRe);
        nc_j1zc_im.putVar(&tmpJzIm);

        if (iX == write_iX) {
            std::cout<<"write_iX"<<std::endl;
        }

        std::cout<<"j1xc[iX]: "<<j1xc[iX]<<std::endl;
        std::cout<<"j1yc[iX]: "<<j1yc[iX]<<std::endl;
        std::cout<<"j1zc[iX]: "<<j1zc[iX]<<std::endl;

    }

    // ProfilerStop();

    std::cout << "DONE" << std::endl;

#if CLOCK >= 1
    clock_t ProgramTime_ = clock();
    double ProgramTimeInSeconds = (ProgramTime_ - ProgramTime) / (double)CLOCKS_PER_SEC;
#if defined(_OPENMP)
    ProgramTimeInSeconds = ProgramTimeInSeconds / nThreads;
    std::cout << "nThreads: " << nThreads << std::endl;
#endif
    std::cout << "Total Time [s]: " << ProgramTimeInSeconds << std::endl;
    std::cout << "Total Time [m]: " << ProgramTimeInSeconds / 60.0 << std::endl;
    std::cout << "Total Time [h]: " << ProgramTimeInSeconds / 3600.0 << std::endl;
#endif

#endif
    return EXIT_SUCCESS;
}
