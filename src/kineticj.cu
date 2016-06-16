#include "c3vec.hpp"
#ifdef __CUDACC__
#include <thrust/device_vector.h>
#include <thrust/host_vector.h>
#endif


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

using namespace std;
using namespace netCDF;
using namespace exceptions;

// Calculate the jP given some know E and f(v)

int main(int argc, char** argv)
{

    // Make sure the "output/" directory exists

    stringstream outputDirName;
    outputDirName << "output/";

    // check directory exists
    struct stat st;
    int dirTest = stat(outputDirName.str().c_str(), &st);
    if (dirTest != 0) {
        cout << "Had to create output/ directory" << endl;
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
        cout << "ERROR: PAPI Failed to initialize with error code: " << papiReturn
             << endl;
        cout << "ERROR: See papi.h for error code explanations " << endl;
        exit(1);
    }
    printf("Real_time:\t%f\nProc_time:\t%f\nTotal flpins:\t%lld\nMFLOPS:\t\t%f\n",
        realTime, cpuTime, flpIns, mFlops);

    papiReturn = PAPI_flops(&realTime, &cpuTime, &flpIns, &mFlops);
    if (papiReturn < 0) {
        cout << "ERROR: PAPI Failed to initialize with error code: " << papiReturn
             << endl;
        cout << "ERROR: See papi.h for error code explanations " << endl;
        exit(1);
    } else {
        cout << "PAPI called successfully with return code: " << papiReturn << endl;
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

    // Read E
    string eField_fName = cfg.lookup("eField_fName");
    vector<C3Vec> e1Re_CYL, e1Im_CYL, b1Re_CYL, b1Im_CYL;
    vector<C3VecI> e1_CYL, b1_CYL;
    vector<C3Vec> b0_CYL, b0_XYZ;
    vector<float> r, n_m3;
    float freq;
    int eReadStat = read_e_field(eField_fName, species_number, freq, r, n_m3, e1_CYL, b1_CYL,
        e1Re_CYL, e1Im_CYL, b1Re_CYL, b1Im_CYL, b0_CYL);

    // Read GC terms
    string gc_fName = cfg.lookup("gc_fName");
    vector<C3Vec> curv_CYL, grad_CYL;
    std::vector<float> r_gc, bDotGradB;
    int gcReadStat = read_gc_file(gc_fName, r_gc, curv_CYL, grad_CYL, bDotGradB);

    float wrf = freq * 2 * physConstants::pi;
    float xGridMin = cfg.lookup("xGridMin");
    float xGridMax = cfg.lookup("xGridMax");
    int nXGrid = cfg.lookup("nXGrid");
    cout << "nXGrid: " << nXGrid << endl;

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
        int iStat;
        density_m3[iX] = kj_interp1D(xGrid[iX], r, n_m3, iStat);
        C3Vec this_b0 = kj_interp1D(xGrid[iX], r, b0_CYL, iStat);
        bMag_kjGrid[iX] = mag(this_b0);
        T_keV[iX] = 2.0; // kj_interp1D(xGrid[iX],r,n_m3);
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
    float nStepsPerCycle = cfg.lookup("nStepsPerCycle");
    float tRF = (2 * physConstants::pi) / wrf;
    int nJpCycles = cfg.lookup("nJpCycles");
    int nJpPerCycle = cfg.lookup("nJpPerCycle");
    int nPhi = cfg.lookup("nPhi");
    int nJp = nJpCycles * nJpPerCycle;
    float dtJp = tRF / nJpPerCycle;
    int istat = 0;
    int nPx = cfg.lookup("nPx");
    int nPy = cfg.lookup("nPy");
    int nPz = cfg.lookup("nPz");
    float amu = cfg.lookup("species_amu");
    float Z = cfg.lookup("species_Z");
    int nThermal = cfg.lookup("nThermal");
    long int nP = nPx * nPy * nPz;
    float wc = Z * physConstants::e * MaxB0 / (amu * physConstants::mi);
    float cyclotronPeriod = 2 * physConstants::pi / wc;
    float dtMin = -cyclotronPeriod / nStepsPerCycle;
    int nSteps = nRFCycles * tRF / abs(dtMin) + 1;

    for (int iX = 0; iX < nXGrid; iX++) {
        float this_wc = Z * physConstants::e * bMag_kjGrid[iX] / (amu * physConstants::mi);
        wrf_wc[iX] = wrf / this_wc;
    }

#if PRINT_INFO >= 1
    cout << "dtMin [s]: " << dtMin << endl;
    cout << "Cyclotron Period: " << cyclotronPeriod << endl;
    cout << "RF Period: " << tRF << endl;
    cout << "nSteps: " << nSteps << endl;
    cout << "nStepsPerCycle: " << nStepsPerCycle << endl;
    cout << "freq: " << freq << endl;
    cout << "Max B0: " << MaxB0 << endl;
#endif

    vector<float> thisT;
    try {
        thisT.resize(nSteps);
    } catch (const std::bad_alloc& error) {
        cout << "Allocation error at " << __FILE__ << __LINE__ << endl;
        cout << error.what();
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
        // expWeight[i] = 1.0;//abs(exp(-_i*wrf_c*thisT[i]));
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
                nPz, nThermal, dv, r, b0_CYL));

        particleWorkList.insert( particleWorkList.end(), moreWork.begin(), moreWork.end() );
    }

#ifdef __CUDACC__
    //thrust::device_vector<CParticle> particleWorkList_device = particleWorkList;
#endif

    // Create the vx,vy,vz iterators

    vector<float> vx(nWork,0);
    vector<float> vy(nWork,0);
    vector<float> vz(nWork,0);

    transform( vx.begin(), vx.end(), particleWorkList.begin(), vx.begin(), set_vx() );
    transform( vy.begin(), vy.end(), particleWorkList.begin(), vy.begin(), set_vy() );
    transform( vz.begin(), vz.end(), particleWorkList.begin(), vz.begin(), set_vz() );

    // Move particles
    cout << "Moving particles with for_each ..." << endl;

    // Velocity space calculation

    vector<C3Vec> df0_dv_XYZ(nWork,0);
    vector<C3VecI> E1(nWork,0);
    vector<C3VecI> B1(nWork,0);
    vector<C3VecI> vCrossB(nWork,0);
    vector<C3VecI> vCrossB_E1(nWork,0);
    vector<complex<float> > forceDotGradf0(nWork,0);
    vector<complex<float> > dtIntegral(nWork,0);
    vector<complex<float> > f1(nWork,0);
    vector<complex<float> > vxf1(nWork,0);
    vector<complex<float> > vyf1(nWork,0);
    vector<complex<float> > vzf1(nWork,0);

    for (int i = 0; i < nSteps; i++) {

        float dtIntFac = 1;
        if (i > 0) dtIntFac = 2;

        dtIntFac = dtMin / 2.0 * dtIntFac;

        // Move particle
        for_each( particleWorkList.begin(), particleWorkList.end(), moveParticle(dtMin, r, b0_CYL) ); 
       
        // df0(v)/dv 
        transform( particleWorkList.begin(), particleWorkList.end(), df0_dv_XYZ.begin(), get_df0_dv() ); 

        // E1(x) 
        transform( particleWorkList.begin(), particleWorkList.end(), E1.begin(), getPerturbedField(r,e1_CYL,nPhi,hanningWeight[i]) ); 

        // B1(x) 
        transform( particleWorkList.begin(), particleWorkList.end(), B1.begin(), getPerturbedField(r,b1_CYL,nPhi,hanningWeight[i]) ); 

        // v x B1 
        transform( particleWorkList.begin(), particleWorkList.end(), B1.begin(), vCrossB.begin(), vCross() );

        // E1 + v x B1
        transform( E1.begin(), E1.end(), vCrossB.begin(), vCrossB_E1.begin(), std::plus<C3VecI>() );

        //  (E1 + v x B1) . grad_v(f0(v))
        transform( vCrossB_E1.begin(), vCrossB_E1.end(), df0_dv_XYZ.begin(), forceDotGradf0.begin(), doDotProduct() );

        // int( (E1 + v x B1) . grad_v(f0(v)), dt ) via running dt integral
        transform( dtIntegral.begin(), dtIntegral.end(), forceDotGradf0.begin(), dtIntegral.begin(), runningIntegral(dtIntFac) );

        // f1(v) = -q/m * int( (E1 + v x B1) . grad_v(f0(v)), dt )
        transform( dtIntegral.begin(), dtIntegral.end(), particleWorkList.begin(), f1.begin(), multiplyByChargeOverMass() ); 

        // q . f1(v) // first step in velocity momemnt for current calculation 
        transform( f1.begin(), f1.end(), particleWorkList.begin(), f1.begin(), multiplyByCharge() ); 

        // q . vx . f1(v) 
        transform( f1.begin(), f1.end(), vx.begin(), vxf1.begin(), std::multiplies< complex<float> >() ); 

        // q . vy . f1(v) 
        transform( f1.begin(), f1.end(), vy.begin(), vyf1.begin(), std::multiplies< complex<float> >() ); 

        // q . vz . f1(v) 
        transform( f1.begin(), f1.end(), vz.begin(), vzf1.begin(), std::multiplies< complex<float> >() ); 
    }

    // Reduce velocity space to current via the first velocity moment

    for (int i=0;i<nXGrid;i++) {
        j1xc[i] = dv * accumulate( vxf1.begin()+nP*i, vxf1.begin()+nP*i+nP, complex<float>(0) );
        j1yc[i] = dv * accumulate( vyf1.begin()+nP*i, vyf1.begin()+nP*i+nP, complex<float>(0) );
        j1zc[i] = dv * accumulate( vzf1.begin()+nP*i, vzf1.begin()+nP*i+nP, complex<float>(0) );
        cout << j1xc[i].real() << "  " << j1xc[i].imag() << endl;
    }

    stringstream ncjPFileName2("jP2.nc");

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

    cout << "DONE" << endl;

#if CLOCK >= 1
#if not defined(_OPENMP)
        clock_t endTimeFunctor = clock();
        double timeInSecondsFunctor = (endTimeFunctor - startTimeFunctor) / (double)CLOCKS_PER_SEC;
        cout << "Time for this spatial point: " << timeInSecondsFunctor << endl;
        cout << "Time per particle: " << timeInSecondsFunctor / nWork << endl;
#endif
#endif

cout << "Continuing with non functor approach ..." << endl;

#pragma omp parallel for private(istat, tid, spoken)
    for (int iX = 0; iX < nXGrid; iX++) {

#if defined(_OPENMP)
        nThreads = omp_get_num_threads();
        tid = omp_get_thread_num();
        if (tid == 0 && spoken == 0) {
            cout << "tid : " << tid << endl;
            cout << "OMP_NUM_THREADS: " << nThreads << endl;
            spoken = 1;
        }
#endif
        vector<CParticle> ThisParticleList(
            create_particles(xGrid[iX], amu, Z, T_keV[iX], density_m3[iX], nPx, nPy,
                nPz, nThermal, dv, r, b0_CYL));

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
        int f1_write_iX = 75;
        ofstream f1File;
        if (iX == f1_write_iX) {
            f1File.open("output/f1.txt", ios::out | ios::trunc);
            f1File << " vx  vy  vz  re(f1) im(f1) " << endl;
        }
#endif

        vector<float> f1(nP);
        vector<complex<float> > f1c(nP);

        for (int iP = 0; iP < nP; iP++) {

            vector<C3Vec> thisOrbitE1_re_XYZ(nSteps, C3Vec(0, 0, 0));
            vector<C3Vec> thisOrbitE1_im_XYZ(nSteps, C3Vec(0, 0, 0));

            vector<C3Vec> thisOrbitB1_re_XYZ(nSteps, C3Vec(0, 0, 0));
            vector<C3Vec> thisOrbitB1_im_XYZ(nSteps, C3Vec(0, 0, 0));

            CParticle thisParticle_XYZ(ThisParticleList[iP]);

            float qOverm = thisParticle_XYZ.q / thisParticle_XYZ.m;

            float Ze = thisParticle_XYZ.q;
#if LOWMEM_ORBIT_WRITE >= 1
            ofstream OrbitFile;
            ofstream v1File;
            ofstream e1_dot_grad_File;
            ofstream df0dv_File;

            int write_iX = 75;
            int write_iP = 33;
            if (iX == write_iX && iP == write_iP) {
                cout << "Write Particle Properties:" << endl;
                cout << " vTh: " << thisParticle_XYZ.vTh << endl;
                cout << " v1: " << thisParticle_XYZ.v_c1 << endl;
                cout << " v2: " << thisParticle_XYZ.v_c2 << endl;
                cout << " v3: " << thisParticle_XYZ.v_c3 << endl;

                OrbitFile.open("output/orbit.txt", ios::out | ios::trunc);
                OrbitFile << "wc / wrf: " << wrf_wc[iX] << endl;
                OrbitFile << " t  x  y  z  re(e1)  im(e1)  re(e2)  im(e2)  re(e3)  "
                             "im(e3)  re(b1)  im(b1)  re(b2)  im(b2)  re(b3)  im(b3) "
                             "status"
                          << endl;
                v1File.open("output/orbit_v1.txt", ios::out | ios::trunc);
                v1File << " t  re(v11)  im(v11)  re(v12)  im(v12)  re(v13)  im(v13)"
                       << endl;
                e1_dot_grad_File.open("output/orbit_e1_dot_grad_df0_dv.txt",
                    ios::out | ios::trunc);
                e1_dot_grad_File << " t  re(v1xb01)  im(v1xb01)  re(v1xb02)  "
                                    "im(v1xb02)  re(v1xb03)  im(v1xb03)"
                                 << endl;
                df0dv_File.open("output/df0dv.txt", ios::out | ios::trunc);
                df0dv_File << " t  vx  vy  vz  valp  vbet  vpar  vper  gyroAngle  "
                              "df0dv_x  df0dv_y  df0dv_z"
                           << endl;
            }
#endif
            // generate orbit and get time-harmonic e along it

            vector<C3Vec> thisOrbit_XYZ(nSteps);
            vector<C3VecI> thisE1c_XYZ(nSteps, C3VecI());
            vector<C3VecI> thisB1c_XYZ(nSteps, C3VecI());
            C3VecI thisV1c_(0, 0, 0), dVc(0, 0, 0), crossTerm(0, 0, 0);
            vector<complex<float> > this_e1_dot_gradvf0(nSteps);
            vector<C3VecI> this_vCrossB1(nSteps);

            for (int i = 0; i < nSteps; i++) {
#if DEBUG_MOVE >= 1
                cout << "Position Before Move: " << thisParticle_XYZ.c1 << "  "
                     << thisParticle_XYZ.c2 << "  " << thisParticle_XYZ.c3 << endl;
                cout << "p.status: " << thisParticle_XYZ.status << endl;
#endif
                thisOrbit_XYZ[i] = C3Vec(thisParticle_XYZ.c1, thisParticle_XYZ.c2,
                    thisParticle_XYZ.c3);
#if GC_ORBITS >= 1
                int MoveStatus = rk4_move_gc(thisParticle_XYZ, dtMin, thisT[i], r, b0_CYL, r_gc,
                    curv_CYL, grad_CYL, bDotGradB, wrf);
#else
                int MoveStatus = rk4_move(thisParticle_XYZ, dtMin, r, b0_CYL);
#endif
                int OverallStatus = max(thisParticle_XYZ.status, MoveStatus);
#if DEBUG_MOVE >= 1
                if (MoveStatus > 0) {
                    cout << "Position After Move: " << thisParticle_XYZ.c1 << "  "
                         << thisParticle_XYZ.c2 << "  " << thisParticle_XYZ.c3 << endl;
                    cout << "ERROR: rk4_move* threw an error" << endl;
                    cout << "MoveStatus: " << MoveStatus << endl;
                    exit(1);
                }
#endif

                C3Vec thisPos(thisParticle_XYZ.c1, thisParticle_XYZ.c2,
                    thisParticle_XYZ.c3);
                C3Vec thisVel_XYZ(thisParticle_XYZ.v_c1, thisParticle_XYZ.v_c2,
                    thisParticle_XYZ.v_c3);
                C3Vec thisB0 = kj_interp1D(thisOrbit_XYZ[i].c1, r, b0_CYL, istat);
#if GC_ORBITS >= 1
                thisVel_XYZ = thisB0 / mag(thisB0) * thisParticle_XYZ.vPar; // vPar vector in XYZ
                cout << thisParticle_XYZ.vPar << "  " << thisParticle_XYZ.vPer << endl;
                kj_print(thisVel_XYZ, "thisVel_XYZ");
                C3Vec gradv_f0_XYZ = maxwellian_df0_dv(thisVel_XYZ, T_keV[iX], density_m3[iX],
                    thisParticle_XYZ.amu, thisParticle_XYZ.Z);
#else
                C3Vec gradv_f0_XYZ = maxwellian_df0_dv(thisVel_XYZ, T_keV[iX], density_m3[iX],
                    thisParticle_XYZ.amu, thisParticle_XYZ.Z);
#endif

                C3VecI E1_XYZ;
                complex<float> _i(0, 1);
                // why is this exp(-iwt) here? surely it's not required for the freq domain calc?
                //E1_XYZ = hanningWeight[i] * exp(-_i * wrf * thisT[i]) * getE1orB1_XYZ(thisParticle_XYZ, r, e1_CYL, nPhi);
                E1_XYZ = hanningWeight[i] * getE1orB1_XYZ(thisParticle_XYZ, r, e1_CYL, nPhi);
                thisE1c_XYZ[i] = E1_XYZ * (1 - thisParticle_XYZ.status);

                C3VecI B1_XYZ;
                //B1_XYZ = hanningWeight[i] * exp(-_i * wrf * thisT[i]) * getE1orB1_XYZ(thisParticle_XYZ, r, b1_CYL, nPhi);
                B1_XYZ = hanningWeight[i] * getE1orB1_XYZ(thisParticle_XYZ, r, b1_CYL, nPhi);
                thisB1c_XYZ[i] = B1_XYZ * (1 - thisParticle_XYZ.status);

#if DEBUG_MOVE >= 2
                cout << "thisE1c[i].c1: " << thisE1c_XYZ[i].c1 << endl;
                cout << "thisE1c[i].c2: " << thisE1c_XYZ[i].c2 << endl;
                cout << "thisE1c[i].c3: " << thisE1c_XYZ[i].c3 << endl;

                cout << "thisB1c[i].c1: " << thisB1c_XYZ[i].c1 << endl;
                cout << "thisB1c[i].c2: " << thisB1c_XYZ[i].c2 << endl;
                cout << "thisB1c[i].c3: " << thisB1c_XYZ[i].c3 << endl;
#endif
#if DEBUG_FORCE_TERM >= 1
                cout << "thisE1c[i].c1: " << thisE1c_XYZ[i].c1 << endl;
                cout << "thisE1c[i].c2: " << thisE1c_XYZ[i].c2 << endl;
                cout << "thisE1c[i].c3: " << thisE1c_XYZ[i].c3 << endl;

                cout << "thisB1c[i].c1: " << thisB1c_XYZ[i].c1 << endl;
                cout << "thisB1c[i].c2: " << thisB1c_XYZ[i].c2 << endl;
                cout << "thisB1c[i].c3: " << thisB1c_XYZ[i].c3 << endl;

                cout << "thisVel_XYZ.c1: " << thisVel_XYZ.c1 << endl;
                cout << "thisVel_XYZ.c2: " << thisVel_XYZ.c2 << endl;
                cout << "thisVel_XYZ.c3: " << thisVel_XYZ.c3 << endl;

#endif

#if GC_ORBITS >= 1
                // For GC orbits (electrons) use only the orbit parallel piece,
                // since the perp peice will cancel due to E not varying within
                // a cyclotron period.
                //
                // NO - WE DO NOT HAVE thisVel_XYZ for GC!!!!!
                C3Vec orbitParallelUnitVector_XYZ = thisB0 / mag(thisB0);
                // kj_print(orbitParallelUnitVector_XYZ,"unit");
                complex<float> orbitParallel_E = dot(thisE1c_XYZ[i], orbitParallelUnitVector_XYZ);
                float orbitParallel_gradv_f0 = dot(gradv_f0_XYZ, orbitParallelUnitVector_XYZ);
                // cout<<"E : "<<orbitParallel_E<<" gf0:
                // "<<orbitParallel_gradv_f0<<endl;
                this_e1_dot_gradvf0[i] = orbitParallel_E * orbitParallel_gradv_f0;
                // cout<<this_e1_dot_gradvf0[i]<<" e1dotgrad"<<endl;
                complex<float> _full = dot(thisE1c_XYZ[i], gradv_f0_XYZ);
                // cout<<"_full : "<<_full<<endl;
#else
                this_vCrossB1[i] = cross(thisVel_XYZ, thisB1c_XYZ[i]);
                C3VecI this_force = this_vCrossB1[i] + thisE1c_XYZ[i];

                // C3VecI this_force_CYL;
                // float this_t =
                // sqrt(pow(thisParticle_XYZ.c1,2)+pow(thisParticle_XYZ.c2,2));
                // this_force_CYL = rot_CYL_to_XYZ ( this_t, this_force, -1);
                // this_force_CYL.c1 = 0;
                // this_force_CYL.c3 = 0;
                // this_force = rot_CYL_to_XYZ ( this_t, this_force_CYL, +1);

                this_e1_dot_gradvf0[i] = dot(this_force, gradv_f0_XYZ);

                // C3Vec  this_gradv_f0_CYL;
                // this_force_CYL = rot_CYL_to_XYZ ( this_t, this_force, -1);
                // this_gradv_f0_CYL = rot_CYL_to_XYZ ( this_t, gradv_f0_XYZ, -1);
                // this_e1_dot_gradvf0[i] = dot(this_force_CYL, this_gradv_f0_CYL);
#endif

#if LOWMEM_ORBIT_WRITE >= 1
                if (iX == write_iX && iP == write_iP) {
                    df0dv_File << scientific;
                    df0dv_File << thisT[i] << "    " << thisVel_XYZ.c1 << "    "
                               << thisVel_XYZ.c2 << "    " << thisVel_XYZ.c3 << "    "
                               << thisParticle_XYZ.vAlp << "    " << thisParticle_XYZ.vBet
                               << "    " << thisParticle_XYZ.vPar << "    "
                               << thisParticle_XYZ.vPer << "    " << thisParticle_XYZ.phs
                               << "    " << gradv_f0_XYZ.c1 << "    " << gradv_f0_XYZ.c2
                               << "    " << gradv_f0_XYZ.c3 << endl;
                }

                if (iX == write_iX && iP == write_iP) {
                    OrbitFile << scientific;
                    OrbitFile << thisT[i] << "    " << thisPos.c1 << "    " << thisPos.c2
                              << "    " << thisPos.c3 << "    " << real(thisE1c_XYZ[i].c1)
                              << "    " << imag(thisE1c_XYZ[i].c1) << "    "
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
                              << thisParticle_XYZ.status << endl;
                }
                if (iX == write_iX && iP == write_iP) {
                    e1_dot_grad_File << scientific;
                    e1_dot_grad_File << thisT[i] << "    " << real(this_e1_dot_gradvf0[i])
                                     << "    " << imag(this_e1_dot_gradvf0[i]) << endl;
                }
#endif
            }
#if LOWMEM_ORBIT_WRITE >= 1
            if (iX == write_iX && iP == write_iP) {
                OrbitFile.close();
            }
#endif
            complex<float> this_f1c = -qOverm * intVecArray(thisT, this_e1_dot_gradvf0);

#if LOWMEM_ORBIT_WRITE >= 1
            if (iX == write_iX && iP == write_iP) {

                complex<float> tmp = 0.0;
                for (int i = 0; i < nSteps; i++) {
                    tmp += -qOverm * this_e1_dot_gradvf0[i] * dtMin;
                    v1File << thisT[i] << "    " << real(tmp) << "    " << imag(tmp)
                           << endl;
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
            }

#if F1_WRITE >= 1
            if (iX == f1_write_iX) {
                f1File << scientific;
                f1File << showpos;
                f1File << v0x_i << "    " << v0y_i << "    " << v0z_i << "    "
                       << real(f1c[iP]) << "    " << imag(f1c[iP]) << endl;
            }
#endif
        }

#if CLOCK >= 1
#if not defined(_OPENMP)
        clock_t endTime = clock();
        double timeInSeconds = (endTime - startTime) / (double)CLOCKS_PER_SEC;
        cout << "Time for this spatial point: " << timeInSeconds << endl;
        cout << "Time per particle: " << timeInSeconds / nP << endl;
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
            eT_realTime, eT_cpuTime, eT_flpIns, eT_mFlops / (nJp - 1));
        printf("\nGet v(t) and integrate performance ...\n");
        printf(
            "Real_time:\t%f\nProc_time:\t%f\nTotal flpins:\t%lld\nMFLOPS:\t\t%f\n",
            vT_realTime, vT_cpuTime, vT_flpIns, vT_mFlops / (nJp - 1));

        cpuTime0 = cpuTime;
        realTime0 = realTime;
        flpIns0 = flpIns;
        papiReturn = PAPI_flops(&realTime, &cpuTime, &flpIns, &mFlops);
        printf("\nj(t) performance ...\n");
        printf(
            "Real_time:\t%f\nProc_time:\t%f\nTotal flpins:\t%lld\nMFLOPS:\t\t%f\n",
            realTime - realTime0, cpuTime - cpuTime0, flpIns - flpIns0, mFlops);
#endif

#if __SAVE_ORBITS__ >= 1
        // Write orbits to file

        cout << "Writing orbits to file ... " << endl;

        stringstream ncOrbitsFileName;
        ncOrbitsFileName << "output/orbits_";
        ncOrbitsFileName << setw(3) << setfill('0') << iX;
        ncOrbitsFileName << ".nc";

        try {
            // Really need to fix this but I don't know how to
            // write a vector of structures using netCDF yet.

            NcFile ncOrbitsFile(ncOrbitsFileName.str().c_str(), NcFile::replace);

            NcDim nc_nP = ncOrbitsFile.addDim("nP", this_particles_XYZ.size());
            NcDim nc_nSteps = ncOrbitsFile.addDim("nSteps", nSteps);
            NcDim nc_nJp = ncOrbitsFile.addDim("nJp", nJp);

            vector<NcDim> nc_nPxnSteps(2);
            nc_nPxnSteps[0] = nc_nP;
            nc_nPxnSteps[1] = nc_nSteps;

            vector<NcDim> nc_nPxnJpxnSteps(3);
            nc_nPxnJpxnSteps[0] = nc_nP;
            nc_nPxnJpxnSteps[1] = nc_nJp;
            nc_nPxnJpxnSteps[2] = nc_nSteps;

            NcVar nc_t = ncOrbitsFile.addVar("t", ncFloat, nc_nSteps);

            NcVar nc_x = ncOrbitsFile.addVar("x", ncFloat, nc_nPxnSteps);
            NcVar nc_y = ncOrbitsFile.addVar("y", ncFloat, nc_nPxnSteps);
            NcVar nc_z = ncOrbitsFile.addVar("z", ncFloat, nc_nPxnSteps);

            NcVar nc_vx = ncOrbitsFile.addVar("vx", ncFloat, nc_nPxnSteps);
            NcVar nc_vy = ncOrbitsFile.addVar("vy", ncFloat, nc_nPxnSteps);
            NcVar nc_vz = ncOrbitsFile.addVar("vz", ncFloat, nc_nPxnSteps);

            NcVar nc_e1_x = ncOrbitsFile.addVar("e1_x", ncFloat, nc_nPxnSteps);
            NcVar nc_e1_y = ncOrbitsFile.addVar("e1_y", ncFloat, nc_nPxnSteps);
            NcVar nc_e1_z = ncOrbitsFile.addVar("e1_z", ncFloat, nc_nPxnSteps);

            NcVar nc_e1_x_re = ncOrbitsFile.addVar("e1_x_re", ncFloat, nc_nPxnSteps);
            NcVar nc_e1_y_re = ncOrbitsFile.addVar("e1_y_re", ncFloat, nc_nPxnSteps);
            NcVar nc_e1_z_re = ncOrbitsFile.addVar("e1_z_re", ncFloat, nc_nPxnSteps);

            NcVar nc_e1_x_im = ncOrbitsFile.addVar("e1_x_im", ncFloat, nc_nPxnSteps);
            NcVar nc_e1_y_im = ncOrbitsFile.addVar("e1_y_im", ncFloat, nc_nPxnSteps);
            NcVar nc_e1_z_im = ncOrbitsFile.addVar("e1_z_im", ncFloat, nc_nPxnSteps);

            NcVar nc_v1_x = ncOrbitsFile.addVar("v1x", ncFloat, nc_nPxnJpxnSteps);
            NcVar nc_v1_y = ncOrbitsFile.addVar("v1y", ncFloat, nc_nPxnJpxnSteps);
            NcVar nc_v1_z = ncOrbitsFile.addVar("v1z", ncFloat, nc_nPxnJpxnSteps);

            NcVar nc_v1_x_re = ncOrbitsFile.addVar("v1x_re", ncFloat, nc_nPxnJpxnSteps);
            NcVar nc_v1_y_re = ncOrbitsFile.addVar("v1y_re", ncFloat, nc_nPxnJpxnSteps);
            NcVar nc_v1_z_re = ncOrbitsFile.addVar("v1z_re", ncFloat, nc_nPxnJpxnSteps);

            NcVar nc_v1_x_im = ncOrbitsFile.addVar("v1x_im", ncFloat, nc_nPxnJpxnSteps);
            NcVar nc_v1_y_im = ncOrbitsFile.addVar("v1y_im", ncFloat, nc_nPxnJpxnSteps);
            NcVar nc_v1_z_im = ncOrbitsFile.addVar("v1z_im", ncFloat, nc_nPxnJpxnSteps);

            vector<size_t> startpA(2);
            vector<size_t> countpA(2);
            for (int iP = 0; iP < this_particles_XYZ.size(); iP++) {

                startpA[0] = iP;
                startpA[1] = 0;
                countpA[0] = 1;
                countpA[1] = nSteps;

                vector<float> tmpData(nSteps, 0);
                for (int iS = 0; iS < nSteps; iS++) {
                    tmpData[iS] = orbits_XYZ[iP][iS].c1;
                }
                nc_x.putVar(startpA, countpA, &tmpData[0]);
                for (int iS = 0; iS < nSteps; iS++) {
                    tmpData[iS] = orbits_XYZ[iP][iS].c2;
                }
                nc_y.putVar(startpA, countpA, &tmpData[0]);
                for (int iS = 0; iS < nSteps; iS++) {
                    tmpData[iS] = orbits_XYZ[iP][iS].c3;
                }
                nc_z.putVar(startpA, countpA, &tmpData[0]);

                for (int iS = 0; iS < nSteps; iS++) {
                    tmpData[iS] = orbits_v_XYZ[iP][iS].c1;
                }
                nc_vx.putVar(startpA, countpA, &tmpData[0]);
                for (int iS = 0; iS < nSteps; iS++) {
                    tmpData[iS] = orbits_v_XYZ[iP][iS].c2;
                }
                nc_vy.putVar(startpA, countpA, &tmpData[0]);
                for (int iS = 0; iS < nSteps; iS++) {
                    tmpData[iS] = orbits_v_XYZ[iP][iS].c3;
                }
                nc_vz.putVar(startpA, countpA, &tmpData[0]);

                for (int iS = 0; iS < nSteps; iS++) {
                    tmpData[iS] = e1[iP][iS].c1;
                }
                nc_e1_x.putVar(startpA, countpA, &tmpData[0]);
                for (int iS = 0; iS < nSteps; iS++) {
                    tmpData[iS] = e1[iP][iS].c2;
                }
                nc_e1_y.putVar(startpA, countpA, &tmpData[0]);
                for (int iS = 0; iS < nSteps; iS++) {
                    tmpData[iS] = e1[iP][iS].c3;
                }
                nc_e1_z.putVar(startpA, countpA, &tmpData[0]);

                for (int iS = 0; iS < nSteps; iS++) {
                    tmpData[iS] = real(e1c[iP][iS].c1);
                }
                nc_e1_x_re.putVar(startpA, countpA, &tmpData[0]);
                for (int iS = 0; iS < nSteps; iS++) {
                    tmpData[iS] = real(e1c[iP][iS].c2);
                }
                nc_e1_y_re.putVar(startpA, countpA, &tmpData[0]);
                for (int iS = 0; iS < nSteps; iS++) {
                    tmpData[iS] = real(e1c[iP][iS].c3);
                }
                nc_e1_z_re.putVar(startpA, countpA, &tmpData[0]);

                for (int iS = 0; iS < nSteps; iS++) {
                    tmpData[iS] = imag(e1c[iP][iS].c1);
                }
                nc_e1_x_im.putVar(startpA, countpA, &tmpData[0]);
                for (int iS = 0; iS < nSteps; iS++) {
                    tmpData[iS] = imag(e1c[iP][iS].c2);
                }
                nc_e1_y_im.putVar(startpA, countpA, &tmpData[0]);
                for (int iS = 0; iS < nSteps; iS++) {
                    tmpData[iS] = imag(e1c[iP][iS].c3);
                }
                nc_e1_z_im.putVar(startpA, countpA, &tmpData[0]);
            }

            vector<size_t> startpB(3);
            vector<size_t> countpB(3);
            for (int iP = 0; iP < this_particles_XYZ.size(); iP++) {
                for (int iJ = 0; iJ < nJp; iJ++) {

                    startpB[0] = iP;
                    startpB[1] = iJ;
                    startpB[2] = 0;
                    countpB[0] = 1;
                    countpB[1] = 1;
                    countpB[2] = nSteps;

                    vector<float> tmpData(nSteps, 0);

                    for (int iS = 0; iS < nSteps; iS++) {
                        tmpData[iS] = v1[iP][iJ][iS].c1;
                    }
                    nc_v1_x.putVar(startpB, countpB, &tmpData[0]);
                    for (int iS = 0; iS < nSteps; iS++) {
                        tmpData[iS] = v1[iP][iJ][iS].c2;
                    }
                    nc_v1_y.putVar(startpB, countpB, &tmpData[0]);
                    for (int iS = 0; iS < nSteps; iS++) {
                        tmpData[iS] = v1[iP][iJ][iS].c3;
                    }
                    nc_v1_z.putVar(startpB, countpB, &tmpData[0]);

                    for (int iS = 0; iS < nSteps; iS++) {
                        tmpData[iS] = real(v1c[iP][iJ][iS].c1);
                    }
                    nc_v1_x_re.putVar(startpB, countpB, &tmpData[0]);
                    for (int iS = 0; iS < nSteps; iS++) {
                        tmpData[iS] = real(v1c[iP][iJ][iS].c2);
                    }
                    nc_v1_y_re.putVar(startpB, countpB, &tmpData[0]);
                    for (int iS = 0; iS < nSteps; iS++) {
                        tmpData[iS] = real(v1c[iP][iJ][iS].c3);
                    }
                    nc_v1_z_re.putVar(startpB, countpB, &tmpData[0]);

                    for (int iS = 0; iS < nSteps; iS++) {
                        tmpData[iS] = imag(v1c[iP][iJ][iS].c1);
                    }
                    nc_v1_x_im.putVar(startpB, countpB, &tmpData[0]);
                    for (int iS = 0; iS < nSteps; iS++) {
                        tmpData[iS] = imag(v1c[iP][iJ][iS].c2);
                    }
                    nc_v1_y_im.putVar(startpB, countpB, &tmpData[0]);
                    for (int iS = 0; iS < nSteps; iS++) {
                        tmpData[iS] = imag(v1c[iP][iJ][iS].c3);
                    }
                    nc_v1_z_im.putVar(startpB, countpB, &tmpData[0]);
                }
            }

            vector<size_t> startp(1, 0);
            vector<size_t> countp(1, nSteps);

            nc_t.putVar(startp, countp, &thisT[0]);

        } catch (exceptions::NcException& e) {
            cout << "NetCDF: unknown error" << endl;
            e.what();
            exit(1);
        }

// cout << "DONE" << endl;
#endif

    } // End of xGrid loop

    // Write current(s) to file

    // cout << "Writing jP to file ... ";

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
        cout << ncjPFileName.str().c_str() << endl;
#endif

        NcFile ncjPFile(ncjPFileName.str().c_str(), NcFile::replace);

        NcDim nc_nJp = ncjPFile.addDim("nJp", nJp);
        NcDim nc_scalar = ncjPFile.addDim("scalar", 1);

        NcVar nc_t = ncjPFile.addVar("t", ncFloat, nc_nJp);

        NcVar nc_x = ncjPFile.addVar("x", ncFloat, nc_scalar);
        NcVar nc_freq = ncjPFile.addVar("freq", ncFloat, nc_scalar);

        NcVar nc_j1x = ncjPFile.addVar("j1x", ncFloat, nc_nJp);
        NcVar nc_j1y = ncjPFile.addVar("j1y", ncFloat, nc_nJp);
        NcVar nc_j1z = ncjPFile.addVar("j1z", ncFloat, nc_nJp);

        NcVar nc_j1xc_re = ncjPFile.addVar("j1xc_re", ncFloat, nc_scalar);
        NcVar nc_j1xc_im = ncjPFile.addVar("j1xc_im", ncFloat, nc_scalar);

        NcVar nc_j1yc_re = ncjPFile.addVar("j1yc_re", ncFloat, nc_scalar);
        NcVar nc_j1yc_im = ncjPFile.addVar("j1yc_im", ncFloat, nc_scalar);

        NcVar nc_j1zc_re = ncjPFile.addVar("j1zc_re", ncFloat, nc_scalar);
        NcVar nc_j1zc_im = ncjPFile.addVar("j1zc_im", ncFloat, nc_scalar);

        nc_x.putVar(&xGrid[iX]);
        nc_freq.putVar(&freq);

        vector<size_t> startp(1, 0);
        vector<size_t> countp(1, nJp);

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
    }

    // ProfilerStop();

    cout << "DONE" << endl;

#if CLOCK >= 1
    clock_t ProgramTime_ = clock();
    double ProgramTimeInSeconds = (ProgramTime_ - ProgramTime) / (double)CLOCKS_PER_SEC;
#if defined(_OPENMP)
    ProgramTimeInSeconds = ProgramTimeInSeconds / nThreads;
    cout << "nThreads: " << nThreads << endl;
#endif
    cout << "Total Time [s]: " << ProgramTimeInSeconds << endl;
    cout << "Total Time [m]: " << ProgramTimeInSeconds / 60.0 << endl;
    cout << "Total Time [h]: " << ProgramTimeInSeconds / 3600.0 << endl;
#endif
    return EXIT_SUCCESS;
}
