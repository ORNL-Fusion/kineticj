// Converted this to accept pointers to the x, y(x) arrays since thrust 
// functors can only use pointers

PRAGMA
template <typename T>
HOST DEVICE
T kj_interp1D ( const float &x, const float *xVec, const T *yVec, int n, int &status ) {

    int x0, x1;
	float _x;

#if _PARTICLE_BOUNDARY == 1
	if(x<xVec[0]||x>xVec[n-1]) {
			// Particle absorbing walls
#if DEBUG_INTERP >= 1
            if(n!=n) {
#if DEBUG_INTERP >= 2
                std::cout<<"ERROR: xVec and yVec are not the same size for interpolation!"<<std::endl;
#endif
                status = 1;
                return T(0);
            }
#if DEBUG_INTERP >= 2
			std::cout<<"Particle absorbed at "<<x<<std::endl;
            std::cout<<"x:"<<x<<std::endl;
            std::cout<<"xVec[0]:"<<xVec[0]<<std::endl;
            std::cout<<"xVec[n-1]:"<<xVec[n-1]<<std::endl;
#endif
#endif
			status = 1;
			return T(0);
	}
#elif _PARTICLE_BOUNDARY == 2
			// Periodic 
			if(x<xVec[0]) x = xVec[n-1]-(xVec[0]-x);			
			if(x>xVec[n-1]) x = xVec[0]+(x-xVec[n-1]);			
#elif _PARTICLE_BOUNDARY == 3
			// Particle reflecting walls
			if(x<xVec[0]) x = xVec[0]+(xVec[0]-x);			
			if(x>xVec[n-1]) x = xVec[n-1]-(x-xVec[n-1]);			
#endif

#if DEBUG_INTERP >= 1    
	if(status>0){
#if DEBUG_INTERP >= 2
			std::cout<<"ERROR: Should never get here with _PARTICLE_BOUNDARY ==2|3"<<std::endl;
#endif
			return T(0);
	}
#endif
	_x = (x-xVec[0])/(xVec[n-1]-xVec[0])*(n-1);

	x0 = floor(_x);
	x1 = ceil(_x);

	// Catch for particle at point
	if(x0==x1) {
#if DEBUG_INTERP >= 2
        std::cout << "Particle version of kj_interp1D" << std::endl;
		std::cout << "x0: " << x0 << " x1: " <<x1<< " _x: "<<_x << std::endl;
		std::cout << "Particle at point catch: " << x0/x1 << "  "  << abs(1.0-x0/x1) << std::endl;
#endif
		return yVec[x0];
	}
	else {
		T y0 = yVec[x0];
		T y1 = yVec[x1];
#if DEBUG_INTERP >=1
        if(x0>n-1||x0<0||x1>n-1||x1<1) {
#ifdef __CUDACC__
                printf("ERROR: interpolation point off the end of array\n");
                printf("x.front: %f\n",xVec[0]);
                printf("x.back: %f\n",xVec[n-1]);
                printf("x: %f\n",x);
                printf("n: %f\n",n);
#else
                std::cout<<"ERROR: interpolation point off the end of array"<<std::endl;
                std::cout<<"x.front: "<<xVec[0]<<std::endl;
                std::cout<<"x.back: "<<xVec[n-1]<<std::endl;
                std::cout<<"x: "<<x<<std::endl;
                std::cout<<"n: "<<n<<std::endl;
#endif
                status=1;
                return T(0);
        }
#endif

        T result = y0+(_x-x0)*(y1-y0)/(x1-x0);

#ifndef __CUDACC__
#if DEBUG_INTERP >=1
        if(isnan(result)) {
#if DEBUG_INTERP >= 2
                std::cout<<"ERROR: interpolation produced a NaN"<<std::endl;
#endif
                status=1;
                return T(0);
        } 
        if(isinf(result)) {
#if DEBUG_INTERP >= 2
                std::cout<<"ERROR: interpolation produced a INF"<<std::endl;
#endif
                status = 1;
                return T(0);
        }
#endif
#endif
		return result;
	}

}
