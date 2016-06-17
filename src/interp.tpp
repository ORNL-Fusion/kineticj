// Converted this to accept pointers to the x, y(x) arrays since thrust 
// functors can only use pointers

template<class T>
T kj_interp1D ( const float &x, const float *xVec, const T *yVec, int n, int &status ) {

    int x0, x1;
	float _x;
	float xTmp;
	xTmp = x;

#if _PARTICLE_BOUNDARY == 1
	if(x<xVec[0]||x>xVec[n-1]) {
			// Particle absorbing walls
#if DEBUG_INTERP >= 1
            if(n!=n) {
#if DEBUG_INTERP >= 2
                cout<<"ERROR: xVec and yVec are not the same size for interpolation!"<<endl;
#endif
                status = 1;
                return T(0);
            }
#if DEBUG_INTERP >= 2
			cout<<"Particle absorbed at "<<x<<endl;
            cout<<"Particle number: "<<p.number<<endl;
            cout<<"x:"<<x<<endl;
            cout<<"xVec[0]:"<<xVec[0]<<endl;
            cout<<"xVec[n-1]:"<<xVec[n-1]<<endl;
#endif
#endif
			status = 1;
			return T(0);
	}
#elif _PARTICLE_BOUNDARY == 2
			// Periodic 
			if(xTmp<xVec[0]) xTmp = xVec[n-1]-(xVec[0]-xTmp);			
			if(xTmp>xVec[n-1]) xTmp = xVec[0]+(xTmp-xVec[n-1]);			
#elif _PARTICLE_BOUNDARY == 3
			// Particle reflecting walls
			if(xTmp<xVec[0]) xTmp = xVec[0]+(xVec[0]-xTmp);			
			if(xTmp>xVec[n-1]) xTmp = xVec[n-1]-(xTmp-xVec[n-1]);			
#endif

#if DEBUG_INTERP >= 1    
	if(status>0){
#if DEBUG_INTERP >= 2
			cout<<"ERROR: Should never get here with _PARTICLE_BOUNDARY ==2|3"<<endl;
#endif
			return T(0);
	}
#endif
	_x = (xTmp-xVec[0])/(xVec[n-1]-xVec[0])*(n-1);

	x0 = floor(_x);
	x1 = ceil(_x);
	
	// Catch for particle at point
	if(x0==x1) {
#if DEBUG_INTERP >= 2
        cout << "Particle version of kj_interp1D" << endl;
		cout << "x0: " << x0 << " x1: " <<x1<< " _x: "<<_x << endl;
		cout << "Particle at point catch: " << x0/x1 << "  "  << abs(1.0-x0/x1) << endl;
#endif
		return yVec[x0];
	}
	else {
		T y0 = yVec[x0];
		T y1 = yVec[x1];
#if DEBUG_INTERP >=2
        //cout << "kj_interp1D: " << endl;
        //if(Tid(T)==Tid(C3Vec)) cout << "T is C3Vec" << endl;
        //if(Tid(T)==Tid(float)) cout << "T is float" << endl;
        //kj_print(x0,"x0");
        //kj_print(x1,"x1");
        //kj_print(y0,"y0");
        //kj_print(y1,"y1");
        //cout << endl;
        if(x0>n-1||x0<0||x1>n-1||x1<1) {
                cout<<"ERROR: interpolation point off the end of array"<<endl;
                cout<<"x.front: "<<xVec[0]<<endl;
                cout<<"x.back: "<<xVec[n-1]<<endl;
                cout<<"x: "<<x<<endl;
                cout<<"n: "<<n<<endl;
                cout<<"n: "<<n<<endl;
                status=1;
                return T(0);
        }
#endif
        T result = y0+(_x-x0)*(y1-y0)/(x1-x0);

#if DEBUG_INTERP >=1
        if(isnan(result)) {
#if DEBUG_INTERP >= 2
                cout<<"ERROR: interpolation produced a NaN"<<endl;
#endif
                status=1;
                return T(0);
        } 
        if(isinf(result)) {
#if DEBUG_INTERP >= 2
                cout<<"ERROR: interpolation produced a INF"<<endl;
#endif
                status = 1;
                return T(0);
        }
#endif
		return result;
	}

}
