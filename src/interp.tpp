template<class TYPE>
TYPE kj_interp1D ( const float &x, const vector<float> &xVec, const vector<TYPE> &yVec, int &status ) {

	float _x, x0, x1;
	float xTmp;
	xTmp = x;

#if _PARTICLE_BOUNDARY == 1
	if(x<xVec.front()||x>xVec.back()) {
			// Particle absorbing walls
#if DEBUG_INTERP >= 1
            if(xVec.size()!=yVec.size()) {
#if DEBUG_INTERP >= 2
                cout<<"ERROR: xVec and yVec are not the same size for interpolation!"<<endl;
#endif
                status = 1;
                return TYPE(0);
            }
#if DEBUG_INTERP >= 2
			cout<<"Particle absorbed at "<<x<<endl;
            cout<<"Particle number: "<<p.number<<endl;
            cout<<"x:"<<x<<endl;
            cout<<"xVec.front():"<<xVec.front()<<endl;
            cout<<"xVec.back():"<<xVec.back()<<endl;
#endif
#endif
			status = 1;
			return TYPE(0);
	}
#elif _PARTICLE_BOUNDARY == 2
			// Periodic 
			if(xTmp<xVec.front()) xTmp = xVec.back()-(xVec.front()-xTmp);			
			if(xTmp>xVec.back()) xTmp = xVec.front()+(xTmp-xVec.back());			
#elif _PARTICLE_BOUNDARY == 3
			// Particle reflecting walls
			if(xTmp<xVec.front()) xTmp = xVec.front()+(xVec.front()-xTmp);			
			if(xTmp>xVec.back()) xTmp = xVec.back()-(xTmp-xVec.back());			
#endif

#if DEBUG_INTERP >= 1    
	if(status>0){
#if DEBUG_INTERP >= 2
			cout<<"ERROR: Should never get here with _PARTICLE_BOUNDARY ==2|3"<<endl;
#endif
			return TYPE(0);
	}
#endif
	_x = (xTmp-xVec.front())/(xVec.back()-xVec.front())*(xVec.size()-1);

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
		TYPE y0 = yVec[x0];
		TYPE y1 = yVec[x1];
#if DEBUG_INTERP >=2
        //cout << "kj_interp1D: " << endl;
        //if(typeid(TYPE)==typeid(C3Vec)) cout << "Type is C3Vec" << endl;
        //if(typeid(TYPE)==typeid(float)) cout << "Type is float" << endl;
        //kj_print(x0,"x0");
        //kj_print(x1,"x1");
        //kj_print(y0,"y0");
        //kj_print(y1,"y1");
        //cout << endl;
        if(x0>yVec.size()-1||x0<0||x1>yVec.size()-1||x1<1) {
                cout<<"ERROR: interpolation point off the end of array"<<endl;
                cout<<"x.front: "<<xVec.front()<<endl;
                cout<<"x.back: "<<xVec.back()<<endl;
                cout<<"x: "<<x<<endl;
                cout<<"xVec.size(): "<<xVec.size()<<endl;
                cout<<"yVec.size(): "<<yVec.size()<<endl;
                status=1;
                return TYPE(0);
        }
#endif
        TYPE result = y0+(_x-x0)*(y1-y0)/(x1-x0);

#if DEBUG_INTERP >=1
        if(isnan(result)) {
#if DEBUG_INTERP >= 2
                cout<<"ERROR: interpolation produced a NaN"<<endl;
#endif
                status=1;
                return TYPE(0);
        } 
        if(isinf(result)) {
#if DEBUG_INTERP >= 2
                cout<<"ERROR: interpolation produced a INF"<<endl;
#endif
                status = 1;
                return TYPE(0);
        }
#endif
		return result;
	}

}
