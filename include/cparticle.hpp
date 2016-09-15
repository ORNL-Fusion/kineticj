#ifndef CPARTICLE_HPP
#define CPARTICLE_HPP

#include "cspecies.hpp"

class CParticle: public CSpecies {
		public:
				float c1, c2, c3, v_c1, v_c2, v_c3;
				int number;
				float weight;
				int status;
                float dvx, dvy, dvz, d3v;
                float vPar, vPer, gyroPhase, u, vTh;
                float vAlp, vBet, phs;
                float T, n; // Temp and density of the Maxwellian at this point

				CParticle ();
				CParticle ( double _amu, int _Z);
				CParticle (float c1, float c2, float c3, 
								float v_c1, float v_c2, float v_c3, 
								double _amu, int _Z, float _weight );
				CParticle (float c1, float c2, float c3, 
								float v_c1, float v_c2, float v_c3, 
								double _amu, int _Z, float _weight, float _T, float _n );
				CParticle (CSpecies _species);
};

#endif
