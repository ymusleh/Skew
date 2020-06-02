#include "Skewpoly.h"
#include "SkewRing.h"

class SkewRing;
class Skewpoly;

class Drinfeld {

public:
	Drinfeld(Skewpoly& phi);
	Drinfeld(SkewRing* ring, Skewpoly phi);
	SkewRing* getRing();
	Skewpoly getPhi_T();
	ZZ_pX char_poly_det_rank2();
	ZZ_pEX charPoly2();
	

private:
	SkewRing* ring;
	// The Drinfeld map
	Skewpoly phi_T;
	int rank, n, m;
	// The inclusion map \gamma(T)
	ZZ_pE inclusion;
	// A-Characteristic
	ZZ_pX frakp;
};
