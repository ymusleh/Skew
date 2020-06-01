#include "Skewpoly.h"
#include "SkewRing.h"

class Drinfeld {

public:
	Drinfeld(SkewRing* ring);
	Drinfeld(SkewRing* ring, Skewpoly phit);
	SkewRing* getRing();
	Skewpoly getPhi_T();
	ZZ_pX charPoly();
	ZZ_pEX charPoly2();
	

private:
	SkewRing* ring;
	Skewpoly phi_T;
	ZZ_pE inclusion;
};
