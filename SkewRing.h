#pragma once
#include <NTL/ZZ_pE.h>
#include <NTL/ZZ_pEX.h>
#include <NTL/vec_ZZ_p.h>
#include <NTL/mat_ZZ_pE.h>

#include <NTL/ZZ_pXFactoring.h>
#include <NTL/ZZ_pEXFactoring.h>
//#include <math.h>
#include <iostream>
#include <algorithm>
#include <cmath>
#include <fstream>

using namespace NTL;
// using ZZ_pE = NTL::ZZ_pE;
// using ZZ_pEX = NTL::ZZ_pEX;

class SkewRing
{
public:
	SkewRing(const long, const ZZ_pX&, const ZZ_pEX&);
	// automatically constructs the frobenius 
	SkewRing(const long chr, const ZZ_pX &modulus);
	bool operator==(const SkewRing&);
	long baseFieldChar();
	ZZ_pX getMod();
	void setBaseFieldChar(const long);
	void setBaseField(ZZ_pX&);
	void setBaseField(int);
	void setSkewRelation(ZZ_pEX);
	ZZ_pEX getSkewRelation();
	void init();
	// Application of the Skew relation
	void applySkewRel(ZZ_pE& ret, const ZZ_pE& val);
	void fastSkewRel(ZZ_pE& ret, const ZZ_pE& val);
	void powerSkewRel(ZZ_pE& ret, const ZZ_pE& val, const long e);
	// Inverse Skew relation
	void fastInvRel(ZZ_pE& ret, const ZZ_pE& val);
	//void promote(ZZ_pEX& target, const ZZ_pX& source);
	//void initFrobenius(const ZZ &chr, const ZZ_pX &modulus);
	
//private:
	long chr;
	ZZ_pX modulus;
	ZZ_pX skewrel;
	ZZ_pEX skewRelation;
	// Modulus used for fast modular operations for the frobenius
	ZZ_pXModulus fastmod;
	// Polynomial representation of skew operator modulo modulus
	ZZ_pX polyRep;
	ZZ_pEX pRep;
	// Polynomial representation of the inverse skew operator
	ZZ_pX polyInvRep;
	ZZ_pEX piRep;

};

void promote(ZZ_pEX& target, const ZZ_pX& source);
