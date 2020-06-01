#include "SkewRing.h"

bool SkewRing::operator==(const SkewRing &rhs) {
	return  (this->chr == rhs.chr) && 
		(this->modulus == rhs.modulus) &&
		(this->skewRelation == rhs.skewRelation);
}

// Initialize the base field (ZZ_p) with a prime power characteristic
void SkewRing::setBaseFieldChar(const long characteristic) {
	this->chr = characteristic;
	ZZ_p::init(ZZ(this->chr));
}

// Given an irreducible polynomial, initialize the extension field
void SkewRing::setBaseField(ZZ_pX &modulus) {
	this->modulus = modulus;
	ZZ_pE::init(this->modulus);
}


// Given the extension degree, initialize the extension field
void SkewRing::setBaseField(int degree) {
	NTL::BuildIrred(this->modulus, degree);
	ZZ_pE::init(this->modulus);
}

long SkewRing::baseFieldChar() {
	return this->chr;
}

ZZ_pX SkewRing::getMod() {
	return this->modulus;
}

void SkewRing::init() {
	ZZ_p::init(ZZ(this->chr));
	ZZ_pE::init(this->modulus);
}


// Applies the SkewRelation to a field element
void SkewRing::applySkewRel(ZZ_pE& ret, const ZZ_pE& val) {
	eval(ret, this->skewRelation, val);
	//power(ret, val, this->chr);
}

void SkewRing::fastSkewRel(ZZ_pE& ret, const ZZ_pE& val) {
	ZZ_pX repl = rep(val);
	ZZ_pX res = CompMod(repl, this->polyRep, this->fastmod);
	//eval(rel, this->pRep, val);
	ret = conv<ZZ_pE>(res);
}

void SkewRing::fastInvRel(ZZ_pE& ret, const ZZ_pE& val) {
	ZZ_pX repl = rep(val);
	ZZ_pX res = CompMod(repl, this->polyInvRep, this->fastmod);
	ret = conv<ZZ_pE>(res);
}

// computes \sigma^e(val)
void SkewRing::powerSkewRel(ZZ_pE& ret, const ZZ_pE& val, const long e) {
	long i = e;
	ZZ_pX pwr, repl = rep(val), res, comp;
	SetX(comp);
	if ( e == 0 ) {
		ret = val;
		return;
	}
	else if (e > 0) {
		pwr = this->polyRep;
	}
	else if (e < 0) {
		pwr = this->polyInvRep;
		i = -e;
	}
	// compute the modular representation of \sigma^e via repeated modular composition
	while (i > 1) {
		if (i % 2 == 0) {
			pwr = CompMod(pwr, pwr, this->fastmod);
			i /= 2;
		}
		else {
			comp = CompMod( pwr, comp, this->fastmod);
			pwr = CompMod( pwr, pwr, this->fastmod);
			i = (i - 1) / 2;
			
		}
	}
	pwr = CompMod(pwr, comp, this->fastmod);
	res = CompMod(repl, pwr, this->fastmod);
	ret = conv<ZZ_pE>(res);
}

// promotes a ZZ_pX to a ZZ_pEX
void promote(ZZ_pEX& target, const ZZ_pX& source) {
	for (int i = 0; i < deg(source); ++i) {
		SetCoeff(target, i, coeff(source, i));
	}
}


SkewRing::SkewRing(const long chr, const ZZ_pX &modulus, 
			const ZZ_pEX &skewRelation): chr(chr), modulus(modulus), skewRelation(skewRelation) {};

SkewRing::SkewRing(const long chr, const ZZ_pX &modulus): chr(chr), modulus(modulus) {
	SetCoeff(this->skewrel, this->chr);
	SetCoeff(this->skewRelation, this->chr);
	build(this->fastmod, this->modulus);
	// set up the polynomial representation x^q mod g, and its inverse
	rem(this->polyRep, this->skewrel, this->fastmod);
	PowerMod(this->polyInvRep, this->polyRep, power_long(chr, deg(modulus) - 2), this->fastmod); 
	promote(this->pRep, this->polyRep);
	promote(this->piRep, this->polyInvRep);
	
}
