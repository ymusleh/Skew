#pragma once
#include "SkewRing.h"
#include <vector>
#include <NTL/vector.h>
#include <cmath>

using namespace NTL;

class SkewRing;

class Skewpoly
{
public:
	Skewpoly();
	Skewpoly(SkewRing* ring);
	Skewpoly(int degree, SkewRing* ring);
	Skewpoly operator+(const Skewpoly& rhs);
	Skewpoly operator-(const Skewpoly& rhs);
	Skewpoly operator*(const Skewpoly& rhs);
	void setCoeff(const int pos, const ZZ_pE& val);
	void addCoeff(const int pos, const ZZ_pE& val);
	ZZ_pE getCoeff(const int pos) const;
	int deg() const;
	SkewRing* getRing() const;
	void setRing(SkewRing* ring);
	void print(std::ostream& out);
	void add(const Skewpoly& rhs);
	void leftScalarMult(const ZZ_pE& scalar);
	void rightScalarMult(const ZZ_pE& scalar);
	// equivalent to multiplying on the right by X^i
	void rightShift(const int i);
	void applySkewRel(ZZ_pE& ret, const ZZ_pE& val);
	//void fastMult(Skewpoly& res, Skewpoly& lhs, Skewpoly& rhs);

private:
	int degree;
	SkewRing* ring;
	//Vec<ZZ_pE> coeffs;
	std::vector<ZZ_pE> coeffs;	

};

//Skewpoly operator+(const Skewpoly& lhs, const Skewpoly& rhs);
//Skewpoly operator-(const Skewpoly& lhs, const Skewpoly& rhs);
//Skewpoly operator*(const Skewpoly& lhs, const Skewpoly& rhs);
Skewpoly operator*(const ZZ_pE& scalar, const Skewpoly& rhs);
Skewpoly naiveMult(const Skewpoly& lhs, const Skewpoly& rhs);
void fastMult(Skewpoly& ret, const Skewpoly& lhs, const Skewpoly& rhs);



