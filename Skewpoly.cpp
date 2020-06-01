#include "Skewpoly.h"


Skewpoly Skewpoly::operator+(const Skewpoly& rhs) {

	//if ( !(*(this->ring) == *(rhs.ring)) ) {
	//	return NULL;
	//}

	Skewpoly sumpoly(this->ring);
	//sumpoly.ring = this->ring;
	sumpoly.degree = std::max(this->degree, rhs.degree);
	sumpoly.coeffs.resize(sumpoly.degree + 1);
	//std::vector<ZZ_pE>::iterator i;
	for (int i = 0; i <= sumpoly.degree; ++i) {
		sumpoly.coeffs[i] = ZZ_pE(0);
		if ( i <= this->degree ) {
			sumpoly.coeffs[i] += this->coeffs[i];
		}
		if ( i <= rhs.degree ) {
			sumpoly.coeffs[i] += rhs.coeffs[i];
		}
	}
	// Remove leading zeros
	for (int i = sumpoly.degree; i > 0; --i) {
		if (IsZero(sumpoly.coeffs[i])) {
			sumpoly.coeffs.pop_back();
			sumpoly.degree--;
		} else return sumpoly;
	}
	return sumpoly;

}

// Sets to coefficient at pos to val; if pos > degree, fills in the zeros
void Skewpoly::setCoeff(const int pos, const ZZ_pE &val) {
	if (this->degree < pos) {
		this->coeffs.resize(pos + 1);
		for (int i = this->degree + 1; i < pos; ++i) {
			this->coeffs[i] = ZZ_pE::zero();
		}
		this->degree = pos;
	}
	this->coeffs[pos] = val;
}

void Skewpoly::addCoeff(const int pos, const ZZ_pE& val) {
	if (this->degree >= pos) {
		this->coeffs[pos] += val;
	}
	else this->setCoeff(pos, val);
}

ZZ_pE Skewpoly::getCoeff(const int pos) const {
	if (pos <= this->degree) {
		return this->coeffs[pos];
	} else return ZZ_pE::zero();
}

int Skewpoly::deg() const {
	return this->degree;
}

SkewRing* Skewpoly::getRing() const {
	return this->ring;
}

void Skewpoly::print(std::ostream& out) {

	for (int i = 0; i <= this->degree; ++i) {
		out << this->coeffs[i] << " ";
	}
	out << " end" << std::endl;
	

}

void Skewpoly::leftScalarMult(const ZZ_pE &scalar) {
	for (int i = 0; i <= degree; ++i) {
		this->coeffs[i] = this->coeffs[i]*scalar;
	}
}

void Skewpoly::rightScalarMult(const ZZ_pE& scalar) {
	ZZ_pE temp, accum = scalar;
	//ZZ_pE temp = scalar;
	for (int i = 0; i <= degree; ++i) {
		this->coeffs[i] = accum*this->coeffs[i];
		this->applySkewRel(temp, accum);
		accum = temp;
	}

}



void Skewpoly::add(const Skewpoly& rhs) {
	int ndegree = std::max(this->degree, rhs.degree);
	this->coeffs.resize(ndegree + 1);
	for (int i = 0; i <= ndegree; ++i) {
		if (i > this->degree) this->coeffs[i] = ZZ_pE::zero();
		if (i <= rhs.degree) this->coeffs[i] += rhs.coeffs[i];
	}
        // Remove leading zeros
        for (int i = ndegree; i > 0; --i) {
                if (IsZero(this->coeffs[i])) {
                        this->coeffs.pop_back();
                        ndegree--;
                } else break;
        }

	this->degree = ndegree;
		
}

void Skewpoly::applySkewRel(ZZ_pE& rel, const ZZ_pE& val) {
	this->ring->applySkewRel(rel, val);
}


void Skewpoly::rightShift(const int i) {
	Skewpoly temp(this->degree + i, this->ring);
	for (int j = 0; j <= this->degree; ++j) {
		temp.coeffs[j + i] = this->coeffs[j];
	}
	*this = temp; 
}
 

Skewpoly operator*(const ZZ_pE& scalar, const Skewpoly& rhs) {
	Skewpoly result = rhs;
	result.leftScalarMult(scalar);
}


/*
 *  Non-method operators
 */




Skewpoly naiveMult(const Skewpoly& lhs, const Skewpoly& rhs) {
	Skewpoly acc(lhs.deg() + rhs.deg(), lhs.getRing());
	Skewpoly res(lhs.deg(), lhs.getRing());
	for (int i = 0; i <= rhs.deg(); ++i) {
		res = lhs;
		res.rightScalarMult(rhs.getCoeff(i));
		res.rightShift(i);
		acc.add(res);
	}
	return acc;
	 	
}

void fastMult(Skewpoly& ret, const Skewpoly& lhs, const Skewpoly& rhs) {
	int s = std::max(lhs.deg(), rhs.deg());
	int sstar = (int)ceil(sqrt(s + 1));
	mat_ZZ_pE A, B, C;
	ZZ_pE temp, temp2;
	A.SetDims(sstar, sstar);
	B.SetDims(sstar, s + sstar);
	C.SetDims(sstar, s + sstar);
	for (int i = 0; i < sstar; ++i) {
		for (int j = 0; j < sstar; ++j) {
			lhs.getRing()->powerSkewRel(temp, lhs.getCoeff(i*sstar + j),-1*i*sstar);
			A[i][j] = temp;
		}
		for (int j = 0; j < s + sstar; ++j) {
			
			if ( ((j - i) >= 0) && ((j - i) <= s) ) {
				rhs.getRing()->powerSkewRel(temp, rhs.getCoeff(j - i), i);
				B[i][j] = temp;
			} else {
				B[i][j] = ZZ_pE::zero();
			}
		}
		
	}
	mul(C, A, B);
	for (int i = 0; i < sstar; ++i) {
		for (int j = 0; j < s + sstar; ++j) {
			temp2 = C[i][j];
			lhs.getRing()->powerSkewRel(temp, temp2, i*sstar);
			ret.addCoeff(i*sstar + j, temp);
		}
	}
}




/*
 *  Constructors
 */

Skewpoly::Skewpoly() : degree(0), coeffs({ZZ_pE::zero()}), ring(NULL) {
}

Skewpoly::Skewpoly(SkewRing* ring) : degree(0), coeffs({ZZ_pE::zero()}), ring(ring) {}

Skewpoly::Skewpoly(int degree, SkewRing* ring) : degree(degree), ring(ring) {
	this->coeffs.resize(degree + 1);
	for (int i = 0; i <= this->degree; ++i) {
		this->coeffs[i] = ZZ_pE::zero();
	}
}
