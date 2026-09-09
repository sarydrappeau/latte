#ifndef NEWINTEGRATION_H
#define NEWINTEGRATION_H
#include <NTL/vec_vec_ZZ.h>
#include <NTL/vec_ZZ.h>
#include <NTL/ZZ.h>
#include "PolyTrie.h"
#include "PolyRep.h"
//#include "PolyRep.cpp"
//#include "PolyTrie.cpp"
#include "residue.h"
#include "cone.h"
#include "LattException.h"

NTL_CLIENT

struct simplexZZ
{

	int d;			// dimension of the space. s.length = d + 1
	vec_vec_ZZ s;	// s[i] = vector of the ith vertex
	ZZ v; 			//volume of the parallelepiped. We take v = det(simplex rays) and we do NOT divide by d!.

	void print(ostream & out)
	{
		out << "d = " << d << endl;
		out << "v = " << v << endl;
		int i;
		for( i = 0; i < s.length(); ++i)
		{
			out << "s[" << i << "] = ";
			for(int k = 0; k < s[i].length(); ++k)
				out << s[i][k] << ", ";
			out << endl;
		}

	}//pirnt
};

/**
 * The following structure is relevant for the "polytope triangulation" method.
 * It holds a cache for the numerical values involved in the formula which depend only on
 * the simplex and the linear form coefficients direction, but not on the power of the linear form.
 */
struct SimplexLinFormGeometry
{
	//ASSUMES the polytope has dimension less than maxDimension: the flags below
	//sit in a fixed array rather than a vector, so that filling them in does not
	//go to the heap.
	static const int maxDimension = 1000;

	vec_ZZ innerProduct;		// innerProduct[i] = <l, s_i>
	vec_ZZ denominator;			// denominator[i] = \prod_{j != i} <l, s_i - s_j>; zero when a residue is needed.
	bool repeated[maxDimension];// repeated[i] = 1 iff <l, s_i> == <l, s_j> for some j < i.
	ZZ lcmOfDenominators;		// lcm of the non-zero denominators above.

	void compute(const vec_ZZ &l, const simplexZZ &mySimplex);
};

/**
 * Hold a cache for the values (d + m)! used in the computation.
 */
class DegreeFactorials
{
public:
	DegreeFactorials(): dimension(-1) {}

	const ZZ & get(int d, int m);

private:
	int dimension;				//dimension the table was built for; -1 = no table yet.
	std::vector<ZZ> factorials;	//factorials[m] = (dimension + m)!
};

/**
 * True if the term's direction is the vector l.
 *
 * Linear forms are sorted by direction first and degree last (see
 * BurstTerm::lessThan), so upon iteration we handle every power of one
 * direction before moving to the next direction; the integration loops use this to
 * notice the change and to redo the work that depends on the direction alone.
 */
bool sameDirection(const term<RationalNTL, ZZ>* form, const vec_ZZ &l);
void copyDirection(const term<RationalNTL, ZZ>* form, vec_ZZ &l);

void update(ZZ &a, ZZ &b, const vec_ZZ &l, const simplexZZ &mySimplex,int m, const RationalNTL &coe, const ZZ &de);
void update(ZZ &a, ZZ &b, const SimplexLinFormGeometry &geometry, const simplexZZ &mySimplex, int m, const RationalNTL &coe, const ZZ &de);
void delSpace(string &line);
void convertToSimplex(simplexZZ&, string);
void integrateLinFormSum(ZZ &a, ZZ &b, PolyIterator<RationalNTL, ZZ>* it, const simplexZZ &mySimplex);
RationalNTL integrateLinFormProducts(PolyIterator<RationalNTL, ZZ>* it, const simplexZZ &mySimplex, const int productCount);
void integrateMonomialSum(ZZ &numerator, ZZ &denominator, monomialSum &monomials, const simplexZZ &mySimplex);
void _integrateMonomialSum(ZZ &numerator, ZZ &denominator, _monomialSum &monomials, const simplexZZ &mySimplex);

template <class T>
class FormIntegrateConsumer : public FormSumConsumer<T> {
public:
  FormIntegrateConsumer() { }
  // Take linear form and consume it.
  void ConsumeLinForm(const RationalNTL& coefficient, int degree, const vec_ZZ& coefs);
  void setFormSum(const string& myForms) { linForms = myForms; }
  string getFormSum() { return linForms; }
  void setDimension(int dimension) { mySimplex->d = dimension; }//also stored in mySimplex
  int getDimension() { return mySimplex->d; }
  void setSimplex(simplexZZ& simplex) { mySimplex = &simplex; numerator = to_ZZ(0); denominator = to_ZZ(0); }
  void getResults(ZZ& num, ZZ& den);
  ~FormIntegrateConsumer() {}
private:
  string linForms;
  simplexZZ* mySimplex;
  ZZ numerator, denominator;
};

template <class T>
void FormIntegrateConsumer<T>::ConsumeLinForm(const RationalNTL& coefficient, int degree, const vec_ZZ& coefs)
{
	ZZ de = to_ZZ(1);
	for (int i=1;i<=mySimplex->d+degree;i++)
	{
		de=de*i;
	};
	update(numerator, denominator, coefs, *mySimplex, degree, coefficient, de);
}

template <class T>
void FormIntegrateConsumer<T>::getResults(ZZ& num, ZZ& den)
{
	if (denominator < 0)
	{
		num = to_ZZ(-1) * numerator; den = to_ZZ(-1) * denominator;
	}
	else
	{
		num = numerator; den = denominator;
	}
}

#endif
