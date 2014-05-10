/* file: TGenotype.cpp
 * -------------------
 * Implementation of the TGenotype class
 */
 
#include <algorithm>
#include <set>
#include <assert.h>

#include "TGenotype.h"
#include "dws_random.h"
#include "rec_defs.h"

// static members

// Constructors
TGenotype::TGenotype(){}
TGenotype::TGenotype(const TGenotype& copy) : _bitVector(copy._bitVector) {}
TGenotype::TGenotype(int n, double p)
{
	for(int i=0;i<n;i++) {
		if(DWS::random_chance(p)) {
			this->_bitVector.push_back(1);
		} else {
			this->_bitVector.push_back(0);
		}
	}
}

// Destructor
TGenotype::~TGenotype() {}

// operators
TGenotype& TGenotype::operator=(const TGenotype& other)
{
	if(&other == this) return *this;
	_bitVector = other._bitVector;
	return (*this);
}

// Methods

int TGenotype::NAlleles() const
{
	return this->_bitVector.size();
}

int TGenotype::Compare(const TGenotype& other) const
{
	int count(0);
	int s = this->NAlleles();
	for(int i=0;i<s;i++) if(this->_bitVector[i]==other._bitVector[i]) count++;
	return count;
}


void TGenotype::Mutation(double p)
{
	for(BitVectorT::iterator i = _bitVector.begin(); i!= _bitVector.end();i++){
		if(DWS::random_chance(p)) *i = *i==0?1:0;
	}
}

void TGenotype::MutateNBits(int n, bool flip)
{
	std::set<int> lociSet;
	int temp;
	// get array of loci to mutate
	for(int i=0;i<n;i++) {
		do {
			temp=DWS::random_integer(0,this->NAlleles()-1);
		} while(lociSet.find(temp) != lociSet.end() );
		lociSet.insert(temp);
	}
	// now flip or randomize bits according to flag "flip"
	for(std::set<int>::iterator j = lociSet.begin();j!= lociSet.end();j++) {
		if(flip) _bitVector[*j] = _bitVector[*j]==0?1:0;
		else _bitVector[*j] = DWS::random_bit();
	}
}

//////////////////////////////////////////////////////////////////////////////
/* CreateOffspring()
 * -----------------
 * Static function for creating haploid offspring from two parents.  Caller provides number
 * of crossover events and function distributes them uniformly on genome
 */

void TGenotype::CreateOffspring(const TGenotype& f, const TGenotype& m, double r /*per locus recombination rate*/, TGenotype& result)
{
	using namespace DWS;

#ifdef _DEBUG
	assert(f.NAlleles==m.NAlleles); // should always mate like types!
#endif

	result._bitVector.erase(result._bitVector.begin(), result._bitVector.end());
	result._bitVector.reserve(f.NAlleles());
	const TGenotype * p1, * p2;

	// choose starting parent
	if(random_chance(0.5)) {
		p1 = &f; p2 = &m;
	} else {
		p1=&m; p2=&f;
	}

	// add to offspring switching parents when crossover events occur
	for(int pos=0;pos<f.NAlleles();pos++) {
		if(random_chance(r)) swap(p1 , p2);
		result._bitVector.push_back( p1->_bitVector[pos]);
	}
}

// Extraction for comparison
int TGenotype::BitAt(int i) const
{
	return this->_bitVector[i];
}


// Friends
std::istream& operator>>(std::istream& is, TGenotype& g)
{
	int n, x;
	is >> n;
	g._bitVector.erase(g._bitVector.begin(), g._bitVector.end());
	g._bitVector.reserve(n);
	for(int i=0;i<n;i++) { 
		is >> x ; 
		g._bitVector.push_back(x);
	}
	return is;
}

std::ostream& operator<<(std::ostream& os, const TGenotype& g)
{
	os << g._bitVector.size() << ' ';
	std::for_each(g._bitVector.begin(), g._bitVector.end(), DWS::print<int>(os));
	return os;
}

