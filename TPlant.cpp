/* File: TPlant.cpp
 *-----------------
 */
#include "TPlant.h"
#include "dws_random.h"

// Constructors
TPlant::TPlant()
{}

TPlant::TPlant(const TPlant& copy)
: _flammability(copy._flammability),
  _fitness(copy._fitness)
{}

TPlant::TPlant(double fp /*torch prob*/, int n, double p)
:  _fitness(n,p)
{
	_flammability=DWS::random_chance(fp)?TORCH:DAMP;
}

TPlant::TPlant(const TPlant& f, const TPlant& m, double R /*flam-fit recomb rate*/, double r /*per locus recomb rate*/)
{
	TGenotype::CreateOffspring(f._fitness,m._fitness,r,this->_fitness);
	if(this->_fitness.BitAt(0)==f._fitness.BitAt(0)) { // started with father
		this->_flammability=DWS::random_chance(R)?m._flammability:f._flammability;
	} else { // started with mother
		this->_flammability=DWS::random_chance(R)?f._flammability:m._flammability;
	}
}

TPlant::~TPlant()
{}


TPlant& TPlant::operator=(const TPlant& other)
{
	if(&other == this) return *this;
	_fitness = other._fitness;
	_flammability=other._flammability;
	return (*this);
}



double TPlant::TestFitness(const TGenotype compare) const
{
	return (this->_fitness.Compare(compare) / (double)this->_fitness.NAlleles());
}

double TPlant::TestFitness(const TPlant compare) const
{
	return (this->_fitness.Compare(compare._fitness) / (double)this->_fitness.NAlleles());
}


// static functions
void TPlant::Mate(const TPlant& f, const TPlant& m, double R /*flam-fit recomb rate*/, 
				  double r /*per locus recomb rate*/, TPlant& offspring)
{
	TGenotype::CreateOffspring(f._fitness,m._fitness,r,offspring._fitness);
	if(offspring._fitness.BitAt(0)==f._fitness.BitAt(0)) { // started with father
		offspring._flammability=DWS::random_chance(R)?m._flammability:f._flammability;
	} else { // started with mother
		offspring._flammability=DWS::random_chance(R)?f._flammability:m._flammability;
	}
}


// Friend functions


std::ostream& operator<<(std::ostream& os, const TPlant& g)
{
	os << "FLAMMABILITY=" << (g._flammability==TORCH?"TORCH":"DAMP") << " GENOTYPE=" << g._fitness << std::endl;
	return os;
}
