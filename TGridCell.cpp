/* File TGridCell.cpp
 */

#include "TGridCell.h"
#include "dws_random.h"

#define RESPROUT_AGE 4

//ctors
TGridCell::TGridCell()
:
	_plant(0)
{}

TGridCell::TGridCell(const TGridCell& copy)
{
		if (_plant==0) _plant = new TGenotype;
		*_plant = *(copy._plant);
}


TGridCell::~TGridCell()
{
	delete _plant;
}

	//
TGridCell& TGridCell::operator=(const TGridCell& source) 
{
	if(&source == this) return *this;
	if (_plant==0) _plant = new TGenotype;
	*_plant = *(source._plant);
	return *this;
}

/////////////////////////////////////////////////////////////////
//  Friend fxns
/////////////////////////////////////////////////////////////////

std::ostream& operator<<(std::ostream& os, const TGridCell& g)
{
	os << *g._plant << '\n';
	return os;
}

std::istream& operator>>(std::istream& is, TGridCell& g)
{
	is >> *g._plant;
	return is;
}
