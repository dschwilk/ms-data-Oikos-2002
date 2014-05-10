/* TGridCell.h
 * ------------
 * Interface for the grid cell class
 */

#ifndef _TGridCell_h
#define _TGridCell_h

#include "TGenotype.h"
#include <iostream>


class TGridCell {
public:
	TGridCell();
	TGridCell(const TGridCell& copy);
	~TGridCell();
	TGridCell & operator=(const TGridCell& source);
	TGenotype& Plant() {return *_plant;}

	friend std::ostream& operator<<(std::ostream& os, const TGridCell& g);
	friend std::istream& operator>>(std::istream& is, TGridCell& g);

private:
	TGenotype* _plant;  // set to zero when cell is empty
};


#endif

