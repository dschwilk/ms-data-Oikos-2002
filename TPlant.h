/* file :TPlant.h
 * --------------
 */

#ifndef _plant_h
#define _plant_h

#include "TGenotype.h"


enum FlammabilityT {DAMP, TORCH};

class TPlant {
public:
	TPlant();
	TPlant(const TPlant& copy);
	TPlant(double fp /*torch prob*/, int n, double p); // random TPlant
	TPlant(const TPlant& f, const TPlant& m, double R /*flam-fit recomb rate*/, double r /*per locus recomb rate*/); // create torch or damp with random genome
	~TPlant();
	TPlant& operator=(const TPlant& other);

	FlammabilityT Flammability() const {return _flammability;}
	void ChangeFlammability(FlammabilityT t) {_flammability=t;}
	void Mutation(double p) {this->_fitness.Mutation(p);}

	int GetBirthTime() const {return _birthage;} 
	void SetBirthTime(int t) {_birthage = t;};
	double TestFitness(const TGenotype compare) const;
	double TestFitness(const TPlant compare) const;
	//
	friend std::istream& operator>>(std::istream& is, TPlant& g);
	friend std::ostream& operator<<(std::ostream& os, const TPlant& g);

	//
	static void Mate(const TPlant& f, const TPlant& m, double R /*flam-fit recomb rate*/, 
				  double r /*per locus recomb rate*/, TPlant& offspring);

private:
	TGenotype _fitness;
	FlammabilityT _flammability;
	int _birthage; // birth time step
}; // class TPlant


#endif  // TPlant.h
