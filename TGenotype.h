// file :TGenotype.h
// *****************
// Simple Genotype class for flammability and recombination simulation
// This is a simplified version of the more complete TGenotype class used 
// in serotiny simulation.

#ifndef _genotype_h
#define _genotype_h

#include <iostream>
#include <vector>

typedef std::vector<int> BitVectorT;

class TGenotype {
public:
	TGenotype();
	TGenotype(const TGenotype& copy);
	TGenotype(int n, double p); // creates bitvector of length n, with positions occupied with 1s in prob p)
	~TGenotype();
	TGenotype& operator=(const TGenotype& other);


	int NAlleles() const;
	void Mutation(double p); // give per locus mutation probability
	void MutateNBits(int n, bool flip); // mutates (random or flip according to "flip) flips n bits in bitstring
	int Compare(const TGenotype& other) const; // returns number of alleles that match between two individuals
	int& operator[](int i){return this->_bitVector[i];}
	int BitAt(int i) const;

	// Friend functions
	friend std::istream& operator>>(std::istream& is, TGenotype& g);
	friend std::ostream& operator<<(std::ostream& os, const TGenotype& g);

	// static functions
	static void CreateOffspring(const TGenotype& f, const TGenotype& m, 
		double r /* per locus recombination rate */, TGenotype& result);

private:
	BitVectorT _bitVector;
}; // class TGenotype


#endif  // TGenotype.h