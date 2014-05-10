/* File: TLandscape.h
 * ------------------
 */

#ifndef _TLandscape_h_
#define _TLandscape_h_

#include "TPlant.h"
#include "dws_math.h"
#include "dws_stats.h"
#include <iostream>
#include <set>
#include <map>

typedef DWS::TPos<int> Pos;


class TLandscape {
public:
	TLandscape(unsigned gridsize);
	TLandscape(const TLandscape& copy);
	TLandscape operator=(const TLandscape& source);
	~TLandscape();

	//
	void GetNeighbors(Pos p, int neighborhood, std::set<Pos>& outSet) const;
	void BurnNeighbors(Pos p);
	Pos RandomNeighbor(Pos p, int neighborhood) const;
	int AddTorches(double p);
	int GridSize() const {return _gSize;}
	int Age() const {return _age;}
	void Setr(double r) {_r=r;}
	void SetChangePeriod(int n) {_envChangePeriod = n;}
	void SetChangeLoci(int n) {_nEnvChangeLoci = n;}
	void SetMortalityRate(double m) {this->_mortalityRate = m;}
	void SetMutationRate(double mu) {this->_mutationRate = mu;}
	void SetPollenNeighborhood(int n) {this->_pollenNeighborhood = n;}
	void SetSeedNeighborhood(int s) {this->_seedNeighborhood = s;}
	void SetBurnNeighborhood(int b) {this->_burnNeighborhood=b;}
	void RunForCycles(int n);

	void Print();


	TPlant& PlantAt(Pos p) {return (_grid[p._x][p._y]);}
	double FitnessAt(Pos p) const {return (_fitnessMap.find(p))->second;}
	int PlantAgeAt(Pos p) {return (this->Age() - _grid[p._x][p._y].GetBirthTime());}
	double MeanFitness() {return _meanFit;}
	double TorchProp() {return _torchProp;}
	double SampleSimilarity(int nSamples) const;
	void SampleAgeStructure(DWS::TSample& outS) const;
	void Get1DSnapshot(DWS::TBivSample& FitFlamOutS) const;
	void TorchDampFitness(DWS::TSample& tfit, DWS::TSample& dfit) const;


	std::ostream& WriteFitnessGrid(std::ostream& os) const;
	std::ostream& WriteFlamGrid(std::ostream& os) const;
	friend std::ostream& operator<<(std::ostream& os, const TLandscape& la);
//	friend std::istream& operator>>(std::istream& is, TLandscape& la);
private:
	void BurnFromLightning();
	void BurnNormal();

protected:

	void Burn(); 
	void Recruit();
	void ChangeEnvironment();
	void FillGrid(int newsize);

private:
	bool _gridFull;
	int _gSize;
	TPlant** _grid;  // array of plants
	std::set<Pos> _deadSet;
	std::map<Pos, double> _fitnessMap;
	TGenotype _optimalGenotype;

protected:
	int _nSeeds;
	int _nAlleles;
	bool _lightningBurn;


	int _burnNeighborhood;
	int _seedNeighborhood;
	int _pollenNeighborhood;

	double _meanFit;
	double _torchProp;
	double _mortalityRate;
	double _mutationRate;
	double _r; // recomb rate between fit loci
	double _R;
	double _lightningFreq;

	int _envChangePeriod;
	int _nEnvChangeLoci; // number of loci to change each change period
	int _age;

}; // class TLandscape


#endif // TLandscape.h
