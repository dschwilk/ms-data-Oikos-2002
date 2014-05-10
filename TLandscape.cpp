/* File: TLandscape.cpp
 * ---------------------
 * This class handles the landscape grid and all spatial aspects of the simulation
 *
 */

#include "TLandscape.h"
#include "TGenotype.h"
#include "dws_random.h"
#include "dws_math.h"
#include "rec_defs.h"
#include <map>


TLandscape::TLandscape(unsigned gridsize)
:
 _gSize(gridsize),
 _deadSet(),
 _nSeeds(N_SEEDS),
  _burnNeighborhood(BURN_NEIGHBORHOOD),
  _seedNeighborhood(SEED_NEIGHBORHOOD),
  _pollenNeighborhood(POLLEN_NEIGHBORHOOD),
  _optimalGenotype(N_ALLELES, 0.5),
  _nAlleles(N_ALLELES),
  _envChangePeriod(ENV_CHANGE_PERIOD),
  _nEnvChangeLoci(N_ENV_CHANGE_LOCI),
  _age(0),
  _meanFit(0),
  _torchProp(0),
  _mortalityRate(MORTALITY_RATE),
  _mutationRate(MUTATION_RATE),
  _r(_r_),
  _R(_R_),
  _lightningFreq(LIGHTNING),
  _lightningBurn(false)
{
	 _gridFull=false;
	 _meanFit=0;
	 _torchProp=0;
	 FillGrid(_gSize);
 	_optimalGenotype = TGenotype(_nAlleles, 0.5);
}


TLandscape::TLandscape(const TLandscape& copy) 
: _gSize(copy._gSize),
 _deadSet(copy._deadSet),
 _nSeeds(copy._nSeeds),
 _nAlleles(copy._nAlleles),
 _r(copy._r),
 _R(copy._R),
 _burnNeighborhood(copy._burnNeighborhood),
 _seedNeighborhood(copy._seedNeighborhood),
 _pollenNeighborhood(copy._pollenNeighborhood),
 _lightningFreq(copy._lightningFreq),
 _lightningBurn(copy._lightningBurn),
 _gridFull(copy._gridFull),
 _nEnvChangeLoci(copy._nEnvChangeLoci)
{
	if(copy._gridFull) {
		_grid = new TPlant*[_gSize];
		for(int x=0;x<_gSize;x++)  {
		 _grid[x]= new TPlant[_gSize];
		 for(int y=0;y<_gSize;y++) {
			_grid[x][y] = copy._grid[x][y];
		 }
		}
	}

 
 _optimalGenotype = copy._optimalGenotype;
 _envChangePeriod = copy._envChangePeriod;
 _mortalityRate = copy._mortalityRate;
 _mutationRate = copy._mutationRate;
}


// Destructor
TLandscape::~TLandscape() 
{
	if(_gridFull) {
		for(int x=0;x<_gSize;x++) delete[] _grid[x];
		delete[] _grid;
	}
}



void TLandscape::FillGrid(int newsize)
{
	if(_gridFull){
		for(int x=0;x<_gSize;x++) delete[] _grid[x];
		delete[] _grid;
	}
	_gSize = newsize;
	 double fit;
	_grid = new TPlant*[_gSize];
	for(int x=0;x<_gSize;x++) {
		_grid[x] = new TPlant[_gSize];
		for(int y=0;y<_gSize;y++) {
			_grid[x][y] = TPlant(INITIAL_TORCH_PROP, _nAlleles, 0.5);
			fit = _grid[x][y].TestFitness(_optimalGenotype);
			_fitnessMap[Pos(x,y)] = fit;
			_meanFit+=fit;
			if(_grid[x][y].Flammability()==TORCH) _torchProp++;
			_grid[x][y].SetBirthTime(0);
		}
	}
	_meanFit = _meanFit/(_gSize*_gSize);
	_torchProp = _torchProp/(_gSize*_gSize);
	_gridFull=true;
}

/////////////////////////////////////////////////////////////////////////

/* Burn
 * ----
 * returns number of torches and burned torches and neighbors.
 */

void TLandscape::Burn()
{

	if(this->_lightningBurn) BurnFromLightning();
	else BurnNormal();
}

void TLandscape::BurnFromLightning()
{
	_meanFit=0;
	_torchProp=0;
	double fit;

	for(int x=0;x<_gSize;x++) {
		for(int y=0; y<_gSize;y++){
			Pos P(x,y);
			fit = _fitnessMap[P];
			_meanFit+=fit;
			if (_grid[x][y].Flammability()==TORCH) _torchProp++;  // for keeping stats on torch prop
			// Test for mortality
			if (DWS::random_chance(_mortalityRate)) {
				_deadSet.insert(P);
			}
			// Now test for lightning strike
			if (DWS::random_chance(_lightningFreq)) {
				_deadSet.insert(P);
				if (_grid[x][y].Flammability()==TORCH) {
					BurnNeighbors(P);
				}
			}
		}
	}
	_meanFit=_meanFit / (_gSize*_gSize);
	_torchProp=_torchProp / (_gSize*_gSize);
}

void TLandscape::BurnNormal()
{
	_meanFit=0;
	_torchProp=0;
	double fit;

	for(int x=0;x<_gSize;x++) {
		for(int y=0; y<_gSize;y++){
			Pos P(x,y);
			fit = _fitnessMap[P];
			_meanFit+=fit;
			if (DWS::random_chance(_mortalityRate)) 	_deadSet.insert(P);
			if (_grid[x][y].Flammability()==TORCH) {
				_torchProp++;
				GetNeighbors(P, _burnNeighborhood, _deadSet);
			}
		}
	}
	_meanFit=_meanFit / (_gSize*_gSize);
	_torchProp=_torchProp / (_gSize*_gSize);
}

void TLandscape::Recruit()
{
	using namespace DWS;
	std::map<Pos, TPlant> recruitMap;
	Pos m,f;
	TPlant newSeed,bestSeed;
	double fBest, fNew;
	for(std::set<Pos>::iterator i = _deadSet.begin();i!=_deadSet.end();i++) {
		//  do N_SEEDS for each cell
		// choose mother
		m=RandomNeighbor(*i,_seedNeighborhood);
		// choose father
		do {
			f=RandomNeighbor(m,_pollenNeighborhood);
		} while( f== m);

		TPlant::Mate(PlantAt(f),PlantAt(m), this->_R, this->_r, bestSeed);
		fBest=bestSeed.TestFitness(_optimalGenotype);
		for(int n=0;n<_nSeeds-1;n++) {  // N_SEEDS-1 more
			// choose mother
			m=RandomNeighbor(*i,_seedNeighborhood);
			// choose father
			f=RandomNeighbor(m,_pollenNeighborhood);
			TPlant::Mate(PlantAt(f),PlantAt(m), this->_R, this->_r, newSeed);
			newSeed.Mutation(this->_mutationRate);
			fNew = newSeed.TestFitness(_optimalGenotype);
			if(fNew > fBest) {
				fBest=fNew;
				bestSeed=newSeed;
			}
		}
		recruitMap[*i] = bestSeed;
	}

	// Now place those seeds on grid
	for(std::map<Pos, TPlant>::iterator mi = recruitMap.begin();mi!=recruitMap.end();mi++) {
		PlantAt((*mi).first) = (*mi).second;
		_fitnessMap[(*mi).first] = PlantAt((*mi).first).TestFitness(_optimalGenotype);
		PlantAt((*mi).first).SetBirthTime(this->Age());
	}
	this->_deadSet.erase(this->_deadSet.begin(), this->_deadSet.end());
}


void TLandscape::ChangeEnvironment()
{
	_optimalGenotype.MutateNBits(_nEnvChangeLoci, MUTATE_FLIP);
	for(int x=0;x<_gSize;x++) {
		for(int y=0; y<_gSize;y++){
			Pos p(x,y);
			_fitnessMap[p] = PlantAt(p).TestFitness(_optimalGenotype);
		}
	}
}

			
void TLandscape::RunForCycles(int n)
{
	for(int i=0;i<n;i++) {
		if(_age++ % _envChangePeriod == 0) ChangeEnvironment();
		this->Burn();
		this->Recruit();
	}
}


/* BurnNeighbors
 * --------------
 * Recursive fire-spread procedure used by burnFromLighting
 */
void TLandscape::BurnNeighbors(Pos p)
{
	using namespace DWS;
	std::set<Pos> temp;
	_deadSet.insert(p); // so we don't go in endless loop

	GetNeighbors(p, _burnNeighborhood, temp);
	for(std::set<Pos>::iterator i = temp.begin();i!=temp.end();i++) {
		if (PlantAt(*i).Flammability() == TORCH && (_deadSet.find(*i)==_deadSet.end() )){
			BurnNeighbors(*i);  // call recursively
		}
	}
	_deadSet.insert(temp.begin(), temp.end()); // add this list to the deadSet
}



/* GetNeighbors 
 * ------------
 * produces inclusive neighborhood of cell at position p for a torus landscape.  This includes cell at p itself.
 */
void TLandscape::GetNeighbors(Pos p, int neighborhood, std::set<Pos>& outSet) const
{
	using namespace DWS;
	int temp(0);
	int ox,oy;
	for(int x = p._x - neighborhood; x<= p._x+neighborhood; x++) {
		for(int y= p._y-neighborhood; y <= p._y+neighborhood; y++) {
			ox=(x+_gSize)%_gSize; // mod out for torus
			oy=(y+_gSize)%_gSize;
			outSet.insert(Pos(ox,oy));
			temp++;
		}
	}
}


Pos TLandscape::RandomNeighbor(Pos p, int neighborhood) const
{
	using namespace DWS;
	int x,y;
	x = random_integer(p._x-neighborhood,p._x+neighborhood);
	y = random_integer(p._y-neighborhood, p._y+neighborhood);
	x=(x+_gSize)%_gSize; // mod out for torus
	y=(y+_gSize)%_gSize;
	return Pos(x,y);
}
/* ChooseWeighted
 * --------------
 * Chooses from a set of positions according to their weights
 */

Pos ChooseWeighted(std::map<Pos, double> weightMap)
{
	double line=0;
	std::map<double, Pos> numLineMap;

	// create number line
	for(std::map<Pos, double>::iterator i = weightMap.begin(); i!=weightMap.end();i++) {
		line+=(*i).second;
		numLineMap[line] = (*i).first;
	}
	return (*(numLineMap.lower_bound(DWS::random_real(0,line)) )).second;
}




double TLandscape::SampleSimilarity(int nSamples) const
{
	using namespace DWS;
	double s=0;

	for(int i=0;i<nSamples;i++) {
		s += _grid[random_integer(0,_gSize-1)][random_integer(0,_gSize-1)].TestFitness(_grid[random_integer(0,_gSize-1)][random_integer(0,_gSize-1)]);
	}
	return s/nSamples;
}


void TLandscape::SampleAgeStructure(DWS::TSample& outS) const
{
	using namespace DWS;
	outS.Erase();
	for(int x=0;x<_gSize;x++) {
		for(int y=0; y<_gSize;y++){
			outS.Add(this->Age() - _grid[x][y].GetBirthTime());
		}
	}
}

void TLandscape::TorchDampFitness(DWS::TSample& tfit, DWS::TSample& dfit) const
{
	using namespace DWS;
	tfit.Erase();
	dfit.Erase();
	for(int x=0;x<_gSize;x++) {
		for(int y=0; y<_gSize;y++){
			if(_grid[x][y].Flammability()==TORCH) tfit.Add(FitnessAt(Pos(x,y)));
			else dfit.Add(FitnessAt(Pos(x,y)));
		}
	}
	if(dfit.N() < 1) dfit.Add(0);
	if(tfit.N() < 1) tfit.Add(0);
}

void TLandscape::Get1DSnapshot(DWS::TBivSample& FitFlamOutS) const
{
	using namespace DWS;
	FitFlamOutS.Erase();
	for(int x=0;x<_gSize;x++) {
		FitFlamOutS.Add(_grid[x][0].TestFitness(_optimalGenotype), _grid[x][0].Flammability()==TORCH?1:0);
	}
}

void TLandscape::Print(){
	for(int x=0;x<_gSize;x++) {
		for(int y=0; y<_gSize;y++){
			std::cout << (_grid[x][y].Flammability()==TORCH?'*':'~');
		}
		std::cout<< std::endl;
	}
}

int TLandscape::AddTorches(double p)
{
	int count=0;
	for(int x=0;x<_gSize;x++) {
		for(int y=0; y<_gSize;y++){
			if(_grid[x][y].Flammability()!=TORCH && DWS::random_chance(p)) {
				_grid[x][y].ChangeFlammability(TORCH);
				count++;
			}
		}
	}
	return count;
}


std::ostream& TLandscape::WriteFitnessGrid(std::ostream& os) const
{
	for(int x=0;x<_gSize;x++) {
		for(int y=0; y<_gSize;y++){
			os << (_grid[x][y].TestFitness(_optimalGenotype));
			os << '\n';
		}
	}
	return os;
}


std::ostream& TLandscape::WriteFlamGrid(std::ostream& os) const
{
	for(int x=0;x<_gSize;x++) {
		for(int y=0; y<_gSize;y++){
			os << (_grid[x][y].Flammability()==TORCH)?1:0;
			os << '\n';
		}
	}
	return os;
}


///////////////////////////////////////////////////////////////////////////////
// friend functions
///////////////////////////////////////////////////////////////////////////////
std::ostream& operator<<(std::ostream& os, const TLandscape& la)
{
	os << la._gSize << '\t';
	for(int y=0;y<la._gSize;y++) {
		for(int x=0;x<la._gSize;x++) {
			os << la._grid[x][y] << '\n';
		}
	}
	return os;
}
