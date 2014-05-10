/* File: dwstats.cpp
 * -----------------
 *
 *  Implementation of statistics class
 */

#include "dws_stats.h"
#include "dws_math.h"
#include "dws_random.h"

#include <stdexcept>
#include <algorithm>
#include <numeric>

namespace DWS {

/***************************************************************************************
 * Class TSample Implementation														   *
 * ----------------------------                                                        *
 ***************************************************************************************/
TSample::TSample()
:  _hasMean(false),
  _hasMedian(false),
  _hasVariance(false),
  _hasSOS(false),
  _mean(0),
  _median(0),
  _variance(0),
  _SOS(0)
{
}

TSample::TSample(int size)
: _hasMean(false),
  _hasMedian(false),
  _hasVariance(false),
  _hasSOS(false),
  _mean(0),
  _median(0),
  _variance(0),
  _SOS(0),
  _V(size)
{}

TSample::TSample(const sample_vector& inV)
: _hasMean(false),
  _hasMedian(false),
  _hasVariance(false),
  _hasSOS(false),
  _mean(0),
  _median(0),
  _variance(0),
  _SOS(0)
{
	for(std::vector<double>::const_iterator i = inV.begin(); i<inV.end();++i) {
		_V.push_back(*i);
		_range.Add(*i);
	}

}


TSample::TSample(const TSample& copy)
: _hasMean(copy._hasMean),
  _hasMedian(copy._hasMedian),
  _hasVariance(copy._hasVariance),
  _hasSOS(copy._hasSOS),
  _mean(copy._mean),
  _median(copy._median),
  _variance(copy._variance),
  _SOS(copy._SOS),
  _range(copy._range),
  _V(copy._V)
{}


TSample& TSample::operator=(const TSample& source)
{
	if(this==&source) return *this;
	_hasMean=source._hasMean;
	_hasMedian=source._hasMedian;
	_hasVariance=source._hasVariance;
	_hasSOS=source._hasSOS;
	_mean=source._mean;
	_median=source._median;
	_variance=source._variance;
	_SOS=source._SOS;
	_V=source._V;
	_range=source._range;
	return(*this);
}

// Destructor

TSample::~TSample()
{}

void TSample::Reset()
{
	_hasMean=_hasMedian=_hasSOS=_hasVariance=false;
}

/*********************************************************************************
 * Methods
 *********************************************************************************/

void TSample::Add(double x)
{
	Reset();
	_V.push_back(x);
	this->_range.Add(x);
}


double TSample::Mean()
{
	if(_hasMean) return(_mean);
	if(this->N() < 1) throw std::runtime_error(std::string("Empty sample!"));
	_mean = 0;
	for(std::vector<double>::const_iterator i = _V.begin(); i!= _V.end();) _mean += *(i++);
	_hasMean=true;
	return (_mean=_mean/double(this->N()));
}

double TSample::Median()
{
	if(_hasMedian) return(_median);
	if(this->N() < 1) throw std::runtime_error(std::string("Empty sample!")); 
	sample_vector tempV(_V);
	std::sort(tempV.begin(),tempV.end());
	_median = (N()%2 != 0) ? tempV[(N()/2)] : ((tempV[N()/2] + tempV[(N()/2)-1]) / 2.0);
	_hasMedian=true;
	return(_median);
}

double TSample::Variance()
{
	if(_hasVariance) return(_variance);
	_variance = this->SumOfSquares() / double(this->N()-1);
	_hasVariance=true;
	return(_variance);
}

inline double TSample::Deviate(double x) {return(x-Mean());}

double TSample::SumOfSquares()
{
	if (_hasSOS) return (_SOS);
	add_squared_diff sq(this->Mean());
	_SOS = std::accumulate(_V.begin(), _V.end(), double(0.0), sq);
	_hasSOS=true;
	if(_SOS<EPS) _SOS=0.0;
	return _SOS;
}
	
inline double TSample::SumValueSquares()
{
	return (std::accumulate(_V.begin(), _V.end(), double(0.0), add_square()));
}

void TSample::Sort() 
{
	std::sort(_V.begin(), _V.end());
}

void TSample::Shuffle()
{
	std::random_shuffle(_V.begin(), _V.end());
}

/******************************************************************************
 * Class TBivSample
 ******************************************************************************/

TBivSample::TBivSample(){}

TBivSample::TBivSample(sample_vector v1, sample_vector v2) 
{
	_xv = v1;
	_yv = v2;
}

TBivSample::~TBivSample()
{}

TSample& TBivSample::X()
{
	return(_xv);
} 

TSample& TBivSample::Y()
{
	return(_yv);
}

void TBivSample::Erase()
{
	_xv.Erase();
	_yv.Erase();
}
void TBivSample::Reserve(int n)
{
	_xv.Reserve(n);
	_yv.Reserve(n);
}

void TBivSample::AddObservation(double x, double y)
{
	_xv.Add(x);
	_yv.Add(y);
}

void TBivSample::Shuffle()
{
	_xv.Shuffle();
	_yv.Shuffle();
}


void TBivSample::Correlate(TCorrelation& outCorr, bool one_tail /*=false*/)
{
	outCorr.R( SumCrossProducts() / (sqrt ((_xv.SumOfSquares()*_yv.SumOfSquares()) ) ) );
	outCorr.V(_xv.N()-2);
	outCorr.Z(0.5*log((1.0+outCorr.R()+EPS)/(1.0-outCorr.R()+EPS)));
	outCorr.T(outCorr.R()*sqrt(outCorr.V()/((1.0-outCorr.R()+EPS)*(1.0+outCorr.R()+EPS))));
	outCorr.P(betai(0.5*outCorr.V(), 0.5, outCorr.V() / ( outCorr.V()+(outCorr.T()*outCorr.T()) )));
	if(one_tail) outCorr.P(outCorr.P()*0.50);
}

/* OriginCorrelate 
 * ---------------
 * Correlation forced through the origin
 */
void TBivSample::OriginCorrelate(TCorrelation& outCorr, bool one_tail /*= false*/)
{
	double psum(0);
	for(int i=0;i<_xv.N();i++) {
		psum += (_xv[i] * _yv[i]);
	}
	outCorr.R( psum / (sqrt ((_xv.SumValueSquares()*_yv.SumValueSquares()) ) ) );
	outCorr.V(_xv.N()-2);
	outCorr.Z(0.5*log((1.0+outCorr.R()+EPS)/(1.0-outCorr.R()+EPS)));
	outCorr.T(outCorr.R()*sqrt(outCorr.V()/((1.0-outCorr.R()+EPS)*(1.0+outCorr.R()+EPS))));
	outCorr.P(betai(0.5*outCorr.V(), 0.5, outCorr.V() / ( outCorr.V()+(outCorr.T()*outCorr.T()) )));
	if(one_tail) outCorr.P(outCorr.P()*0.50);
	outCorr.P(std::abs(outCorr.P()));
}

/* Monte Carlo (randomization) method for testing correlation significance 
 * (diff from zero)
 * -----------------------------------------------------------------------
 */
double TBivSample::OriginMCCorrelate(TCorrelation& outCorr, int nReps, bool one_tail /*= false*/)
{
	TSample outBoot;
	outBoot.Reserve(nReps);
	int ct(0), size(this->N()), i;
	double psum(0);
	double r;
	double expVal(0);
	for(int i=0;i<size;i++) {
		psum += (_xv[i] * _yv[i]);
	}
	r = psum / (sqrt ((_xv.SumValueSquares()*_yv.SumValueSquares()) ) );
	outCorr.R(r);
	outCorr.Z(0.5*log((1.0+outCorr.R()+EPS)/(1.0-outCorr.R()+EPS)));
	// now do reps
	TBivSample rep(*this);
	for(i=0;i<nReps;i++) {
		rep.Shuffle();
		psum=0;
		for(int j=0;j<size;j++) 	psum += (rep._xv[j] * rep._yv[j]);
		outBoot.Add(psum / (sqrt ((rep._xv.SumValueSquares() * (rep._yv.SumValueSquares())) ) ));
	}

	expVal = outBoot.Mean();

	// one tail:
	if(one_tail) {
		if(r > expVal) {
			for(i=0;i<nReps;i++) {
				if(outBoot._V[i]>r) ct++;
			}
		} else {
			for(i=0;i<nReps;i++) {
				if(outBoot._V[i]<r) ct++;
			}
		}
	} else {
		double diff = ((r>expVal)?r-expVal:expVal-r);
			for(i=0;i<nReps;i++) {
				if((outBoot._V[i]>(expVal+diff)) || (outBoot._V[i]<(expVal-diff))) ct++;
			}
	}
	outCorr.P(double(ct) / double(nReps));
	return expVal;
}


/* Bootstrap method for testing significance -- doesn't work
 * I need to check this as of 12/5/00 - Dylan
 */
void TBivSample::OriginBootstrapCorrelate(TCorrelation& outCorr, int nReps, bool one_tail /*= false*/)
{
	TSample outBoot;
	outBoot.Reserve(nReps);
	int ct(0), size(this->N()), rand, i, j;
	double psum(0);
	double r;
	double expVal(0);
	for(i=0;i<size;i++) {
		psum += (_xv[i] * _yv[i]);
	}
	r = psum / (sqrt ((_xv.SumValueSquares()*_yv.SumValueSquares()) ) );
	outCorr.R(r);
	outCorr.Z(0.5*log((1.0+outCorr.R()+EPS)/(1.0-outCorr.R()+EPS)));
	// now do reps
	TBivSample rep;
	rep.Reserve(size);
	for( i=0;i<nReps;i++) {
		for(j=0;j<size;j++) {
			rand=random_integer(0,size-1);
			rep.Add(this->_xv[rand],this->_yv[rand]);
		}
		psum=0;
		for( j=0;j<size;j++) 	psum += (rep._xv[j] * rep._yv[j]);
		outBoot.Add(psum / (sqrt ((rep._xv.SumValueSquares() * (rep._yv.SumValueSquares())) ) ));
	}

	expVal = outBoot.Mean();

	// one tail:
	if(one_tail) {
		if(expVal > 0) {
			for(i=0;i<nReps;i++) {
				if(outBoot._V[i]<0) ct++;
			}
		} else {
			for(i=0;i<nReps;i++) {
				if(outBoot._V[i]>0) ct++;
			}
		}
	} else {
			for(i=0;i<nReps;i++) {
				if(expVal>0 && (outBoot._V[i]<0 )) ct++;
				else if(expVal<0 && (outBoot._V[i]>0 )) ct++; 
			}
	}
	outCorr.P(double(ct) / double(nReps));
}


double TBivSample::SignTest(bool one_tail/*=false*/)
{
	int k(0), zeros(0);

	for (int i=0;i<this->N();i++) {
		if(_xv[i] == _yv[i]) zeros++;
		else if(_xv[i]>_yv[i]) k++;
	}

	int n = (this->N()-zeros);
	if(k<(n-k)) k = n-k;

	return (one_tail? betai(k,n-k+1, 0.5):2.0*betai(k,n-k+1, 0.5));
}

double TBivSample::SumCrossProducts()
{
	double sum(0);
	for(int i=0;i<_xv.N();i++) sum+= (_xv.Deviate(_xv[i]) * _yv.Deviate(_yv[i]));
	return sum;
}



/* TBivSample Input / output
 * -------------------------
 */

std::ostream& operator<<(std::ostream& theStream, TBivSample& theS)
{
	theStream << theS.N() << '\t';
	for(int i=0;i<theS.N();i++) {
		theStream << theS.XObs(i) << '\t' << theS.YObs(i) << '\n';
	}
	return theStream;
}


std::istream& operator>>(std::istream& theStream, TBivSample& theS)
{
	int n;
	double tx, ty;
	theStream >> n;
	for(int i=0;i<n;i++) {
		theStream >> tx >> ty;
		theS._xv.Add(tx);
		theS._yv.Add(ty);
	}

	return theStream;
}




/******************************************************************************
 * Class THistogram
 ******************************************************************************/

//ctor
THistogram::THistogram(const int n, const double xlow, const double xhigh)
: _xHigh(xhigh),
  _xLow(xlow),
  _N(n),
  _underflow(0),
  _overflow(0)
{}

//dtor
THistogram::~THistogram()
{}


bool THistogram::Add(double x) 
{
	if(x<_xLow) {
		_underflow++;
		return false;
	} else if(x>_xHigh) {
		_overflow++;
		return false;
	} else {
		(_counts[Bin(x)])++;
		return true;
	}
}


void THistogram::Fill(const TSample& S)
{
	_counts.reserve(this->_N);
	for(int i=0;i<S.N();i++) this->Add(S.Obs(i));
}


/***************************************************************************
 *  Statistics
 **************************************************************************/


// Ttest
// -----
// Tests two samples for significanlty different means.  Assumes equal variance
double TT_Test::Ttest(TSample& s1, TSample& s2)  
{ 
	using namespace DWS;

	double v1, v2, svar;
	v1 = s1.N()-1.0;
	v2 = s2.N()-1.0;
	this->_v = v1+v2;
	svar = (v1*s1.Variance() + v2*s2.Variance())/this->_v;

	this->_t = (s1.Mean() - s2.Mean()) / sqrt(svar * (1.0/s1.N() + 1.0/s2.N()) );
	this->_p = betai(0.5*this->V(), 0.5, this->V()/(this->V()+this->T()*this->T()) );
	return (this->P());
}


// Tutest
// -----
// Tests two samples for significanlty different means.  Does not assume equal variance
double TT_Test::Tutest(TSample& s1, TSample& s2)  
{ 
	using namespace DWS;

	double v1, v2;
	v1 = s1.N()-1.0;
	v2 = s2.N()-1.0;
	this->_v = DWS::sqr(s1.Variance()/s1.N() + s2.Variance()/s2.N()) / (DWS::sqr(s1.Variance()/s1.N())/v1 + DWS::sqr(s2.Variance()/s2.N())/v2);
	this->_t = (s1.Mean() - s2.Mean()) / std::sqrt(s1.Variance()/s1.N() + s2.Variance()/s2.N() );
	this->_p = betai(0.5*this->V(), 0.5, this->V()/(this->V()+DWS::sqr(this->T()) ) );
	return (this->P());
}
}; // namespace DWS
