/* File: dwsstats.h
 * ---------------
 * Dylan Schwilk 1998-2002
 * Stanford University
 *
 *  The DWS statistics module.
 *
 *  Classes in this file:
 * 
 *  1) TSample 
 *  2) and TBivSample classes are for general statistical
 *     analysis
 *
 *  3) THistogram class
 */

#ifndef dwstats_h
#define dwstats_h

#include <cmath>
#include <functional>
#include <vector>
#include <iostream>

#include "dws_math.h"

namespace DWS {

// repeated sampling types
enum RSamplingT {eWithReplace, eWithoutReplace};

typedef std::vector<double> sample_vector;
typedef std::pair<double, double> BivariateT;
typedef std::vector<BivariateT> BivariateVectorT;

// forward declarations
class TSample;


/* Statistic classes
 ********************************/

class TStatistic
{
public:
	TStatistic() : _p(1){};
	inline double P() {return _p;};
	inline double P(double theP) {_p=theP;return _p;};
protected:
	double _p;
};

class TSignTest
{
public:
	TSignTest(){};
};

class TT_Test : public TStatistic
{
public :

	TT_Test() : TStatistic(), _v(0)	{};
	double V() { return _v; };
	double T() {return _t; };
	// ttest assuming equal variances
	double Ttest(TSample& s1, TSample& s2);
	// ttest not assuming equal variances
	double Tutest(TSample& s1, TSample& s2);  
private:
	double _v;  // degrees of freedom
	double _t;
};

 
class TCorrelation : public TStatistic
{
public :
	TCorrelation() : TStatistic(), _r(0), _v(0), _z(0)
	{ }
	double R(){return _r;};
	double V(){return _v;};
	double Z() {return _z;};
	double T() {return _t;};

	double R(double r){_r=r;return _r;};
	double V(int v){_v=v;return _v;};
	double Z(double z) {_z=z;return _z;};
	double T(double t) {_t=t;return _t;}
private:
	double _r, _z, _t;
	double _v;
};


/******************************************************************************
 * Class TSample
 ******************************************************************************/
class TSample 
{
public:
	//constructors
	TSample();
	TSample(int size);
	TSample(const sample_vector& inV);
	//
	TSample(const TSample& copy);
	TSample& operator=(const TSample& source);
	~TSample();
	//
	// Methods
	inline void Erase(){_V.erase(_V.begin(), _V.end());Reset();}
	inline void Reserve(int n){	_V.reserve(n);}
	void Add(double x);
	inline double Min() const {return this->_range._min;}
	inline double Max() const {return this->_range._max;}
	int N() const { return this->_V.size();}
 	double Obs(int i) const {return _V[i];}
	double Mean();
	double Median();
	double Variance();
	inline double SD() {return(this->Variance()<=0.0?0:sqrt(Variance()));}
	double SumOfSquares();     //
	double SumValueSquares();
	inline double Deviate(double x);  // returns (x-xbar)
	double& operator[](int i) {return _V[i];};
	void Sort();
	void Shuffle();

	friend class TBivSample;
	friend std::ostream& operator<<(std::ostream& theStream, TSample& theS);
	friend std::istream& operator>>(std::istream& theStream, TSample& theS);

protected:
	void Reset();
private:

	sample_vector _V;
	bool _hasMean;
	bool _hasMedian;
	bool _hasVariance;
	bool _hasSOS;

	double _mean;
	double _median;
	double _variance;
	double _SOS;
	TRange<double> _range;

}; // class TSample


/***************************************************************************************
 * Class TBivSample                                                                    *
 * ------------------                                                                  *
 ***************************************************************************************/

class TBivSample {
public:
	TBivSample();
	TBivSample(sample_vector v1, sample_vector v2);
	~TBivSample();
	void Erase();
	void Reserve(int n);
	int N() const {return _xv.N();}
	TSample& X(); // returns sample of x vals
	TSample& Y();  
	//
	friend std::ostream& operator<<(std::ostream& theStream, const TBivSample& theS);
	friend std::istream& operator>>(std::istream& theStream, TBivSample& theS);
	void AddObservation(double x, double y);
	inline void Add(double x, double y) {AddObservation(x,y);}
	double XObs(int i) const {return this->_xv.Obs(i);}
	double YObs(int i) const {return this->_yv.Obs(i);}
	void Shuffle();

	/* Method: Correlate()
	 * ------------------
	 * produces a TCorrelation object (statistic) which contains the information about
	 * the correlation performed.  r, degrees of freedom(v), the z-value, and the p-value
	 * obtained by simply using the T-Test.
	 */
	void Correlate(TCorrelation& outCorr, bool one_tail = false);
	void OriginCorrelate(TCorrelation& outCorr, bool one_tail = false);

	/* Method: OriginMCCorrelate(nReps)
	 * ---------------------------------
	 * This method gives the 2-tailed p-value for a correlation using randomization
	 * techniques.  The value returned is the expected coefficient of correspondance
	 */
	double OriginMCCorrelate(TCorrelation& outCorr, int nReps, bool one_tail = false);
	/* */
	void OriginBootstrapCorrelate(TCorrelation& outCorr, int nReps, bool one_tail = false);

	/* Method: SignTest
	 * returns value for a sign test of a bivariate sample
	 */
	double SignTest(bool one_tail=false);
	double SumCrossProducts();
	
private:
	TSample _xv;
	TSample _yv;  // 
};


/***************************************************************************************
 * Class THistogram                                                                    *
 * ------------------                                                                  *
 ***************************************************************************************/

class THistogram {
public:
	THistogram(const int n, const double xlow, const double xhigh);
	~THistogram();
	void Clear();

	inline int NBins() const {return _N;}
	inline int Bin(double x) const {return int(((x-_xLow)/Width()));}
	inline double Width() const {return (_xHigh-_xLow)/((double)_N);}
	inline int Contents(int i) const {return _counts[i];}
	inline int Contents(double x) const {return Contents(Bin(x));}
	bool Add(double x);
	void Fill(const TSample& S);

	//
	friend std::ostream& operator<<(std::ostream& theStream, const THistogram& theS);
	friend std::istream& operator>>(std::istream& theStream, THistogram& theS);

	
private:
	int _N; // number of bins
	double _xLow;
	double _xHigh;


	std::vector<int> _counts;  // zero indexed array of bin counts (first bin is _counts[0])
	int _overflow;
	int _underflow;
};





/*********************************************************************************/
// helper function objects

class add_squared_diff : public std::binary_function<double, double, double>
{
public:
	double xbar;
	add_squared_diff(double m)
	{
		xbar=m;
	}
	double operator()(double y, double x) {
		return (y + (xbar-x)*(xbar-x));
	}
};


class add_square : public std::binary_function<double, double, double>
{
public:
	add_square(){}
	double operator()(double y, double x) {
		return (y + x*x);
	}
};





/*******************************************
 *  Functions declarations
 *******************************************/

double TTest(TSample& s1, TSample& s2, TT_Test& result);
std::ostream& operator<<(std::ostream& theStream, const TBivSample& theS);
std::istream& operator>>(std::istream& theStream, TBivSample& theS);


}; // namespace DWS


#endif  // TSample_h
