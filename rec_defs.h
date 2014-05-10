/* rec_defs.h
 * D. Schwilk 2000
 * Definitions and helper functions for the Recombination project
 */


#include <algorithm>
#include <iostream>

#define GRID_SIZE 150
#define N_SEEDS 4
#define N_ALLELES 25
#define MORTALITY_RATE 0.1
#define _R_ 0.5
#define _r_ 0.05
#define INITIAL_TORCH_PROP 0.00
#define MUTATION_RATE 0.0001

#define BURN_NEIGHBORHOOD 1
#define SEED_NEIGHBORHOOD 1
#define POLLEN_NEIGHBORHOOD 1
#define ENV_CHANGE_PERIOD 8
#define LIGHTNING	0.01
#define MUTATE_FLIP 1  // determines environment mutation type, bit flip or random
#define N_ENV_CHANGE_LOCI 1 // number of loci to change each period

namespace DWS {


template<class T> struct print : public std::unary_function<T, void>
{
  print(std::ostream& out) : os(out), count(0) {}
  void operator() (T x) { os << x << ' '; ++count; }
  std::ostream& os;
  int count;
};


template <class T>
inline void swap(T &a, T &b)
{
	T c = a;
	a=b;
	b=c;
}

}; // namespace