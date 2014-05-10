// Without stability criterion
#include "TPlant.h"
#include "dws_random.h"
#include "TLandscape.h"
#include <fstream>
#include <ctime>
#include <valarray>

//#define N 10000
#define NREPS 1
#define GRID 150
#define TORCH_START 0.005
#define CYCLES 1

double average(double *a , int n);
void WriteSample(const DWS::TBivSample& Sam);

int main(int argc, char **argv)
{
    DWS::set_randomizer_seed(-100);
	
    time_t   start;
	
	time( &start );

	std::ofstream ofile("out.txt");
	std::ofstream flamgrids("flamgrids.txt");
	std::ofstream fitgrids("fitgrids.txt");

	DWS::TSample TFitS;
	DWS::TSample DFitS;

	int changePeriod=8;
	double mortRate=0.1;
	int pneigh = 3;
	int sneigh = 3;
	int bneigh = 3;
	double mu = 0.00001;
	double r = 0.08;


	for(int reps=0;reps < NREPS; reps++) {

					double fit, tprop=0;
					DWS::TSample AgeSample;
					TLandscape LS(GRID);
					LS.SetChangePeriod(changePeriod);
					LS.SetMortalityRate(mortRate);
					LS.SetSeedNeighborhood(sneigh);
					LS.SetPollenNeighborhood(pneigh);
					LS.SetBurnNeighborhood(bneigh);
					LS.Setr(r);
					LS.SetMutationRate(mu);

				
					std::cout << "new landscape" << std::endl;	

					int count=0;
					while(true) {
						// without torches
						LS.RunForCycles(CYCLES);
						fit = LS.MeanFitness();
						LS.TorchDampFitness(TFitS, DFitS);
						std::cout << reps << '\t' << LS.Age() << '\t' << pneigh << '\t' << changePeriod << '\t' << mortRate << '\t' << LS.MeanFitness() << '\t' << fit << '\t' << LS.TorchProp() << std::endl;
						ofile << LS.Age() << '\t' << tprop << '\t' << fit << '\t' << TFitS.Mean() << '\t' << DFitS.Mean() << '\n';
						if (LS.Age() > 500) break;
					}

					// now add torches
					std::cout << "Adding torches" << std::endl;
					LS.AddTorches(TORCH_START);
					LS.WriteFitnessGrid(fitgrids);
					LS.WriteFlamGrid(flamgrids);
	
					while(true) {
						// with torches
						LS.RunForCycles(CYCLES);
						fit = LS.MeanFitness();
						tprop= LS.TorchProp();
						LS.TorchDampFitness(TFitS, DFitS);
						std::cout << reps << '\t' << LS.Age() << '\t' << pneigh << '\t' << changePeriod << '\t' << mortRate << '\t' << LS.MeanFitness() << '\t' << fit << '\t' << tprop << std::endl;
						ofile << LS.Age() << '\t' << tprop << '\t' << fit << '\t' << TFitS.Mean() << '\t' << DFitS.Mean() << '\n';

					
				// 3d grids -------------------------------------
						if(LS.Age()==511) {
							LS.WriteFitnessGrid(fitgrids);
							LS.WriteFlamGrid(flamgrids);
						}

						if(LS.Age()==527) {
							LS.WriteFitnessGrid(fitgrids);
							LS.WriteFlamGrid(flamgrids);
						}

						if(LS.Age()==543) {
							LS.WriteFitnessGrid(fitgrids);
							LS.WriteFlamGrid(flamgrids);
						}
	
						if(LS.Age()==735) {
							LS.WriteFitnessGrid(fitgrids);
							LS.WriteFlamGrid(flamgrids);
						}

						if(LS.Age()==903) {
							LS.WriteFitnessGrid(fitgrids);
							LS.WriteFlamGrid(flamgrids);
						}	
						if(LS.Age()==1503) {
							LS.WriteFitnessGrid(fitgrids);
							LS.WriteFlamGrid(flamgrids);
						}	
						if(LS.Age()==4007) {
							LS.WriteFitnessGrid(fitgrids);
							LS.WriteFlamGrid(flamgrids);
						}
						if(LS.Age()==10007) {
							LS.WriteFitnessGrid(fitgrids);
							LS.WriteFlamGrid(flamgrids);
						}
		
						if(tprop == 1.00) break;
					
					}
	}
	return 0;
}



double average(double *a , int n){
	double r=0;
	for(int i=0;i<n;i++) r+=a[i];
	return(r/double(n));
}


