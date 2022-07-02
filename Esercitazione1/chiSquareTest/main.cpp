/*In this file we generate a set of  values of the chi squared statistics for 
the given random number generator, which generates data from a uniform pdf. 
We will generate N values of the chi squared statistics, and we will store it in a csv file*/

#include <iostream>
#include <fstream>
#include <cmath>
#include <string>
#include <vector>
#include "random.h"


double error(std::vector<double> avg, std::vector<double> avg2, int n){
	if(n==0){
		return 0;
	}else{
		return sqrt((avg2[n]-std::pow(avg[n],2))/n);
	}
}

int main (int argc, char *argv[]){

   Random rnd;
   int seed[4];
   int p1, p2;
   std::ifstream Primes("Primes");
   if (Primes.is_open()){
      Primes >> p1 >> p2 ;
   } else std::cerr << "PROBLEM: Unable to open Primes" << std::endl;
   Primes.close();

   std::ifstream input("seed.in");
   std::string property;
   if (input.is_open()){
      while ( !input.eof() ){
			input >> property;
			if( property == "RANDOMSEED" ){
				input >> seed[0] >> seed[1] >> seed[2] >> seed[3];
				rnd.SetRandom(seed,p1,p2);
			}
		}
		input.close();
	} else std::cerr << "PROBLEM: Unable to open seed.in" << std::endl;
	/*definition of variable:
	-N is the number of independent experiment we want to perform
	-M is the number of bin we will divide the data into
	-n is the number of numbers thegerated for every experiment
	-observed is a temporary histogram in which we will store the number of observation in every bin
	-expected is the number of expected observations in every bin
	*/
	int M=100;
	int observed [M]={};
	int n=1000;
	int N=100;
	double expected=n/M;
	double next=0;
	double chi=0;
	std::vector<double> chiStatistics;
	
	for(int i=0; i<N;i++){
		chi=0;
		std::fill(observed,observed+M,0);
		//Generate n numbers, and fill the histogram 'observed'
		for(int j=0; j<n; j++){
			next=rnd.Rannyu();
			observed[(int)std::floor(next*M)+1]++;
		}
		//Evaluate the chi statistics
		for(int j=0;j<M;j++){
			chi+=std::pow((expected-observed[j]),2)/expected;
		}
		chiStatistics.push_back(chi);
	}
	
	//Write the N values in a csv file
	std::ofstream output("chiSquaredStatistics.csv");
	output<<"chiSquaredStatistics"<<std::endl;
	for(int i=0;i<N;i++){
		output<<chiStatistics[i]<<std::endl;
	}
	
	
	rnd.SaveSeed();
	return 0;
}
