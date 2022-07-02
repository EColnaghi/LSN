

#include <iostream>
#include <fstream>
#include <algorithm>
#define _USE_MATH_DEFINES
#include <cmath>
#include <string>
#include <vector>
#include "random.h"

double error(std::vector<double> avg, std::vector<double> avg2, int n){
	if(n==0){
		return 0;
	}else{
		return sqrt((avg2[n]-pow(avg[n],2))/n);
	}
}

double func(double x){
	return M_PI*std::cos(M_PI*x/2)/2;
}



double p(double x){
	return 3*(1-std::pow(x,2))/2;
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


	
	/*Declaration of variables and inizializations of the output stream*/
	std::ofstream output("monteCarloIntegration.csv");
	output<<"uniformSampling,uniformError,importanceSampling,importanceError"<<std::endl;
	std::ofstream output2("acceptanceRate.csv");
	output2<<"nBlocks,acceptanceRate"<<std::endl;
	double x;
	double y;
	double integral=0;
	double integral2=0;
	std::vector<double> meanUniformIntegral;
	std::vector<double> meanISIntegral;
	std::vector<double> meanSquaredUniformIntegral;
	std::vector<double> meanSquaredISIntegral;	
	std::vector<double> msdUniformIntegral;	
	std::vector<double> msdISIntegral;
	int M=1e4;
	int N=1e2;
	double attempted=0;
	double accepted=0;
	
	/*Computate the integral*/
	for(int j=0; j<N; j++){	
		integral=0;
		integral2=0;
		for(int i=0; i<M;i++){
			x=rnd.Rannyu();
			y=rnd.Rannyu(0,1.5);
			attempted++;
			//use sampled the point to calculate the integral with uniform sampling
			integral+=func(x);
			//use reject method to sample from a distribution with pdf p(x)
			while(y>p(x)){
				x=rnd.Rannyu();
				y=rnd.Rannyu(0,1.5);
				attempted++;
			}
			accepted++;
			//use sampled point to calculate the integral with importance sampling
			integral2+=func(x)/p(x);
		}

		output2<<j<<','<<accepted/attempted<<std::endl;
		meanUniformIntegral.push_back(integral/M);
		meanSquaredUniformIntegral.push_back(pow(integral/M,2));
		
		meanISIntegral.push_back(integral2/M);
		meanSquaredISIntegral.push_back(pow(integral2/M,2));
	}
	
	for(int i=0; i<N; i++){
		
		meanUniformIntegral[i+1]+=meanUniformIntegral[i];
		meanSquaredUniformIntegral[i+1]+=meanSquaredUniformIntegral[i];
		
		meanISIntegral[i+1]+=meanISIntegral[i];
		meanSquaredISIntegral[i+1]+=meanSquaredISIntegral[i];
	}	
	
	for(int i=0; i<N; i++){	
		meanUniformIntegral[i]/=(i+1);
		meanSquaredUniformIntegral[i]/=(i+1);
		
		meanISIntegral[i]/=(i+1);
		meanSquaredISIntegral[i]/=(i+1);
		
		
		msdUniformIntegral.push_back(error(meanUniformIntegral,meanSquaredUniformIntegral,i));
		msdISIntegral.push_back(error(meanISIntegral,meanSquaredISIntegral,i));
		output<<meanUniformIntegral[i]<<","<<msdUniformIntegral[i]<<","<<meanISIntegral[i]<<","<<msdISIntegral[i]<<std::endl;

	}		
	
		
	
	rnd.SaveSeed();
	return 0;
}
