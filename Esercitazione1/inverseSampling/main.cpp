/*In this program we want to generate data to verify the generalized central limit theorem
for three different distributions: an exponential distribution,*/

#include <iostream>
#include <fstream>
#include <algorithm>
#define _USE_MATH_DEFINES
#include <cmath>
#include <string>
#include <vector>
#include "random.h"


//return the mean of the sum of n exponential random variables extracted from the distribution e^(-lambda*x)
double exponential(Random& rnd, double lambda, int n){
	double sum=0;
	
	for(int i = 0; i<n; i++){
		sum+=rnd.exponential(lambda);
	}
	return sum/n;
}


//return the mean of the sum of n random variables extracted from the pdf of a standard, six sided dice
double standard(Random& rnd, double mu, double gamma, int n){
	double sum=0;
	
	for(int i = 0; i<n; i++){
		sum+=std::floor(rnd.Rannyu(1,7));
	}
	return sum/n;
}


//return the mean of the sum of n random variables extracted from a cauchy distribution centered in mu with FWHM gamma
double cauchy(Random& rnd, double mu, double gamma, int n){
	double sum=0;
	
	for(int i = 0; i<n; i++){
		sum+=rnd.cauchy(mu, gamma);
	}
	return sum/n;
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
	//declare the variables
	int n=10000;			//size of the sampling
	int N[]={1,2,10,100};	//size of the sample of which we will take the mean
	double lambda=1;		//constant of the exponential distribution
	double mu=0;			//center of the cauchy distribution
	double gamma=1;			//gamma parameter of the cacuhy distribution
	
	std::ofstream output("distributions.csv");
	output<<"e1,c1,s1,e2,c2,s2,e10,c10,s10,e100,c100,s100,"<<std::endl;
	
	for(int i=0; i<n;i++){
		for(int j:N){
			output<<exponential(rnd,lambda,j)<<',';
			output<<cauchy(rnd,mu, gamma,j)<<',';
			output<<standard(rnd,mu, gamma,j)<<',';
			
		}
		output<<std::endl;
	}
	

	
	
	rnd.SaveSeed();
	return 0;
}

