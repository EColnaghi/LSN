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



double exponential(Random& rnd, double lambda, int n){
	double sum=0;
	
	for(int i = 0; i<n; i++){
		sum+=rnd.exponential(lambda);
	}
	return sum/n;
}

double standard(Random& rnd, double mu, double gamma, int n){
	double sum=0;
	
	for(int i = 0; i<n; i++){
		sum=+std::floor(rnd.Rannyu()*6)+1;
	}
	return sum/n;
}

double cauchy(Random& rnd, double mu, double gamma, int n){
	double sum=0;
	
	for(int i = 0; i<n; i++){
		sum=+rnd.cauchy(mu, gamma);
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

	int M=4;
	int n=10000;
	int N[]={1,2,10,100};
	double lambda=1;
	double mu=0;
	double gamma=1;
	
	std::ofstream output("distributions.csv");
	output<<"e1,c1,u1,e2,c2,u2,e10,c10,u10,e100,c100,u100,"<<std::endl;
	
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

