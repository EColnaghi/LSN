

#include <iostream>
#include <fstream>
#include <algorithm>
#define _USE_MATH_DEFINES
#include <cmath>
#include <string>
#include <vector>
#include "random.h"

double error(std::vector<double> A, std::vector<double> A2, int n){
	if(n==0){
		return 0;
	}else{
		return std::sqrt((A2[n]-pow(A[n],2))/n);
	}
}

void progError(std::vector<double> &A, std::vector<double> &A2, std::vector<double> &err){
	int N= A.size();
	for(int i = 1; i<N;i++){
		A[i]+=A[i-1];
		A[i-1]/=i;
		
		A2[i]+=A2[i-1];
		A2[i-1]/=i;
		
	}
	
	A[N-1]=A[N-1]/N;
	
	A2[N-1]=A2[N-1]/N;
	
	
	for(int i = 0; i<N;i++){
		err.push_back(error(A,A2, i));	
	}
	
}

double N(double x){
	return 0.5 * (1 + std::erf(x / std::sqrt(2)));
}


/*double d1(double sigma, double T, double s0, double K, double r){
	return 1/(sigma * std::sqrt(T)) * (std::log(s0 / K) + (r + (std::pow(sigma,2)) / 2) * T);
}*/

/*double d2(double d1, double sigma, double T){
	return d1 - sigma * std::sqrt(T);
}*/

/*double call(double s0, double K, double r, double T,double sigma){
	double d1 = 1/(sigma * std::sqrt(T)) * (std::log(s0 / K) + (r + (std::pow(sigma,2)) / 2) * T);
	double d2 = d1 - sigma * std::sqrt(T);
	return s0 * N(d1) - K * std::exp(-r * T) * N(d2);
}*/


/*double put(double s0, double K, double r, double T,double sigma){
	double d1 = 1/(sigma * std::sqrt(T)) * (std::log(s0 / K) + (r + (std::pow(sigma,2)) / 2) * T);
	double d2 = d1 - sigma * std::sqrt(T);
	return 0;
}*/

double incrementGBM(Random &rnd, double t, double mu, double sigma){
	return std::exp((mu-0.5*std::pow(sigma,2))*t+sigma*rnd.Gauss(0,1)*std::sqrt(t));
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
	std::ofstream output("blackScholesPricing.csv");
	//output<<"nBlocks, directCall, discreteCall, directPut,discretizedPut,errDirectCall, errDiscretizedCall, errDirectPut,errDiscretizedPut"<<std::endl;
	output<<"nBlocks,directCall,discreteCall,directCallErr,discreteCallErr,directPut,discretePut,directPutErr,discretePutErr"<<std::endl;

	//declare and initialize the variables
	double s0=100;
	double deliveryTime=1;
	double strikePrice=100;
	double interest=0.1;
	double volatility=0.25;
	int nStep=100;
	std::vector<double> directCall;
	std::vector<double> discreteCall;
	std::vector<double> directCallSqr;
	std::vector<double> discreteCallSqr;
	std::vector<double> directCallErr;
	std::vector<double> discreteCallErr;
	
	std::vector<double> directPut;
	std::vector<double> discretePut;
	std::vector<double> directPutSqr;
	std::vector<double> discretePutSqr;
	std::vector<double> directPutErr;
	std::vector<double> discretePutErr;
	
	
	
	int M=1e2;//number of 
	int N=1e3;//number of 
	double callAccDirect=0;
	double callAccDiscrete=0;
	
	double putAccDirect=0;
	double putAccDiscrete=0;
	
	
	double temp=0;
	double path=0;
	
	/*Computation of the prices*/
	for(int k=0; k<N;k++){
		callAccDirect=0;
		callAccDiscrete=0;

		putAccDirect=0;
		putAccDiscrete=0;
	/*In the next for the underlying geometric brownian motion is computed, both by taking a sile step of length 1 and
	by dividing my time in nStep time steps*/
		for(int i = 0; i<M; i++){
			temp=s0*incrementGBM(rnd, deliveryTime, interest, volatility);
			callAccDirect+=std::exp(-interest*deliveryTime)*std::max(temp-strikePrice, 0.0);
			putAccDirect+=std::exp(-interest*deliveryTime)*std::max(strikePrice-temp, 0.0);


			path=s0;		
			for(int j=0; j<nStep;j++){
				path*=incrementGBM(rnd, deliveryTime/nStep, interest, volatility);

			}
			callAccDiscrete+=std::exp(-interest*deliveryTime)*std::max(path-strikePrice, 0.0);
			putAccDiscrete+=std::exp(-interest*deliveryTime)*std::max(strikePrice-path, 0.0);
		}
		
	/*we add to the vector containing the prices the results. We also store the square of the average of M tries, which will
	later be used to compute statistical uncertainties*/
	directCall.push_back(callAccDirect/M);
	discreteCall.push_back(callAccDiscrete/M);
	
	directCallSqr.push_back(std::pow(callAccDirect/M,2));
	discreteCallSqr.push_back(std::pow(callAccDiscrete/M,2));

	directPut.push_back(putAccDirect/M);
	discretePut.push_back(putAccDiscrete/M);
	
	directPutSqr.push_back(std::pow(putAccDirect/M,2));
	discretePutSqr.push_back(std::pow(putAccDiscrete/M,2));			
	}
	
	/*In the next block I average the averages obatining on M results of the experiment and their squares to compute the
	uncertainties*/
	
	progError(directCall, directCallSqr, directCallErr);
	progError(discreteCall, discreteCallSqr, discreteCallErr);
	progError(directPut, directPutSqr, directPutErr);
	progError(discretePut, discretePutSqr, discretePutErr);
	/*Writing the results in the output file*/
	for(int i =0; i<N;i++){
		output<<i<<','<<directCall[i]<<','<<discreteCall[i]<<','<<directCallErr[i]<<','<<discreteCallErr[i]<<','<<directPut[i]<<','<<discretePut[i]<<','<<directPutErr[i]<<','<<discretePutErr[i]<<std::endl;
	}
	rnd.SaveSeed();
	return 0;
}
