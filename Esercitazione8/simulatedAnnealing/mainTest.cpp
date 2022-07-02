#include <iostream>
#include <fstream>
#include <algorithm>
#define _USE_MATH_DEFINES
#include <cmath>
#include <string>
#include <vector>
#include "random.h"
#include "mainTest.h"
#include "params.h"
#include "metropolis.h"

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


	
	std::vector<double> integrals;
	int blockLength		 = 1000;
	int nBlocks	     = 1000;
	double mu		 = 1.65917;
	double sigma	 = 0.052317;
	double muJump	 = 0.05;
	double sigmaJump = 0.05;
	double sampleJump=0.5;
	double exponent=-0.5;
	int SAStep=25;
	double T_0=10;
	Params params(mu, sigma, &rnd, sampleJump, nBlocks, blockLength);
	
	std::cout<<"Integrale con mu: "<<params.getMu()<<" e sigma: "<<params.getSigma()<<" : "<<params.getIntegral()<<std::endl;
	/*int nStep=100;
	std::cout<<"using The metropolis to sample from the wavefunction..."<<std::endl;
	std::cout<<"Using "<<blockLength<<" blocks, "<<nBlocks<<" steps for each block"<<std::endl;
		
	
	
	std::ofstream output("positions.csv");
	output<<"x"<<std::endl;
	
	
	std::ofstream output2("params.csv");
	output2<<"mu,sigma,I,errs"<<std::endl;
	
	
	double T=2;
	std::vector<Params> paramsVector;
	std::cout<<"Inizio a temperatura: "<<T_0<<std::endl;
	for(int i = 0; i < SAStep;i++){	
		T=temp(T_0,i, exponent);
		std::cout<<"Iterazione: "<<i<<"/"<<SAStep<<"   Temperatura: "<<T<<" ";
		auto toOptimize = [T](Params params)->double{ return std::exp(-params.getIntegral()/T);};
		auto jump       = [&rnd, muJump, sigmaJump, nBlocks, blockLength, sampleJump](Params p)->Params{ Params p2;
			p2.setRnd(&rnd);
			p2.setSampleJump(sampleJump);
			p2.setnBlocks(nBlocks);
			p2.setBlockLength(blockLength);
			p2.setParams(p.getMu()+muJump*(rnd.Rannyu()-0.5),p.getSigma()+sigmaJump*(rnd.Rannyu()-0.5));
			return p2;};
		std::cout<<"entro in verbose"<<std::endl;
		metropolisSamplingVerbose(paramsVector,toOptimize, jump, params,nStep,rnd);
		std::cout<<"Integral: "<<paramsVector.back().getIntegral()<<"Error: "<<paramsVector.back().getError()<<"  mu: "<<paramsVector.back().getMu()<<"  sigma: "<<paramsVector.back().getSigma()<<std::endl<<std::endl;;
		params=paramsVector.back();
	}
	
	for(int i =0; i<(int)paramsVector.size();i++){output2<<paramsVector[i].getMu()<<","<<paramsVector[i].getSigma()<<","<<paramsVector[i].getIntegral()<<","<<paramsVector[i].getError()<<std::endl;}
	/*for(int j =0; j<blockLength;j++){
		std::vector<double> positions=paramsVector.back().samplePositions();
		for(int i =0; i<(int)positions.size();i++){output<<positions[i]<<std::endl;}
	}*/
	//std::cout<<"mu, sigma: "<<paramsVector.back().getMu()<<","<<paramsVector.back().getSigma()<<std::endl;
	rnd.SaveSeed();
	return 0;
}


//return the temperature at the i-th step of the simulated annealing algorithm
double temp(double T_0,int i, double exponent){
	return T_0*(double)std::pow((double)i+1.0, exponent);
}

