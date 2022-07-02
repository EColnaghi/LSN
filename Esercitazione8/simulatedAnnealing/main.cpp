#include <iostream>
#include <fstream>
#include <algorithm>
#define _USE_MATH_DEFINES
#include <cmath>
#include <string>
#include <vector>
#include "random.h"
#include "main.h"
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


	int blockLength		 = 1000;
	int nBlocks	     = 100;
	double mu		 = 1.41405;
	double sigma	 = 0.047;
	double muJump	 = 0.02;
	double sigmaJump = 0.02;
	double sampleJump=0.5;
	double exponent=-1;
	int SAStep=100;
	double T_0=1;
	Params params(mu, sigma, &rnd, sampleJump, nBlocks, blockLength);
	std::vector<double> positions, integrals,errors;
	int nStep=100;
	std::cout<<"using The metropolis to sample from the wavefunction..."<<std::endl;
	std::cout<<"Using "<<blockLength<<" blocks, "<<nBlocks<<" steps for each block"<<std::endl; 
		
	
	
	std::ofstream output("positions.csv");
	output<<"x"<<std::endl;
	
	
	std::ofstream output2("params.csv");
	output2<<"mu,sigma,integral,error"<<std::endl;
	
	std::ofstream output3("integrals.csv");
	output3<<"nBlock,integral,error"<<std::endl;
	
	
	double T=2;
	std::vector<Params> paramsVector;
	std::cout<<"Inizio a temperatura: "<<T_0<<std::endl;
	for(int i = 0; i < SAStep;i++){	
		T=temp(T_0,i, exponent);
		std::cout<<"Iterazione: "<<i<<"/"<<SAStep<<"   Temperatura: "<<T<<" ";
		//create the boltzmann distribution in the parameters space we want to optimize and the jumping function for the metropolis algorithm
		auto toOptimize = [T](Params params)->double{ return std::exp(-params.getIntegral()/T);};
		auto jump       = [&rnd, muJump, sigmaJump, nBlocks, blockLength, sampleJump](Params p)->Params{ Params p2;
			p2.setRnd(&rnd);
			p2.setSampleJump(sampleJump);
			p2.setnBlocks(nBlocks);
			p2.setBlockLength(blockLength);
			p2.setParams(p.getMu()+(rnd.Rannyu(-muJump,muJump)),p.getSigma()+(rnd.Rannyu(-sigmaJump,sigmaJump)));
			return p2;};
		//perform nStep metropolis steps in the parameters space
		metropolisSampling(paramsVector,toOptimize, jump, params,nStep,rnd);
		std::cout<<"SA step: "<<i<<" temp: "<<T<<std::endl;
		std::cout<<"Integral: "<<paramsVector.back().getIntegral()<<"Error: "<<paramsVector.back().getError()<<"  mu: "<<paramsVector.back().getMu()<<"  sigma: "<<paramsVector.back().getSigma()<<std::endl<<std::endl;;
		params=paramsVector.back();
	}
	
	//save the parameters in an output file
	for(int i =0; i<(int)paramsVector.size();i++){output2<<paramsVector[i].getMu()<<","<<paramsVector[i].getSigma()<<","<<paramsVector[i].getIntegral()<<","<<paramsVector[i].getError()<<std::endl;}
	
	//save in an output file the positions sampled with the final parameters
	for(int j =0; j<blockLength;j++){
		positions=paramsVector.back().samplePositions();
		for(int i =0; i<(int)positions.size();i++){output<<positions[i]<<std::endl;}
	}

	//save and print in an output file the integrals with the final parameters
	paramsVector.back().printIntegral(integrals, errors);
	for(int i=0;i<nBlocks;i++){
		output3<<i<<","<<integrals[i]<<","<<errors[i]<<std::endl;
	}
	
	rnd.SaveSeed();
	return 0;
}


//return the temperature at the i-th step of the simulated annealing algorithm
double temp(double T_0,int i, double exponent){
	return T_0*(double)std::pow((double)i+1.0, exponent);
	//return T_0*(double)std::exp((double)i*exponent);
}

