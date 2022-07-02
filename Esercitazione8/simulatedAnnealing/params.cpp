#include <iostream>
#include <fstream>
#include <algorithm>
#define _USE_MATH_DEFINES
#include <cmath>
#include <string>
#include <vector>
#include <iomanip>
#include <list>
#include "params.h"
#include "random.h"
#include "metropolis.h"

Params::Params(double mu, double sigma, Random *rnd, double sampleJump, double N){
	this->mu=mu;
	this->sigma=sigma;
	this->nBlocks=N;
	this->rnd=rnd;
	this->sampleJump=sampleJump;	
	this->blockLength=1;
	this->computeIntegral();
}


Params::Params(double mu, double sigma, Random *rnd, double sampleJump, double nBlocks, double blockLength){
	this->mu=mu;
	this->sigma=sigma;
	this->nBlocks=nBlocks;
	this->blockLength=blockLength;
	this->rnd=rnd;
	this->sampleJump=sampleJump;
	this->computeIntegral();
}


Params::Params(Random *rnd, double sampleJump, double nBlocks, double blockLength){
	this->nBlocks=nBlocks;
	this->blockLength=blockLength;
	this->rnd=rnd;
	this->sampleJump=sampleJump;
	this->computeIntegral();
}

Params::Params(Random *rnd, double nBlocks, double blockLength){
	this->nBlocks=nBlocks;
	this->blockLength=blockLength;
	this->rnd=rnd;
	this->sampleJump=1;
	this->computeIntegral();
}

Params::Params(){
	this->mu=0;
	this->sigma=0;
	this->nBlocks=1000;
	this->integral=0;
	this->sampleJump=1;
	
	this->blockLength=1;
}

void Params::setRnd(Random *rnd){
	this->rnd=rnd;
}

void Params::setMu(double mu){
	this->mu=mu;
	this->computeIntegral();
}

void Params::setSigma(double sigma){
	this->sigma=sigma;
	this->computeIntegral();
}

void Params::setParams(double mu,double sigma){
	this->mu=mu;
	this->sigma=sigma;
	this->computeIntegral();
}

void Params::setSampleJump(double sampleJump){
	this->sampleJump=sampleJump;
}

void Params::setnBlocks(int nBlocks){
	this->nBlocks=nBlocks;
}

void Params::setBlockLength(int blockLength){
	this->blockLength=blockLength;
}

double Params::getMu(){
	return this->mu;
}

double Params::getSigma(){
	return this->sigma;
}


double Params::getIntegral(){
	return this->integral;
}

double Params::getError(){
	return this->error;
}

/*return a vector with a sample of positions of the wavefunction. The length of the vector is the length of a single block used to calculate the integral*/
std::vector<double> Params::samplePositions(){
	std::vector<double> positions;
	double x_0			=0;
	double mu2         	=this->mu;
	double sigma2      	=1/(2*std::pow(this->sigma,2));
	auto psi       	  	=[mu2, sigma2](double x)->double{return std::pow(std::exp(-std::pow((x-mu2),2)*sigma2)+std::exp(-std::pow((x+mu2),2)*sigma2),2);};
	auto jump         	=[this](double x)->double{	
		double random=this->sampleJump*(this->rnd->Rannyu()-0.5);
		double toReturn=x+random;	
		return toReturn;
		};
	metropolisSampling(positions, psi, jump, x_0,this->blockLength, *rnd);
    return positions;
}

/*Calculate the integral via Monte Carlo integration. The statistical uncertainties are obtained through datablocking. To do so are used nBlocks blocks
of length blockLength. The positions are sampled with the metropolis*/
void Params::computeIntegral(){
	std::vector<double> positions;
	double meanIntegral =0;
	double sqrIntegral  =0;
	double integral		=0;
	double mu2         	=this->mu;
	//double sigma2      	=1/(2*std::pow(this->sigma,2));
	double sigma2      	=std::pow(this->sigma,2);
	
	auto hamiltonian  	=[sigma2, mu2](double x){
		double potentialEnergy=(std::pow(x,4)-5*std::pow(x,2)/2);
		double exponential=std::exp(2*x*mu2/(sigma2));
		double kineticEnergy = -(x*x-((2*x*mu2*(-1+exponential))/(1+exponential))+mu2*mu2-sigma2)/(2*std::pow(sigma2,2));
		return kineticEnergy+potentialEnergy;
		};
	
	
	for(int i = 0; i<nBlocks; i++){
		meanIntegral=0;
		positions=this->samplePositions();
		for(int j = 0; j < (int)positions.size();j++){
			meanIntegral += hamiltonian(positions[j]);
		}
		meanIntegral/=(double) positions.size();
		integral+=meanIntegral;
		sqrIntegral+=meanIntegral*meanIntegral;
		
	}
	
	sqrIntegral   /=this->nBlocks;
	integral      /=this->nBlocks;
	this->integral =integral;
	this->error    =std::sqrt((sqrIntegral-std::pow(integral,2))/(this->blockLength-1));
}

/*save the progressive means of the integrals and the errors obtained through datablocking in the vectors integral and errors*/
void Params::printIntegral(std::vector<double>& integrals, std::vector<double>& errors){
	std::vector<double> positions;
	double meanIntegral =0;
	std::vector<double> sqrIntegral;
	double mu2         	=this->mu;
	//double sigma2      	=1/(2*std::pow(this->sigma,2));
	double sigma2      	=std::pow(this->sigma,2);
	
	auto hamiltonian  	=[sigma2, mu2](double x){
		double potentialEnergy=(std::pow(x,4)-5*std::pow(x,2)/2);
		double exponential=std::exp(2*x*mu2/(sigma2));
		double kineticEnergy = -(x*x-((2*x*mu2*(-1+exponential))/(1+exponential))+mu2*mu2-sigma2)/(2*std::pow(sigma2,2));
		return kineticEnergy+potentialEnergy;
		};
	
	
	for(int i = 0; i<nBlocks; i++){
		meanIntegral=0;
		positions=this->samplePositions();
		for(int i = 0; i < (int)positions.size();i++){
			meanIntegral += hamiltonian(positions[i]);
		}
		meanIntegral/=(double) positions.size();
		integrals.push_back(meanIntegral);
		sqrIntegral.push_back(meanIntegral*meanIntegral);
		
	}
	
	for(int i = 1; i<nBlocks; i++){
		integrals[i]+=integrals[i-1];
		sqrIntegral[i]+=sqrIntegral[i-1];	
	}
	
	errors.push_back(0);
	for(int i = 1; i<nBlocks; i++){
		integrals[i]/=(i+1);
		sqrIntegral[i]/=(i+1);
		errors.push_back(std::sqrt((sqrIntegral[i]-std::pow(integrals[i],2))/((double)i)));
	}
	
}


double error(std::vector<double> A, std::vector<double> A2, int n){
	if(n==0){
		return 0;
	}else{
		return std::sqrt((A2[n]-pow(A[n],2))/(double)n);
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