

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



//Questa funzione prende in input un oggetto random e un numero di step e restituisce un vettore contenendo le deviazioni standard a ogni
//tempo per un random walk che si muove su un reticolo cubico tridimensionale di passo reticolare 1
std::vector<double> discreteRealization(Random* rnd, int nStep){
	double x[]={0,0,0};
	double verse=0;
	int direction=0;
	double msd=0;
	std::vector<double> msds;
	
	for(int i=0;i<nStep;i++){
		direction=(int)std::floor(rnd->Rannyu(0,3));
		verse=rnd->Rannyu(-1,1);
		verse=verse/std::abs(verse);
		x[direction]+=verse;
		msd=std::pow(x[0],2)+std::pow(x[1],2)+std::pow(x[2],2);
		msds.push_back(msd);
	}
	return msds;
}


//Questa funzione prende in input un oggetto random e un numero di step e restituisce un vettore contenendo le deviazioni standard a ogni
//tempo per un random walk che si muove a salti di lunghezza unitaria in una direzione scelta uniformemente su S2
std::vector<double> continueRealization(Random* rnd, int nStep){
	double x[]={0,0,0};
	double phi=0;
	double theta=0;
	double rho=0;
	double msd=0;
	std::vector<double> msds;
	
	for(int i=0;i<nStep;i++){
		//rho=rnd->Rannyu(-1,1);
		rho=1;
		phi=rnd->Rannyu(0,M_PI);
		theta=rnd->Rannyu(0,2*M_PI);
		x[0]+=rho*std::cos(theta)*std::sin(phi);
		x[1]+=rho*std::cos(theta)*std::cos(phi);
		x[2]+=rho*std::sin(theta);
		msd=std::pow(x[0],2)+std::pow(x[1],2)+std::pow(x[2],2);
		msds.push_back(msd);
	}
	return msds;
}

void sumVectors(std::vector<double> &a, std::vector<double> &b){
	for(int i=0;i<(int)a.size();i++){
		a[i]+=b[i];
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


	

	/*N e M contengono il numero di blocchi e il numero di realizzazioni per blocco, nStep la lunghezza delle realizzazioni,
	singleRealization conterrà la singola realizzazione di un random walk e i restanti vettori contengono le medie e le deviazioni
	standard mediate sui blocchi*/
	int M=1e3;
	int N=1e2;
	int nStep=1e2;
	std::vector<double> singleRealization;	
	std::vector<double> msdsDiscreteMean(nStep, 0);
	std::vector<double> blockDiscreteMean(nStep, 0);
	std::vector<double> blockDiscreteSquare(nStep, 0);
	std::vector<double> msdsDiscreteSquare(nStep, 0);
	
	
	
	std::vector<double> msdsContinueMean(nStep, 0);
	std::vector<double> blockContinueMean(nStep, 0);
	std::vector<double> blockContinueSquare(nStep, 0);
	std::vector<double> msdsContinueSquare(nStep, 0);

	/*Dichiarazione delle variabili e inizializzazione dell'output*/
	std::ofstream output("randomWalkerMSD.csv");
	output<<"t,discreteMSD,discreteError,continueMSD,continueError"<<std::endl;
	
	for(int i=0;i<N;i++){
		std::fill(blockDiscreteMean.begin(), blockDiscreteMean.end(), 0);
		std::fill(blockDiscreteSquare.begin(), blockDiscreteSquare.end(), 0);
		
		std::fill(blockContinueMean.begin(), blockContinueMean.end(), 0);
		std::fill(blockContinueSquare.begin(), blockContinueSquare.end(), 0);
		
		
		//Creo M realizzazioni e sommo gli spostamenti quadratici medi nei due vettori blockDiscreteMean e blockContinueMean
		//la somma degli spostamenti quadratici al tempo t e' contenuta nella posizione t del vettore
		for(int j=0;j<M;j++){
			singleRealization=discreteRealization(&rnd,nStep);
			std::transform (blockDiscreteMean.begin(), blockDiscreteMean.end(), singleRealization.begin(), blockDiscreteMean.begin(), std::plus<double>());
			
			
			singleRealization=continueRealization(&rnd,nStep);
			
			std::transform (blockContinueMean.begin(), blockContinueMean.end(), singleRealization.begin(), blockContinueMean.begin(), std::plus<double>());
			
		}
		
		//Divido per M, in modo da ottenere la media, e prendo i quadrati in modo da avere la media al quadrato, che che usero'
		// per calcolare l'errore
		for(int j=0;j<nStep;j++){
			blockDiscreteMean[j]/=(double)M;
			blockDiscreteSquare[j]+=std::pow(blockDiscreteMean[j],2);
			
			
			blockContinueMean[j]/=(double)M;
			blockContinueSquare[j]+=std::pow(blockContinueMean[j],2);
		}
		
		std::transform (msdsDiscreteMean.begin(), msdsDiscreteMean.end(), blockDiscreteMean.begin(), msdsDiscreteMean.begin(), std::plus<double>());
		std::transform (msdsDiscreteSquare.begin(), msdsDiscreteSquare.end(), blockDiscreteSquare.begin(), msdsDiscreteSquare.begin(), std::plus<double>());
	
		std::transform (msdsContinueMean.begin(), msdsContinueMean.end(), blockContinueMean.begin(), msdsContinueMean.begin(), std::plus<double>());
		std::transform (msdsContinueSquare.begin(), msdsContinueSquare.end(), blockContinueSquare.begin(), msdsContinueSquare.begin(), std::plus<double>());
		
		std::cout<<"Blocco numero: "<<i<<"/"<<N<<std::endl<<std::endl;
	}
	
	for(int i=0;i<nStep;i++){
		msdsDiscreteMean[i]/=(double)N;
		msdsDiscreteSquare[i]/=(double)N;	
		
		
		msdsContinueMean[i]/=(double)N;
		msdsContinueSquare[i]/=(double)N;	
		
		
	}
	
	for(int i=0; i<nStep; i++){	
		//output<<i<<','<<msdsDiscreteMean[i]<<','<< (std::pow(msdsDiscreteMean[i],2)-msdsDiscreteSquare[i])<<std::endl;
		output<<i<<','<<msdsDiscreteMean[i]<<','<<error(msdsDiscreteMean,msdsDiscreteSquare,i)<<','<<msdsContinueMean[i]<<','<<error(msdsContinueMean,msdsContinueSquare,i)<<std::endl;
		//output<<i<<','<<msdsDiscreteMean[i]<<','<<std::sqrt(msdsDiscreteSquare[i]-std::pow(msdsDiscreteMean[i],2))<<','<<msdsContinueMean[i]<<','<<std::sqrt(msdsContinueSquare[i]-std::pow(msdsContinueMean[i],2))<<std::endl;
	}
	
		
	
	rnd.SaveSeed();
	return 0;
}
