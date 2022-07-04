

#include <iostream>
#include <fstream>
#include <algorithm>
#define _USE_MATH_DEFINES
#include <cmath>
#include <string>
#include <vector>
#include <iomanip>
#include "random.h"
#include "salesman.h"
#include "main.h"


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


	
	
	
	
	double random1,random2;
	int size=34;					//number of cities
	int sizePop=100;				//size of the population
	int nGens=1000;				//number of generations
	double pChoice = 0.1;			//related to probability of the selection of the fittest
	City city(0.0,0.0);
	double defaultProb=0.1;			//probability of mutation
	double pCrossOver=0.45;			//probability of crossing over
	double mean=0;
	
	//initialize the output
	std::ofstream outputCity("city.csv");
	outputCity<<"x,y"<<std::endl;
	std::ofstream outputPath("path.csv");
	outputPath<<"gen,path,length,meanLength"<<std::endl;
	std::vector<City> cities;
	std::vector<SalesmanPath> pop;
	std::vector<int> path;
	for(int i = 0; i<size; i++){
		path.push_back(i);
	}
	
	//fill the vector with cities in a square
	for(int i = 0; i<size; i++){
		random1=rnd.Rannyu(-1,1);
		random2=rnd.Rannyu(-1,1);
		city= City(random1,random2);
		cities.push_back(city);
	}
	
	//fill the vector with cities on a circle
	/*for(int i = 0; i<size; i++){
		random1=2*rnd.Rannyu()*M_PI;
		city= City(std::cos(random1),std::sin(random1));
		cities.push_back(city);
	}*/
	//print the cities in the ouput file
	for(int i =0; i< size; i++){
		outputCity<<cities[i].x<<","<<cities[i].y<<std::endl;
	}
	
	//populate the vector of solutions
	for(int i = 0; i<sizePop; i++){
		std::random_shuffle(path.begin()+1, path.end());
		std::vector<int> tempPath =path;
		SalesmanPath sp(tempPath,cities,rnd, defaultProb, pCrossOver);
		pop.push_back(sp);
	}
	
	//sort the solutions by length
	std::sort(pop.begin(),pop.end());
	
	//main cycle
	for(int i =0; i<nGens; i++){
		mean=0;
		
		//select the fittest
		selectBestFit(pop,rnd, pChoice);
		for(auto &k:pop){
			//mutate and cross over
			int index=std::floor(rnd.Rannyu()*size);
			k.mutate(pop[index]);
		}
		
		//sort the mutated population
		std::sort(pop.begin(),pop.end());
	
		outputPath<<i<<",";
		for(int j=0; j<size; j++){
			outputPath<<pop[0].path[j]<<';';
		}
		
		outputPath<<"0,"<<pop[0].length;
		//print the mean in the output fil
		for(int j=0; j<size/2; j++){
			mean+=pop[j].length;
		}
		outputPath<<","<<mean/(size/2)<<std::endl;
		
		std::cout<<"Generazione: "<<i<<" best path length: "<<pop[0].length<<std::endl<<std::endl;
	}
	
	
		
	return 0;
}

void selectBestFit(std::vector<SalesmanPath> &paths, Random &rnd, double p){
	//std::cout<<"entro in selectBestFit"<<std::endl;
	int M = paths.size();
	//std::vector<SalesmanPath> tempPaths;
	int j;
	std::vector<std::vector<int>> tempPaths;
	for(int i=0; i<M; i++){
		j=(int)((M-1)*(1-std::pow(rnd.Rannyu(),p)));
		tempPaths.push_back(paths[j].path);
	}
	for(int i=0; i<M; i++){

		//std::cout<<std::endl;
		//paths[i].setPath(tempPaths[i]);
		paths[i].path=tempPaths[i];
		paths[i].computeLength();
	}
}

void crossOverPop(std::vector<SalesmanPath> &paths, Random &rnd){
	int size=paths.size();
	int index1=std::floor(rnd.Rannyu()*size);
	int index2=std::floor(rnd.Rannyu()*size);
	paths[index1].crossOver(paths[index2]);
	
}

double error(std::vector<double> A, std::vector<double> A2, int n){
	if(n==0){
		return 0;
	}else{
		return std::sqrt((A2[n]-pow(A[n],2))/n);
	}
}







