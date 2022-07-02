#include <math.h>
#include <iostream>
#include <fstream>
#include <algorithm>
#define _USE_MATH_DEFINES
#include <string>
#include <vector>
#include <iomanip>
#include "random.h"
#include "salesman.h"
#include "main.h"
#include <mpi.h>



int main (int argc, char *argv[]){
	int world_size;
    // Initialize the MPI environment
    MPI_Init(&argc, &argv);
	// Get the number of processes
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);
    // Get the rank of the process
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
	Random rnd;
	int seed[4];
	int p1, p2;
	std::ifstream Primes("Primes");
	for(int i = 0; i < world_rank;i++){
		if (Primes.is_open()){
			Primes >> p1 >> p2 ;
			std::cout<<"per il processo "<<world_rank<<" leggo i primi "<<p1<<","<<p2<<std::endl;
		} else std::cerr << "PROBLEM: Unable to open Primes" << std::endl;
	}
	if (Primes.is_open()){
		Primes >> p1 >> p2 ;
		std::cout<<"per il processo "<<world_rank<<" leggo i primi "<<p1<<","<<p2<<std::endl;
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

	
    // Print off a hello world message
	std::cout<<"Hello world from process rank "<< world_rank <<" out of "<< world_size <<" processors\n"<<std::endl;
	
	
	//initialize variables
	double cityX,cityY;
	int size=500;			//number of cities
	int nMigr=10;			//number of generations between migrations
	int nExchange=1;		//number of continents that migrates at every migration event
	int popMigrating=3;		//number of individuals that migrates at each migration
	int sizePop=1000;		//size of the population
	int nGens=1000;			//number of generation
	double pChoice = 0.5;	//related to the probability of reproduction
	City city(0.0,0.0);		
	double defaultProb=0.1;//default probability of mutation
	double pCrossOver=0.5;	//probability of crossover
	double mean=0;
	int readInputFile = 1;	//if 1 read the input file American_capitals.dat, otherwise generate cities 
	//randomly change the parameters, except for the process 0
	/*if(world_rank!=0){
		defaultProb+=rnd.Rannyu(-0.05,0.35);
		pCrossOver +=rnd.Rannyu(-0.2,0.2);
		pChoice    +=rnd.Rannyu(-0.2,0.2);
	}*/
	
	std::vector<City> cities;
	std::vector<SalesmanPath> pop;
	std::vector<int> path;
	
	
	//Intitialize the output
	std::ofstream outputCity("city.csv");
	outputCity<<"x,y"<<std::endl;
	std::string outputPathStr="path";
	outputPathStr+=std::to_string(world_rank);
	outputPathStr+=".csv";
	std::cout<<"Process "<<world_rank<<" writes its data in file"<<outputPathStr<<std::endl;
	std::ofstream outputPath(outputPathStr);
	
	std::string outputParamsStr="params";
	outputParamsStr+=std::to_string(world_rank);
	outputParamsStr+=".csv";
	std::cout<<"Process "<<world_rank<<" writes its parameters in file"<<outputParamsStr<<std::endl;
	std::ofstream outputParams(outputParamsStr);
	
	outputPath<<"gen,path,length,meanLength"<<std::endl;
	outputParams<<"size,nMigr,nExchange,popMigrating,sizePop,nGens,pchoice,pmutation,pcrossover"<<std::endl;
	outputParams<<size<<","<<nMigr<<","<<nExchange<<","<<popMigrating<<","<<sizePop<<","<<nGens<<","<<pChoice<<","<<defaultProb<<","<<pCrossOver<<std::endl;
	
	//Initialize the input
	std::ifstream inputCity("American_capitals.dat");
	std::string s;
	
	//initialize MPI variables
	int* pathToSend = new int[size];
	int* pathReceive;
	int sender, receiver;
	MPI_Status stat;
	pathReceive = new int[size*world_size];
	int itag=1;
	

	
	//read the input files
	//UNCOMMENT THIS FOR US CAPITALS
	
	//process 0 generate the city and send them to the others
	//UNCOMMENT THIS FOR RANDOMLY GENERATED CITIES

	if(readInputFile==1){size=50;}	
	for(int i = 0; i<size; i++){
		

		if(world_rank==0){
			if(readInputFile==1){
				inputCity>>cityX;
				inputCity>>cityY;
				city= City(cityX,cityY);
				cities.push_back(city);
				std::cout<<world_rank<<" :genero la citta "<<i<<" con x: "<<cityX<<" e y: "<<cityY<<std::endl;
			}else{
				cityX=rnd.Rannyu(-1,1);
				cityY=rnd.Rannyu(-1,1);
			}
		}
		
		
		MPI_Bcast(&cityX,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
		MPI_Bcast(&cityY,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
		city= City(cityX,cityY);
		cities.push_back(city);
	}
	
	
	for(int i = 0; i<size; i++){
		path.push_back(i);
	}
	
	//The process 0 print the cities list in the output file	
	/*if(world_rank==0){
		
		for(int i =0; i< size; i++){
			outputCity<<cities[i].x<<","<<cities[i].y<<std::endl;
		}
	}*/
	//initialize the population
	for(int i = 0; i<sizePop; i++){
		
		std::random_shuffle(path.begin()+1, path.end());
		std::vector<int> tempPath =path;
		SalesmanPath sp(tempPath,cities,rnd, defaultProb, pCrossOver);
		pop.push_back(sp);
	}
	std::sort(pop.begin(),pop.end());
	std::cout<<"il processo "<<world_rank<<" ha un vector per la popolazione di: "<<pop.size()<<std::endl;
	
	//main cycle: do nGens iterations of combining, exchanging and reproducing solutions
	for(int i =0; i<nGens; i++){
		mean=0;
		//sort, reproduce the best and mutate
		std::sort(pop.begin(),pop.end());
		selectBestFit(pop,rnd, pChoice);
		for(auto &k:pop){
			int index=std::floor(rnd.Rannyu()*size);
			k.mutate(pop[index]);
		}
		std::sort(pop.begin(),pop.end());
		//print the next generation
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
		if(i%nMigr==0){
			
			
			/*for(int l=0; l<nExchange;l++){
				
				//process 0 decides who sends and who receives
				if(world_rank==0){
					sender=(int)std::floor(rnd.Rannyu()*world_size);
					receiver=(int)std::floor(rnd.Rannyu()*world_size);
					while(receiver==sender){
						receiver=(int)std::floor(rnd.Rannyu()*world_size);
					}
					
					std::cout<<"Il sender dovrebbe essere "<<sender<<" mentre il receiver "<<receiver<<std::endl;
				}
				
				//process 0 broadcast the id of the sender and of the receiver
				MPI_Bcast(&sender, 1, MPI_INTEGER, 0, MPI_COMM_WORLD);
				MPI_Bcast(&receiver, 1, MPI_INTEGER, 0, MPI_COMM_WORLD);
				
				
				for(int k=0; k<popMigrating;k++){
					//the receiver receive the path and put it in the last positions of the array
					if(world_rank==receiver){	
						std::cout<<receiver<<": Devo ricevere il path"<<std::endl;				
						MPI_Recv(pathReceive, size, MPI_INTEGER, sender,itag,MPI_COMM_WORLD, &stat);
						
						std::cout<<receiver<<": Ho ricevuto il path, i primi tre elementi sono: "<<pathReceive[1]<<","<<pathReceive[2]<<","<<pathReceive[3]<<std::endl;	
						std::cout<<receiver<<": Ricevuto il path da "<<sender<<" l'ho copiato nel vettore"<<std::endl;
						
						//I can't  use the setPath method, because MPI_Recv store the data in an array
						//std::copy(pathReceive, pathReceive+size, pop[sizePop-k-1].path.begin());
						//pop[sizePop-k].computeLength();
						pop[sizePop-k-1].setPath(pathReceive);
					}
					
					//tjhe sender sends his k-th best path
					if(world_rank==sender){
						std::cout<<world_rank<<": Devo inviare il path"<<std::endl;	
						std::copy(pop[k].path.begin(), pop[k].path.end(), pathToSend);
						MPI_Send(pathToSend, size, MPI_INTEGER, receiver,itag, MPI_COMM_WORLD);
						std::cout<<world_rank<<": Ho inviato il path, i cui primi tre elementi sono: "<<pathToSend[1]<<","<<pathToSend[2]<<","<<pathToSend[3]<<std::endl;
						std::cout<<world_rank<<": mentre i primi tre elementi del primo path sono:   "<<pop[0].path[1]<<","<<pop[0].path[2]<<","<<pop[0].path[3]<<std::endl;
					}
				}
			}*/
		}
		//std::cout<<"Nuova generazione"<<std::endl;
	}
	
    MPI_Finalize();	
	
	
	return 0;
}

void selectBestFit(std::vector<SalesmanPath> &paths, Random &rnd, double p){
	int M = paths.size();
	int j;
	std::vector<std::vector<int>> tempPaths;
	for(int i=0; i<M; i++){

		
		j=(int)((M-1)*(1-std::pow(rnd.Rannyu(),p)));
		//if(world_rank==0){std::cout<<tempPaths.size()<<" "<<paths[j].path.size()<<" "<<tempPaths.capacity()<<std::endl;
		//for(int k=0;k<paths[j].path.size();k++){std::cout<<paths[j].path[k]<<",";}std::cout<<std::endl;}
		tempPaths.push_back(paths[j].path);
	}
	for(int i=0; i<M; i++){
		//std::cout<<std::endl;
		paths[i].setPath(tempPaths[i]);
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






void pProcess(const char* s){
	if(world_rank==0){
	std::cout<<"Processo "<<world_rank<<" "<<s<<std::endl;
	}
}