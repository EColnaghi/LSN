

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

bool isInInterval(double x, double xstart,double xend){
	return (xstart<x&&xend>x)||(xend<x&&xstart>x);
	
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

	//declare variables
	int M=10000;
	int N=100;
	double distance=2;
	double length=1;
	double xstart,ystart,xend,yend;
	int intersect;
	double norm;
	double piEstimate;
	std::vector<double> cumSum;
	std::vector<double> cumSqrSum;
	std::vector<double> avgErrs;
	int crossing=0;
	
	std::ofstream output("buffon.csv");
	output<<"piEstimate,error"<<std::endl;
	
	
	std::ofstream output2("needles.csv");
	output2<<"xstart,ystart,xend,yend,crossing"<<std::endl;
	
	
	for(int j=0;j<N;j++){
		intersect=0;
		for(int i =0;i<M;i++){
			//generate random number between -1 and 1
			xend=rnd.Rannyu(-1,1);
			yend=rnd.Rannyu(-1,1);
			
			//reject the result if the length is bigger than one
			//the point picked will be distributed uniformly in the unitary circle
			while(std::pow(xend,2)+std::pow(yend,2)>1){
				xend=rnd.Rannyu(-1,1);
				yend=rnd.Rannyu(-1,1);
			}
			
			//normalize the vector
			norm=std::sqrt(std::pow(xend,2)+std::pow(yend,2));
			xend*=(length/norm);
			yend*=(length/norm);
			
			//generate a random number in {-1,1}
			//double rndSign=rnd.Rannyu(-1,1);
			//rndSign/=(std::abs(rndSign));
			//use the number generated to flip the vector
			//yend*=rndSign;
			
			xstart=rnd.Rannyu(-xend/2,distance-xend/2);
			ystart=rnd.Rannyu(-yend/2,distance-yend/2);
			//translate the vector
			xend+=xstart;
			yend+=ystart;
			
			//check if the vector crosses the border
			crossing=int(isInInterval(2,xstart,xend)||isInInterval(0,xstart,xend));
			intersect+=crossing;
			output2<<xstart<<","<<ystart<<","<<xend<<","<<yend<<","<<crossing<<std::endl;
		}
		piEstimate=(2*length * M)/(distance*intersect);
		
		cumSum.push_back(piEstimate);
		cumSqrSum.push_back(std::pow(piEstimate,2));
		
	}
	for(int i=1; i<N; i++){
		cumSum[i]+=cumSum[i-1];
		cumSqrSum[i]+=cumSqrSum[i-1];
	}
	
	for(int i=0; i<N; i++){
		cumSum[i]=cumSum[i]/(i+1);
		cumSqrSum[i]=cumSqrSum[i]/(i+1);
		output<<cumSum[i]<<','<<error(cumSum,cumSqrSum,i)<<std::endl;
	}
	
	
	
	rnd.SaveSeed();
	return 0;
}
