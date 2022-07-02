/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/

#include <iostream>
#include <fstream>
#include <cmath>
#include <string>
#include <vector>
#include "random.h"

using namespace std;

double error(vector<double> avg, vector<double> avg2, int n){
	if(n==0){
		return 0;
	}else{
		return sqrt((avg2[n]-pow(avg[n],2))/n);
	}
}

int main (int argc, char *argv[]){

   Random rnd;
   int seed[4];
   int p1, p2;
   ifstream Primes("Primes");
   if (Primes.is_open()){
      Primes >> p1 >> p2 ;
   } else cerr << "PROBLEM: Unable to open Primes" << endl;
   Primes.close();

   ifstream input("seed.in");
   string property;
   if (input.is_open()){
      while ( !input.eof() ){
			input >> property;
			if( property == "RANDOMSEED" ){
				input >> seed[0] >> seed[1] >> seed[2] >> seed[3];
				rnd.SetRandom(seed,p1,p2);
			}
		}
		input.close();
	} else cerr << "PROBLEM: Unable to open seed.in" << endl;

	int M=10000;
	int N=100;
	int L=M/N;	
	
	vector<double> avgs;
	vector<double> avgs2;
	vector<double> avgsErrs;
	vector<double> msds;
	vector<double> msds2;
	vector<double> msdErrs;
	
	double sum=0;
	double sumMSD=0;
	double next;
	for(int j=0; j<N; j++){
		sum=0;
		sumMSD=0;
		for(int i=0; i<L; i++){
		  next=rnd.Rannyu();
		  sum+=next;
		  sumMSD+=(pow(next-0.5,2));
		}
		avgs.push_back(sum/L);
		avgs2.push_back(pow(sum/L,2));
		msds.push_back(sumMSD/L);
		msds2.push_back(pow(sumMSD/L,2));
		
	}
	
	for(int i=1; i<N; i++){
		avgs[i]+=avgs[i-1];
		avgs2[i]+=avgs2[i-1];
		msds[i]+=msds[i-1];
		msds2[i]+=msds2[i-1];
	}
	for(int i=0; i<N; i++){
		avgs[i]=avgs[i]/(i+1);
		avgs2[i]=avgs2[i]/(i+1);
		msds[i]=msds[i]/(i+1);
		msds2[i]=msds2[i]/(i+1);

		avgsErrs.push_back(error(avgs,avgs2,i));
		msdErrs.push_back(error(msds,msds2,i));
	}
	
	
	ofstream output("meanBlocking.csv");
	output<<"avgs"<<","<<"avgs2"<<","<<"avgErrs"<<","<<"msds"<<","<<"msds2"<<","<<"msdErr"<<endl;
	
	for(int i=0; i<avgs.size(); i++){
	  output<<avgs[i]<<','<<avgs2[i]<<","<<avgsErrs[i]<<","<<msds[i]<<","<<msds2[i]<<","<<msdErrs[i]<<endl;
	}
	
	
	
	rnd.SaveSeed();
	return 0;
}

/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/
