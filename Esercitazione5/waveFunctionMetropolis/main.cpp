#include <iostream>
#include <fstream>
#include <algorithm>
#define _USE_MATH_DEFINES
#include <cmath>
#include <string>
#include <vector>
#include "random.h"

struct Point { double x, y,z;};
typedef struct Point Point;


double error(std::vector<double> A, std::vector<double> A2, int n){
	if(n==0){
		return 0;
	}else{
		return std::sqrt((A2[n]-pow(A[n],2))/n);
	}
}

Point averagePoint(std::vector<Point> positions){
	double sumX=0;
	double sumY=0;
	double sumZ=0;
	for(int i =0; i<(int)positions.size(); i++){
		sumX+=positions[i].x;
		sumY+=positions[i].y;
		sumZ+=positions[i].z;
	}
	sumX/=positions.size();
	sumY/=positions.size();
	sumZ/=positions.size();
	Point toReturn={sumX,sumY,sumZ};
	return toReturn;
}

double averageRadius(std::vector<Point> positions){
	double sum=0;
	for(int i =0; i<(int)positions.size(); i++){
		sum+=std::sqrt(std::pow(positions[i].x,2)+std::pow(positions[i].z,2)+std::pow(positions[i].y,2));
	}
	sum/=positions.size();
	return sum;
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

Point jump100(Point p, Random &rnd){
	Point q = p;
	q.x += (2.5*rnd.Rannyu()-1.25);
	q.y += (2.5*rnd.Rannyu()-1.25);
	q.z += (2.5*rnd.Rannyu()-1.25);
	//std::cout<<q.x<<" "<<q.y<<" "<<q.z<<std::endl;
	return q;
}


Point jump210(Point p, Random &rnd){
	Point q = p;
	q.x += (6*rnd.Rannyu()-3);
	q.y += (6*rnd.Rannyu()-3);
	q.z += (6*rnd.Rannyu()-3);
	//std::cout<<q.x<<" "<<q.y<<" "<<q.z<<std::endl;
	return q;
}


Point jump100Gauss(Point p, Random &rnd){
	Point q = p;
	q.x += (rnd.Gauss(0,1));
	q.y += (rnd.Gauss(0,1));
	q.z += (rnd.Gauss(0,1));
	//std::cout<<q.x<<" "<<q.y<<" "<<q.z<<std::endl;
	return q;
}


Point jump210Gauss(Point p, Random &rnd){
	Point q = p;
	q.x += (rnd.Gauss(0,2));
	q.y += (rnd.Gauss(0,2));
	q.z += (rnd.Gauss(0,2));
	//std::cout<<q.x<<" "<<q.y<<" "<<q.z<<std::endl;
	return q;
}
double psi100(Point p){
	double r=std::sqrt(std::pow(p.x,2)+std::pow(p.y,2)+std::pow(p.z,2));
	//std::cout<<std::pow((1/std::sqrt(M_PI))*std::exp(-r),2)<<std::endl;
	return std::pow((1/std::sqrt(M_PI))*std::exp(-r),2);
}


double psi210(Point p){
	double r     = std::sqrt(std::pow(p.x,2)+std::pow(p.y,2)+std::pow(p.z,2));
	double theta = std::atan((std::sqrt(std::pow(p.x,2)+std::pow(p.y,2)))/p.z);
	return std::pow((1.0/8.0)*std::sqrt(2.0/M_PI)*r*std::exp(-r/2.0)*cos(theta),2);
}


/*this function fill the vectors xs with points sampled from the distribution pdf using transition probability jump, 
starting from point x_0. Length is the number of points we want and rnd is the Random object we use to generate
the numbers*/
template<typename Function, typename Function2>
void metropolisSampling3D(std::vector<Point> &xs, Function pdf, Function2 jump, Point x_0, int length, Random &rnd){
	
	Point currentPos,candidate;
	double probCandidate,probCurrent;
	double reject;
	double treshold;
	probCurrent=pdf(x_0);
	currentPos=x_0;
	int rejectionNumbers=0;
	
	for(int i =0; i<length;i++){
		candidate = jump(currentPos, rnd);
		probCandidate=pdf(candidate);
		treshold=std::min(1.0,(probCandidate/probCurrent));
		reject=rnd.Rannyu();
		//std::cout<<"currentProb, candidateProb: "<<probCurrent<<" " << probCandidate<<std::endl;;
		if(reject<treshold){
			currentPos=candidate;
			probCurrent=pdf(currentPos);
		}else{
			rejectionNumbers++;
		}
		xs.push_back(currentPos);	
	}
	std::cout<<"----------------------"<<std::endl;
	std::cout<<"rejection rate: "<<(double)rejectionNumbers/length<<std::endl;
	std::cout<<"----------------------"<<std::endl;
}

/*this function fill the vectors avgRadius, avgRadiusSqr and err with the progressive means and
the estimated errors for N blocks of length length, starting from point xs0*/

template<typename Function, typename Function2>
void computeMeans(std::vector<double> &avgRadius, std::vector<double> &avgRadiusSqr, std::vector<double> &err, Point xs0, int length, int N, Function pdf, Function2 jump, Random rnd){

	avgRadius.clear();
	avgRadiusSqr.clear();
	errors.clear(); 

	std::vector<Point> avgPos;
	std::vector<Point> positions;
	
	
	//metropolisSampling3D(xs,pdf,jump,x_0,length,rnd){
	for(int i = 0;i<N;i++){
		metropolisSampling3D(positions, pdf, jump, xs0, length, rnd);
		double temp=averageRadius(positions);
		avgRadius.push_back(temp);
		avgRadiusSqr.push_back(std::pow(temp,2));
		xs0=positions.back();
		positions.clear();
	}
	progError(avgRadius,avgRadiusSqr,err);
	
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
	double x0=0;
	double y0=0;
	double z0=15;
	int length=1000;
	double N=100;
	std::vector<Point> avgPos;
	std::vector<double> avgRadius;
	std::vector<double> avgRadiusSqr;
	std::vector<double> errors;
	
	std::vector<Point> positions;
	
	Point xs0={x0,y0,z0};
	
	std::ofstream output("psi100.csv");
	output<<"r2,r2^2,err"<<std::endl;
	std::ofstream output2("psi210Averages.csv");
	output2<<"r2,r2^2,err"<<std::endl;
	std::ofstream output3("psi100AveragesGaussianJumps.csv");
	output3<<"r2,r2^2,err"<<std::endl;	
	std::ofstream output4("psi210AveragesGaussianJumps.csv");
	output4<<"r2,r2^2,err"<<std::endl;
	
	//computate the averages, the errors, and write them to the output files
	computeMeans(avgRadius,avgRadiusSqr,errors,xs0,length,N,psi100,jump100,rnd);
	for(int i =0; i< N;i++){output<<avgRadius[i]<<','<<avgRadiusSqr[i]<<','<<errors[i]<<std::endl;}
	computeMeans(avgRadius,avgRadiusSqr,errors,xs0,length,N,psi210,jump210,rnd);
	for(int i =0; i< N;i++){output2<<avgRadius[i]<<','<<avgRadiusSqr[i]<<','<<errors[i]<<std::endl;}
	computeMeans(avgRadius,avgRadiusSqr,errors,xs0,length,N,psi100,jump100Gauss,rnd);
	for(int i =0; i< N;i++){output3<<avgRadius[i]<<','<<avgRadiusSqr[i]<<','<<errors[i]<<std::endl;}		
	computeMeans(avgRadius,avgRadiusSqr,errors,xs0,length,N,psi210,jump210Gauss,rnd);
	for(int i =0; i< N;i++){output4<<avgRadius[i]<<','<<avgRadiusSqr[i]<<','<<errors[i]<<std::endl;}		
	
	rnd.SaveSeed();
	return 0;
}
