#include <vector>
#include "random.h"
#ifndef __metropolis__
#define __metropolis__

template<typename Function, typename Function2, typename Point>
void metropolisSampling(std::vector<Point> &xs, Function pdf, Function2 jump, Point x_0, int length, Random &rnd){
	double probCandidate,probCurrent;
	double reject;
	double treshold;
	probCurrent=pdf(x_0);
	Point currentPos=x_0;
	Point candidate=x_0;
	int rejectionNumbers=0;
	
	for(int i =0; i<length;i++){
		candidate = jump(currentPos);
		probCandidate=pdf(candidate);
		treshold=(probCandidate/probCurrent);
		reject=rnd.Rannyu();
		if(reject<treshold){
			currentPos=candidate;
			probCurrent=pdf(currentPos);
		}else{
			rejectionNumbers++;
		}
		xs.push_back(currentPos);
	}
}

template<typename Function, typename Function2, typename Point>
void metropolisSamplingVerbose(std::vector<Point> &xs, Function pdf, Function2 jump, Point x_0, int length, Random &rnd){
	double probCandidate,probCurrent;
	double reject;
	double treshold;
	probCurrent=pdf(x_0);
	Point currentPos=x_0;
	Point candidate=x_0;
	int rejectionNumbers=0;
	
	for(int i =0; i<length;i++){
		candidate = jump(currentPos);
		probCandidate=pdf(candidate);
		std::cout<<"candidate position: "<<candidate.getMu()<<" "<<candidate.getSigma()<<" "<<candidate.getIntegral()<<" prob "<< probCandidate<<std::endl;
		std::cout<<"Current position: "<<currentPos.getMu()<<" "<<currentPos.getSigma()<<" "<<currentPos.getIntegral()<<" prob "<< probCurrent<<std::endl;
		treshold=(probCandidate/probCurrent);
		reject=rnd.Rannyu();
		
		std::cout<<"Treshold: "<<treshold<<" reject: "<<reject<<std::endl;
		if(reject<treshold){
			currentPos=candidate;
			probCurrent=probCandidate;
		}else{
			rejectionNumbers++;
		}
		xs.push_back(currentPos);
	}
	//std::cout<<"Probability current: "<<probCurrent<<"  probability candidate: "<<probCandidate<<std::endl;
}

#endif