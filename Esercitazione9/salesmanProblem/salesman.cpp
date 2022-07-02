#include <iostream>
#include <fstream>
#include <algorithm>
#define _USE_MATH_DEFINES
#include <cmath>
#include <string>
#include <vector>
#include <iomanip>
#include <list>
#include "salesman.h"
#include "random.h"

SalesmanPath::SalesmanPath(){
	
	
}

SalesmanPath::SalesmanPath(std::vector<int> &path, std::vector<City> &cities, Random &rnd){

	this->path      = path;
	this->cities 	= cities;
	this->size      = (int) path.size();
	this->computeLength();
	this->rnd=&rnd;
	this->p_pp=0.1;
	this->p_i=0.1;
	this->p_s=0.1;
	this->p_bp=0.1;
	this->p_co=0.1;
	}

SalesmanPath::SalesmanPath(std::vector<int> &path, std::vector<City> &cities, Random &rnd, double defaultP){

	this->path      = path;
	this->cities 	= cities;
	this->size      = (int) path.size();
	this->computeLength();
	this->rnd=&rnd;
	this->p_pp=defaultP;
	this->p_i=defaultP;
	this->p_s=defaultP;
	this->p_bp=defaultP;
	this->p_co=defaultP;
}


SalesmanPath::SalesmanPath(std::vector<int> &path, std::vector<City> &cities, Random &rnd, double defaultP, double p_co){

	this->path      = path;
	this->cities 	= cities;
	this->size      = (int) path.size();
	this->computeLength();
	this->rnd=&rnd;
	this->p_pp=defaultP;
	this->p_i=defaultP;
	this->p_s=defaultP;
	this->p_bp=defaultP;
	this->p_co=p_co;
}

	
SalesmanPath::SalesmanPath(const SalesmanPath& old){
	
	this->path      = old.path;
	this->cities 	= old.cities;
	this->size      = (int) path.size();
	this->computeLength();
	this->rnd		=old.rnd;
	this->p_pp		=old.p_pp;
	this->p_i		=old.p_i;
	this->p_s		=old.p_s;
	this->p_bp		=old.p_bp;
	this->p_co		=old.p_co;
}

void SalesmanPath::setPath(std::vector<int> path){
	this->path=path;
	this->computeLength();
}

void SalesmanPath::computeLength(){
	double length=0;
	/*for(auto k: this->path){
		std::cout<<k<<" ";
	}
	std::cout<<std::endl;*/
	int size=this->size;
	for(int i =0; i<size; i++)
	{
		//std::cout<<"Entro in distance per "<<path[i]<<" e "<<path[(i+1)]%size<<std::endl;
		
		//std::cout<<"Entro in distance per "<<&(this->cities[path[i]])<<" e "<<&cities[path[(i+1)]%size]<<std::endl;
		double toSum=this->cities[path[i]].distance(cities[path[(i+1)%size]]);
		length+=toSum;
	}
	this->length=length;
}

void SalesmanPath::mutate(){
	//std::cout<<"entro in mutate per il path "<<this<<std::endl;
	double nextRand=rnd->Rannyu();
	
	if(nextRand<this->p_pp){
		this->pairPermutation();
	
		//std::cout<<"pp";
	}
	
	nextRand=rnd->Rannyu();
	if(nextRand<this->p_i){
		this->inversion();
		
		//std::cout<<"i";
	}
	
	nextRand=rnd->Rannyu();
	if(nextRand<this->p_s){
		this->shift();
		//std::cout<<"s";
	}
	
	nextRand=rnd->Rannyu();
	if(nextRand<this->p_bp){
		this->contiguousPermutation();
		//std::cout<<"cp";
	}
	
	
}



void SalesmanPath::mutate(SalesmanPath& c2){
	double nextRand=rnd->Rannyu();
	
	if(nextRand<this->p_pp){
		this->pairPermutation();
	}
	
	nextRand=rnd->Rannyu();
	if(nextRand<this->p_i){
		this->inversion();
	}
	
	nextRand=rnd->Rannyu();
	if(nextRand<this->p_s){
		this->shift();
	}
	
	nextRand=rnd->Rannyu();
	if(nextRand<this->p_bp){
		this->contiguousPermutation();
	}
	
	if(nextRand<this->p_co){
		this->crossOver(c2);
	
	}
		
}

void SalesmanPath::pairPermutation(){
	double random1 =this->rnd->Rannyu();
	double random2 =this->rnd->Rannyu();
	int toPermute1= std::floor(1+(random1*(this->size-1)));
	int toPermute2= std::floor(1+(random2*(this->size-1)));
	int temp= this->path[toPermute1];
	this->path[toPermute1]=this->path[toPermute2];
	this->path[toPermute2]=temp;
	
	//std::cout<<"permuto la coppia "<<toPermute1<<" e "<<toPermute2<<std::endl;
	
	this->computeLength();
}

//exchange m contiguous cities, starting at position N, with othe m contiguous cities
void SalesmanPath::contiguousPermutation(){
	
	int m = std::floor( this->rnd->Rannyu()*(this->size/2));
	int N =std::floor( 1+this->rnd->Rannyu()*(this->size-m));
	int n =std::floor( 1+this->rnd->Rannyu()*(this->size-m));

	std::swap_ranges(this->path.begin()+N,this->path.begin()+N+m,this->path.begin()+n);
	
	this->computeLength();
}

//shift m contiguous cities, starting from N, of n spaces
void SalesmanPath::shift(){
	double m =std::floor( 1+(this->rnd->Rannyu()*this->size-1));
	int n =std::floor( 1+this->rnd->Rannyu()*(this->size-m));
	int N =std::floor( 1+this->rnd->Rannyu()*(this->size-m));
	
	int temp;
	for(int i = 0; i<m;i++){
		temp=this->path[i+N];
		this->path[i+N]=this->path[i+n];
		this->path[i+n]=temp;
	}
	
	this->computeLength();
}

void SalesmanPath::inversion(){
	int n,m;
	m=std::floor(1+(this->size-1)*rnd->Rannyu());
	n=std::floor(1+(this->size-m)*rnd->Rannyu());
	std::reverse(this->path.begin()+m,this->path.begin()+m+n);
	
	this->computeLength();
}

void SalesmanPath::print(){
	std::cout<<" Path corrente: ";
	for(int i=0;i<this->size;i++){
		std::cout<<"  "<<this->path[i];
	}
	std::cout<<"   Lunghezza:   "<<this->length<<std::endl;
}

void SalesmanPath::crossOver(SalesmanPath &p){
	int n,m;
	m=std::floor(this->size*rnd->Rannyu());
	n=std::floor(this->size*rnd->Rannyu());
	std::vector<int> tempPathThis;
	std::vector<int> tempPath; 	
	std::vector<int> toCopy;	
	
	
	if(m>n){
		int temp=m;
		m=n;
		n=temp;
	}
	
	std::copy(p.path.begin()+m,p.path.begin()+n,std::back_inserter(tempPath));
	std::copy(this->path.begin()+m,this->path.begin()+n,std::back_inserter(tempPathThis));	
	
	for(int i =0; i<(int)p.path.size(); i++){
		for(int j=0; j<(int)tempPathThis.size();j++){
			if(p.path[i]==tempPathThis[j]){
				toCopy.push_back(p.path[i]);
			}
		}
	}
	std::copy(toCopy.begin(),toCopy.end(),this->path.begin()+m);
	this->computeLength();
}

bool SalesmanPath::check(){
	bool flag=true;
	std::vector<int> temp=this->path;
	sort(temp.begin(), temp.end());
	for(int i =0;i<size;i++){
		if(i!=temp[i]){flag=false;}
	}
	return flag;
}

bool operator< (const SalesmanPath& c1, const SalesmanPath& c2)
{
    return c1.length < c2.length;
}

bool operator> (const SalesmanPath& c1, const SalesmanPath& c2)
{
    return c1.length > c2.length;
}

City::City(double x, double y){
	this->x=x;
	this->y=y;
}

double City::distance(City c2){
	
	return(std::pow(this->x-c2.x,2)+std::pow(this->y-c2.y,2));
	
}