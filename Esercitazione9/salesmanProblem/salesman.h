#ifndef Salesman
#define Salesman
#include "random.h"
#include <vector>

class City{
	public:
		double x;
		double y;
	
	City(double x, double y);
	double distance(City c2);	
};

class SalesmanPath{
	Random *rnd;
	int size;
	double p_pp, p_s, p_i, p_bp, p_co;//probabilities of pair permutation, shifting, inversion and contiguous block permutation
	std::vector<City> cities;
	
	public:
		std::vector<int> path;
		double length;
		
			
		SalesmanPath();
		SalesmanPath(std::vector<int> &path, std::vector<City> &cities, Random& rnd);
		SalesmanPath(std::vector<int> &path, std::vector<City> &cities, Random& rnd, double defaultp);
		SalesmanPath(std::vector<int> &path, std::vector<City> &cities, Random& rnd, double defaultp, double p_co);
		SalesmanPath(std::vector<int> &path, std::vector<City> &cities, Random& rnd, double pp, double pi, double ps, double pc, double p_co);
		SalesmanPath(const SalesmanPath& old);
			
		
		friend bool operator< (const SalesmanPath& c1, const SalesmanPath& c2);
		friend bool operator> (const SalesmanPath& c1, const SalesmanPath& c2);
		bool check();
		void pairPermutation();
		void contiguousPermutation();
		void inversion();
		void shift();
		void crossOver(SalesmanPath &s);
		void computeLength();
		void print();
		void mutate(SalesmanPath& c2);
		void mutate();
		void setPath(std::vector<int> path);
		std::vector<int> setPath();

};

#endif