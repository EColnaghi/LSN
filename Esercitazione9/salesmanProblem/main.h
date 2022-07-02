#ifndef dataBlocking
#define dataBlocking
#include "random.h"
#include "salesman.h"
#include <vector>


//std::vector<SalesmanPath> 	selectBestFit(std::vector<SalesmanPath> paths, Random &rnd, int N);
void selectBestFit(std::vector<SalesmanPath> &paths, Random &rnd, double p);
void crossOverPop(std::vector<SalesmanPath> &paths);

#endif