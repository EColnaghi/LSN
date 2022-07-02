#ifndef dataBlocking
#define dataBlocking
#include "random.h"
#include "salesman.h"
#include <vector>


//std::vector<SalesmanPath> 	selectBestFit(std::vector<SalesmanPath> paths, Random &rnd, int N);
int world_rank;

void selectBestFit(std::vector<SalesmanPath> &paths, Random &rnd, double p);
void crossOverPop(std::vector<SalesmanPath> &paths);

void pProcess(const char* s);

#endif