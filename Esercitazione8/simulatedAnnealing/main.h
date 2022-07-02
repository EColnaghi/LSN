
#ifndef dataBlocking
#define dataBlocking

//Random numbers
#include "random.h"
#include "params.h"
#include <vector>
int seed[4];
Random rnd;

//functions

double error(std::vector<double> A, std::vector<double> A2, int n);
double annealingStep(Random &rnd, int length, Params params, double x_0, double T, double &integral);
void progError(std::vector<double> &A, std::vector<double> &A2, std::vector<double> &err);
template<typename Function>
double computeIntegral(std::vector<double> positions,Function params);
double temp(double T_0,int i, double exponent);
#endif