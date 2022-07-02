#ifndef __params__
#define __params__
#include "random.h"
#include <vector>


class Params
{
  private:
  double mu;
  double sigma;
  double integral;
  double error;
  double sampleJump;
  int blockLength;
  int nBlocks;
  Random *rnd;
  double hamiltonian();
  void computeIntegral();
  
  public:
  Params(double mu, double sigma, Random *rnd, double sampleJump, double nBlocks, double blockLentgh);
  Params(double mu, double sigma, Random *rnd, double sampleJump, double nBlocks);
  Params(Random *rnd, double nBlocks, double blockLength);
  Params(Random *rnd, double sampleJump, double nBlocks, double blockLength);
  Params();
  void setMu(double mu);
  void setSigma(double sigma);
  void setParams(double mu, double sigma);
  void setRnd(Random *rnd);
  void setSampleJump(double sampleJump);
  void setBlockLength(int blockLength);
  void setnBlocks(int nBlocks);
  void printIntegral(std::vector<double>& integrals, std::vector<double>& errors);
  
  std::vector<double> samplePositions();
  double psi(double position);
  double hamiltonian(double position);
  double getMu();
  double getSigma();
  double getIntegral();
  double getError();
};

#endif
