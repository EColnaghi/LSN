#ifndef __Blocking__
#define __Blocking__

class Blocking {

private:
  int ;

protected:

public:
  // constructors
  Random();
  // destructor
  ~Random();
  // methods
  void SetRandom(int * , int, int);
  void SaveSeed();
  double Rannyu(void);
  double Rannyu(double min, double max);
  double Gauss(double mean, double sigma);
  double cauchy(double mu, double gamma);
  double exponential(double lambda);
};

#endif // __Random__