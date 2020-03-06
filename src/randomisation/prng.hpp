#ifndef __PRNG__
#define __PRNG__

#include<random>

typedef std::mt19937 PRNG;

class prng {
  PRNG gen;
  std::uniform_int_distribution<unsigned int> dis;

public:  
  prng(const unsigned int& MAX_RAND): gen(PRNG(time(NULL))), dis(0,MAX_RAND) {}
  unsigned int rand_max() {return dis.max();}
  unsigned int operator()() { return dis(gen); }
};

#endif
