#ifndef AUX_MATH
#define AUX_MATH
namespace aux_math {
  long double fact(const unsigned int& k) {
    if(k==0)
      return 1;
    else 
      return k*fact(k-1); // Tail recursion!
  }
  
  long double number_of_ordered_subsets(const unsigned int& n, const unsigned int& k) {
    if(k<=n) {
      long double numer = 1;
      for(unsigned int i=0; i<k; ++i)
        numer*=n-i;
      return numer;
    }
    return 0;
  }

  long double binom(const unsigned int& n, const unsigned int& k) {
    if(k<=n) {
      long double numer = 1;
      for(unsigned int i=0; i<k; ++i)
        numer*=(n-i)/(double)(i+1);
      return numer;
    }
    return 0;
  }  
}
#endif
