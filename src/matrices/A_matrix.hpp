#ifndef __A_MATRIX_HPP__
#define __A_MATRIX_HPP__

#include "symm_matrix.hpp"
#include <fstream>
#include <math.h>
#include <iostream>
#include <string>

class A_matrix : public symm_matrix<bool> {//Adjacency matrix
  public:
  A_matrix();
  A_matrix(const char*); //Loading from file
  A_matrix(const unsigned int&);
  A_matrix(const A_matrix&);

  bool* operator[](const unsigned int &i) const {return &array[dim_x * i];};
  const A_matrix& operator=(const A_matrix&);

  unsigned int size() const {return dim_x;};
  bool check_consistency() const;
  unsigned int num_links() const;

  void save(const std::string&) const;

  friend std::ostream& operator<<(std::ostream&, const A_matrix&);

  A_matrix& GB_Metropolis_generator(double(&Hamiltonian)(const A_matrix&), const unsigned int &N_iters, const int &seed, const bool &initialize_randomly=true, const double &temp=1); //Generates random graph from the ensemble defined by Gibbs-Boltzmann distribution, which Hamiltonian is passed as a function of the adjacency matrix.
  A_matrix& MF_GB_Metropolis_generator(double(&H)(const unsigned int& N_links), const unsigned int &N_iters, const int &seed, const bool &initialize_randomly=true, const double &temp=1); //Generates random graph from the ensemble defined by Gibbs-Boltzmann distribution, which Hamiltonian is passed as a function of the adjacency matrix.
};

#endif
