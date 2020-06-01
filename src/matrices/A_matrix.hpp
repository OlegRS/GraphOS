#ifndef __A_MATRIX_HPP__
#define __A_MATRIX_HPP__

#include "symm_matrix.hpp"
#include "col_vector.hpp"
#include "../randomisation/prng.hpp"
#include <fstream>
#include <math.h>
#include <iostream>
#include <string>
#include <functional>

class A_matrix : public symm_matrix<bool> {//Adjacency matrix
  public:
  A_matrix();
  A_matrix(const char*); //Loading from file
  A_matrix(const unsigned int&);
  A_matrix(const A_matrix&);

  bool* operator[](const unsigned int &i) const {return &array[dim_x * i];};
  const A_matrix& operator=(const A_matrix&);
  bool check_consistency() const;
  
  void save(const std::string&) const;

  A_matrix& set_Erdos_Renyi(const double &p, prng& rnd);

  col_vector<col_vector<unsigned int> > random_node_pairs_col_vector(const unsigned int &N_pairs, prng& rnd) const;

  // GENERATION
  A_matrix& single_link_GB_Metropolis_generator(double(&Hamiltonian)(const A_matrix&), const unsigned int &N_iters, prng& rnd, const bool &initialize_randomly=true, const double &temp=1); //Generates random graph from the ensemble defined by Gibbs-Boltzmann distribution, which Hamiltonian is passed as a function of the adjacency matrix.
  A_matrix& single_link_MF_GB_Metropolis_generator(double(&H)(const unsigned int& N_links), const unsigned int &N_iters, prng& rnd, const bool &initialize_randomly=true, const double &temp=1); //Generates random graph from the ensemble defined by Gibbs-Boltzmann distribution, which Hamiltonian is passed as a function of the adjacency matrix.
  A_matrix& MF_GB_Metropolis_generator(double(&H)(const unsigned int& N_links), const unsigned int &N_iters, prng& rnd, const unsigned int &N_pairs_max=1, const bool &initialize_randomly=true, const double &temp=1); //Generates random graph from the ensemble defined by Gibbs-Boltzmann distribution, which Hamiltonian is passed as a function of the adjacency matrix.
  A_matrix& sample_p_star_model(const unsigned int &N_iters, prng& rnd, const col_vector<double> &T, const unsigned int &N_pairs_max=1, const bool &initialize_randomly=true, const double &temp=1);
  A_matrix& sample_p_star_model_with_single_link_Metropolis(const unsigned int &N_iters, prng& rnd, const col_vector<double> &T, const bool &initialize_randomly=true, const double &temp=1);
  A_matrix& sample_triad_model_with_single_link_Metropolis(const unsigned int &N_iters, prng& rnd, const double h, const double sigma, const double tau, const bool &initialize_randomly=true, const double &temp=1);

  // ANALYSIS
  unsigned int size() const {return dim_x;};
  unsigned int num_links() const;
  unsigned int degree(const unsigned int& i) const;
  col_vector<unsigned int> degree_sequence_col_vec() const;
  col_vector<long double> num_p_stars_sequence_col_vec(const unsigned int& p) const;
  long double average_p_stars(const unsigned int& p) const;
  long double num_p_stars(const unsigned int& node, const unsigned int& p) const;
  unsigned long long num_triangles(const unsigned int& node) const; // Number of triangles with a given node
  unsigned long long num_triangles() const; // Full number of triangles in the network (INEFFICIENT!!!)
  double loc_clust_coeff(const unsigned int& node) const;
  col_vector<double> loc_clust_coeff_col_vec() const;
  
  friend std::ostream& operator<<(std::ostream&, const A_matrix&);
};

#endif
