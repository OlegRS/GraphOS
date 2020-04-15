#ifndef __MATRIX_HPP__
#define __MATRIX_HPP__

#include <string>
#include <cmath>
#include <iostream>
#include <fstream>

template <typename T> class matrix {
protected:
  unsigned int dim_y;
  unsigned int dim_x;

  T* array;
public:
  matrix();
  matrix(const matrix<T>&);
  matrix(const unsigned int&, const unsigned int&);
  matrix(const unsigned int&); //Create square matrix
  matrix(const std::string&); //Load matrix from file

  const unsigned int& get_dim_x() const { return dim_x; }
  const unsigned int& get_dim_y() const { return dim_y; }

  T* operator[](const unsigned int&) const;
  matrix<T> &operator+=(const matrix&);
  matrix<T> &operator/=(double);
  const matrix<T>& operator=(const matrix<T>&);

  T tr() const; //Returns trace of the matrix
  
  const matrix<T>& abs(); //returns matrix of absolute values of all elements
  
  T max_entrie() const;
  matrix<T> max_n_entries(const unsigned int&) const;
  const matrix<T>& keep_max_n_entries(const unsigned int&); //set all entries to zero except n largest ones
  matrix<bool> flag_max_n_entries(const unsigned int&) const; //flag n largest entries with 1, set the rest to zero
  
  T min_entrie() const;
  matrix<T> min_n_entries(const unsigned int&) const;
  const matrix<T>& keep_min_n_entries(const unsigned int&); //set all entries to zero except n smallest ones
  matrix<bool> flag_min_n_entries(const unsigned int&) const; //flag n largest entries with 1, set the rest to zero

  const matrix<T>& remove_max_n_entries(const unsigned int&);
  const matrix<T>& remove_min_n_entries(const unsigned int&);

  unsigned int count_nonzero_entries() const;
  matrix<bool> flag_nonzero_entries() const; //flag n largest entries with 1, set the rest to zero
  
  void clear();
  
  void print(std::ostream&) const;
  void save(const std::string&) const;

  ~matrix();

  matrix<T> operator*(const matrix<T>&) const;
  
  template <typename Z> friend matrix<Z> abs(const matrix<Z>&);
  template <typename Z> friend double root_mean_square(const matrix<Z>&, const matrix<Z>&);
  template <typename Z> friend double Euclid_distance(const matrix<Z>&, const matrix<Z>&);
  template <typename Z> friend unsigned int Hamming_distance(const matrix<Z>&, const matrix<Z>&);
  template <typename Z> friend std::ostream& operator<<(std::ostream&, const matrix<Z>&);
};

#endif
