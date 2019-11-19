#ifndef __SYMM_MATRIX_HPP__
#define __SYMM_MATRIX_HPP__

#include "matrix.hpp"

#include <iostream>
#include <string>

class A_matrix;

template <typename T> class symm_matrix : public matrix<T> {//Symmetric matrix
public:
  symm_matrix();
  symm_matrix(const char* file_name);
  symm_matrix(const unsigned int &n);
  symm_matrix(const symm_matrix&);
  symm_matrix(const A_matrix&);
  symm_matrix(const matrix<T>&);

  unsigned int size() const;
  bool check_consistency() const;

  symm_matrix<T> sqr() const;
  symm_matrix<T> quad() const;
  symm_matrix<T> six() const;

  T* operator[](const unsigned int &i) const {return &matrix<T>::array[matrix<T>::dim_x*i];};
};

#endif
