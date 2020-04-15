#ifndef __COL_VECTOR_HPP__
#define __COL_VECTOR_HPP__

#include "matrix.hpp"
#include <list>

template <typename T> class col_vector : public matrix<T> {//Column vector
public:
  col_vector();
  col_vector(const unsigned int& n);
  col_vector(const unsigned int& size, const T*);
  col_vector(const col_vector&);
  col_vector(const matrix<T>&);
  col_vector(const std::list<T>&);

  unsigned int size() const;

  T avrg() const;

  std::list<T> unique_elements() const;

  T& operator[](const unsigned int&) const;
};

#endif
