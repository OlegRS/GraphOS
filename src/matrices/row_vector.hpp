#ifndef __ROW_VECTOR_HPP__
#define __ROW_VECTOR_HPP__

#include "matrix.hpp"

template <typename T> class row_vector : public matrix<T> {
public:
  row_vector();
  row_vector(const unsigned int &n);
  row_vector(const row_vector&);
  row_vector(const matrix<T>&);

  unsigned int size() const;

  T& operator[](const unsigned int&) const;
};
#endif
