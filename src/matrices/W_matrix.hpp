#ifndef __W_MATRIX_HPP__
#define __W_MATRIX_HPP__

#include "matrix.hpp"
#include <string>

class W_matrix : public matrix<double> {
public:
  W_matrix();
  W_matrix(const char*); //Loading from file
  W_matrix(const unsigned int&);
  W_matrix(const W_matrix&);

  double* operator[](const unsigned int&) const;
  const W_matrix& operator=(const W_matrix&);

  unsigned int size() const;
  bool check_consistency() const;

  void save(const std::string&) const;

  //////////// COMPUTATIONAL FUNCTION BEGIN ///////////
  double assortativity();
  ///////////// COMPUTATIONAL FUNCTION END ////////////

  friend std::ostream& operator<<(std::ostream&, const W_matrix&);
};

#endif
