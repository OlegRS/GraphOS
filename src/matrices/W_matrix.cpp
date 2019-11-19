#include "W_matrix.hpp"
#include "col_vector.hpp"
#include <iostream>
#include <fstream>
#include <math.h>

#include "col_vector.tpp"
#include "matrix.tpp"


W_matrix::W_matrix() : matrix() {}
W_matrix::W_matrix(const char* file_name) : matrix(file_name) {
  //Check that matrix is square
  if(dim_x != dim_y) {
    std::cerr << "---------------------------------------\n";
    std::cerr << "ERROR: W(x,x') should be square matrix!\n";
    std::cerr << "---------------------------------------\n";
    exit(1);
  }
  //Check that matrix is symmetric and normalised
  double sum = 0;
  for(unsigned int i=0; i<dim_y; i++)
    for(unsigned int j=0; j<i; j++)
      if(array[i*dim_x+j] == array[j*dim_x+i])
        sum += 2*array[i*dim_x+j];
      else {
        std::cerr << "------------------------------------------\n";
        std::cerr << "ERROR: W(x,x') should be symmetric matrix!\n";
        std::cerr << "------------------------------------------\n";
        exit(2);
      }
  for(unsigned int i=0; i<dim_y; i++)
    sum += array[i*(dim_x+1)];
  if(sum-1 > 0.00001) {
    std::cerr << "----------------------------------------------------------------------------\n";
    std::cerr << "WARNING: Sum of all entries of W(x,x') is " << sum <<" instead of one!\n";
    std::cerr << "----------------------------------------------------------------------------\n";
  }
}

W_matrix::W_matrix(const unsigned int &sz) : matrix<double>(sz) {
  unsigned int size = dim_x*dim_y;
  for(unsigned int i=0; i<size; ++i)
    array[i] = 0;
}

W_matrix::W_matrix(const W_matrix& W)  {
  
  dim_x = dim_y = W.dim_x;
  unsigned int sz = dim_x*dim_y;

  array = new double[sz];
  
  for(unsigned int i=0; i<sz; i++)
    array[i] = W.array[i];
}

unsigned int W_matrix::size() const {
  return dim_x;
}

double* W_matrix::operator[](const unsigned int &i) const {
  return &array[dim_x * i];
}

const W_matrix& W_matrix::operator=(const W_matrix &W) {

  dim_x = dim_y = W.dim_x;
  unsigned int sz = dim_x*dim_y;

  delete[] array;  
  array = new double[sz];
  
  for(unsigned int i=0; i<sz; i++)
      array[i] = W.array[i];
  
  return *this;
}

bool W_matrix::check_consistency() const {
  double sum = 0;
  for(unsigned int i=0; i < dim_y*dim_x; i++)
    sum += array[i];
  
  if(sum-1 > 0.00001 || 1-sum > 0.00001) {
    std::cerr << "----------------------------------------------------------------------------\n";
    std::cerr << "WARNING: Sum of all entries of W_matrix is " << sum <<" instead of one!\n";
    std::cerr << "----------------------------------------------------------------------------\n";
    return false;
  }
  return true;
}

void W_matrix::save(const std::string& file_name) const {
  std::ofstream ofs(file_name.c_str());
  if(!ofs.is_open()) {
    std::cerr << "-------------------------------------------------------------------\n";
    std::cerr << "ERROR: Cannot open the file " << file_name << " to save the matrix!\n";
    std::cerr << "-------------------------------------------------------------------\n";
    return;
  }
  ofs << *this;
  ofs.close();
}

std::ostream& operator<<(std::ostream& os, const W_matrix& W) {
  for(unsigned int i=0; i<W.dim_y; i++) {
    for(unsigned int j=0; j < W.dim_x-1; j++)
      os << W.array[i*W.dim_x + j] << ' ';
    os << W.array[i*W.dim_x + W.dim_x-1] << '\n';
  }
  return os;
};

double W_matrix::assortativity() {
  col_vector<double> marg(dim_x);
  for(unsigned int i=1; i<dim_x; ++i) {
    marg[i] = 0;
    for(unsigned int j=1; j<dim_x; ++j)
      marg[i] += (*this)[i][j];
  }
  
  double sum1=0;
  for(unsigned int i=1; i<dim_x; ++i) {
    for(unsigned int j=1; j<i; ++j)
      sum1 += 2*(*this)[i][j]*i*j;
    sum1 += (*this)[i][i]*i*i;
  }

  double sum2 = 0;
  for(unsigned int i=1; i<dim_x; ++i)
    sum2 += marg[i]*i;
  sum2 *= sum2;

  double sum3 = 0;
  for(unsigned int i=1; i<dim_x; ++i)
    sum3 += marg[i]*i*i;

  if(sum3 - sum2)
    return (sum1 - sum2) / (sum3 - sum2);
  else
    if(sum1 - sum2 == 0) return 1;
    else {
      std::cerr << "------------------------------------------------------\n";
      std::cerr << "ERROR: Cannot compute assortativity! Division by zero!\n";
      std::cerr << "------------------------------------------------------\n";
      return 1000;
    }
}
