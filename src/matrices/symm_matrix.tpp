#ifndef __SYMM_MATRIX_TPP__
#define __SYMM_MATRIX_TPP__

template <typename T> symm_matrix<T>::symm_matrix() : matrix<T>() {}
template <typename T> symm_matrix<T>::symm_matrix(const char* file_name) : matrix<T>(file_name) {
  //Check that matrix is square
  if(matrix<T>::dim_x != matrix<T>::dim_y) {
    std::cerr << "------------------------------------------------\n";
    std::cerr << "ERROR: Symmetric matrix should be square \n";
    std::cerr << "------------------------------------------------\n";
    exit(1);
  }
  //Check that matrix is symmetric
  for(unsigned int i=0; i<matrix<T>::dim_y; i++)
    for(unsigned int j=0; j<i; j++)
      if(matrix<T>::array[i*matrix<T>::dim_x+j] != matrix<T>::array[j*matrix<T>::dim_x+i]) {
        std::cerr << "---------------------------------------------------\n";
        std::cerr << "ERROR: Symmetric matrix should be symmetric \n";
        std::cerr << "---------------------------------------------------\n";
        exit(2);
      }
}
template <typename T> symm_matrix<T>::symm_matrix(const unsigned int &n) : matrix<T>(n,n) {}
template <typename T> symm_matrix<T>::symm_matrix(const symm_matrix &M) : matrix<T>(M) {}
template <typename T> symm_matrix<T>::symm_matrix(const A_matrix &A) //: matrix<T>::dim_x(A.dim_x), matrix<T>::dim_y(A.dim_y) {
{
  matrix<T>::dim_x = A.get_dim_x();
  matrix<T>::dim_y = A.get_dim_y();
  unsigned int sz = matrix<T>::dim_x * matrix<T>::dim_y;
  matrix<T>::array = new T[sz];
  for(unsigned int i=0; i<matrix<T>::dim_y; ++i)
    for(unsigned int j=0; j<=i; ++j)
      matrix<T>::array[i*matrix<T>::dim_x+j] = matrix<T>::array[j*matrix<T>::dim_x+i] = A[i][j];
}
template <typename T> symm_matrix<T>::symm_matrix(const matrix<T> &M) : matrix<T>(M) {
  if(matrix<T>::dim_y != matrix<T>::dim_x) {
    std::cerr << "------------------------------------------------------\n"
              << "ERROR: Attempt to copy asymmetric matrix to symmetric\n"
              << "------------------------------------------------------\n";
    exit(1);
  }
}
template <typename T> unsigned int symm_matrix<T>::size() const {
  return matrix<T>::dim_x;
}

template <typename T> bool symm_matrix<T>::check_consistency() const {
  if(matrix<T>::dim_x != matrix<T>::dim_y)
    return false;
  for(int i=0; i<matrix<T>::dim_y; ++i)
    for(int j=i+1; j<matrix<T>::dim_x; ++j)
      if(matrix<T>::array[i*matrix<T>::dim_x+j] != matrix<T>::array[j*matrix<T>::dim_x+i]) {
        std::cout << "i=" << i << ", j=" << j << '\n';
        return false;
      }
  return true;
}

template <typename T> symm_matrix<T> symm_matrix<T>::sqr() const {
  symm_matrix<T> M_(matrix<T>::dim_x);
  for(int i=0; i< matrix<T>::dim_y; ++i)
    for(int j=i; j< matrix<T>::dim_x; ++j) {
      M_.array[i*matrix<T>::dim_x+j] = M_.array[j*matrix<T>::dim_x+i] = 0;
      for(int l=0; l< matrix<T>::dim_x; ++l)
        M_.array[i*matrix<T>::dim_x+j] = M_.array[j*matrix<T>::dim_x+i]+=matrix<T>::array[i*matrix<T>::dim_x+l]*matrix<T>::array[l*matrix<T>::dim_x+j];
    }
  
  return M_;
}

template <typename T> symm_matrix<T> symm_matrix<T>::quad() const {return (this->sqr()).sqr();}

template <typename T> symm_matrix<T> symm_matrix<T>::six() const {return (this->quad())*(this->sqr());}

#endif
