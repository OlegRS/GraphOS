#ifndef MATRIX_H
#define MATRIX_H

#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <iostream>
#include <fstream>

template <typename T> class matrix {
protected:
  unsigned int dim_x;
  unsigned int dim_y;

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

  template <typename Z> friend matrix<Z> abs(const matrix<Z>&);
  template <typename Z> friend double root_mean_square(const matrix<Z>&, const matrix<Z>&);
  template <typename Z> friend double Euclid_distance(const matrix<Z>&, const matrix<Z>&);
  template <typename Z> friend unsigned int Hamming_distance(const matrix<Z>&, const matrix<Z>&);
  template <typename Z> friend std::ostream& operator<<(std::ostream&, const matrix<Z>&);
};

template <typename T> class col_vector : public matrix<T> {
public:
  col_vector();
  col_vector(const unsigned int &n);
  col_vector(const col_vector&);
  col_vector(const matrix<T>&);

  unsigned int size() const;

  T& operator[](const unsigned int&) const;
};

template <typename T> class row_vector : public matrix<T> {
public:
  row_vector();
  row_vector(const unsigned int &n);
  row_vector(const row_vector&);
  row_vector(const matrix<T>&);

  unsigned int size() const;

  T& operator[](const unsigned int&) const;
};

class W_matrix : public matrix<double> {
public:
  W_matrix();
  W_matrix(const char*);
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

/////////////////////////////////////////////////////////////////////////////////
//////////////////////// MATRIX DEFINITIONS /////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////

template <typename T> matrix<T>::matrix() : dim_x(0), dim_y(0), array(NULL) {}
template <typename T> matrix<T>::matrix(const matrix<T> &M) : dim_x(M.dim_x), dim_y(M.dim_y) {
  unsigned int sz = dim_x * dim_y;
  array = new T[sz];
  for(unsigned int i=0; i<sz; i++)
    array[i] = M.array[i];
}
template <typename T> matrix<T>::matrix(const unsigned int &dim_y_, const unsigned int &dim_x_)
  : dim_x(dim_x_), dim_y(dim_y_) {
  unsigned int sz = dim_x * dim_y;
  array = new T[sz];
}
template <typename T> matrix<T>::matrix(const unsigned int &dim_x_) : matrix(dim_x_, dim_x_) {}
template <typename T> matrix<T>::matrix(const std::string &file_name) {  

  std::ifstream ifs(file_name.c_str());
  
  if( !ifs.is_open() ) {
    std::cerr << "-------------------------------------------------------------------\n"
              << "ERROR: Cannot open the file " << file_name << " to load the matrix!\n"
              << "-------------------------------------------------------------------\n";
    exit(1);
  }

  dim_x = dim_y = 0;
  for( char a = ifs.get(); !ifs.eof(); a = ifs.get() )
    if(a == ' ')
      dim_x++;
    else
      if( a == '\n' )
        dim_y++;

  dim_x += dim_y;
  if(dim_x % dim_y) {
    std::cerr << "--------------------------------------------------\n"
              << "ERROR: Matrix " << file_name << " has uneven rows!\n"
              << "--------------------------------------------------\n";
    exit(2);
  }
  
  ifs.clear();
  ifs.seekg(0);

  array = new T[dim_x];
  dim_x = dim_x/dim_y;

  //dim_x--; //This can be used to load matrices with '\n' in the end of each row
  for(unsigned int i =0; i<dim_y; i++)
    for(unsigned int j=0; j<dim_x; j++) {
      ifs >> array[i*dim_x + j];
    }
  
  ifs.close();
}

template <typename T> matrix<T>::~matrix() {  delete[] array; }

template <typename T> T* matrix<T>::operator[](const unsigned int& i) const {
  return &array[dim_x * i]; //Returns pointer to an element
}

template <typename T> matrix<T>& matrix<T>::operator+=(const matrix<T> &M) {
  if(dim_x == M.dim_x && dim_y == M.dim_y) {
    unsigned int sz = dim_x * dim_y;
    for(unsigned int i=0; i<sz; i++)
      array[i] += M.array[i];
    return *this;
  }
  else {
    std::cerr << "-------------------------------------------------------------------------\n"
              << "ERROR: Addition assignment is called on matrix with different dimensions!\n"
              << "-------------------------------------------------------------------------\n";
    return *this;
  }
}

template <typename T> matrix<T>& matrix<T>::operator/=(double x) {
  if(x) {
    unsigned int sz = dim_x * dim_y;
    for(unsigned int i=0; i<sz; i++)
      array[i] /= x;
    return *this;
  }
  else {
    std::cerr << "----------------------------------------------\n"
              << "ERROR: Elementwise division of matrix by zero!\n"
              << "----------------------------------------------\n";
    exit(1);
    //return *this;
  }
}

template <typename T> const matrix<T>& matrix<T>::operator=(const matrix<T> &M) {
  array = new T[M.dim_x * M.dim_y];
  dim_x = M.dim_x;
  dim_y = M.dim_y;
  unsigned int sz = dim_x * dim_y;
  for(unsigned int i=0; i<sz; i++)
    array[i] = M.array[i];
  return *this;
}

template <typename T> const matrix<T>& matrix<T>::abs() {
  unsigned int sz = dim_x*dim_y;
  for(unsigned int i=0; i<sz; ++i)
    if(array[i]<0) array[i] *= -1;
  return *this;
}

template <typename T> T matrix<T>::max_entrie() const {
  T max = array[0];
  unsigned int sz = dim_x*dim_y;
  for(unsigned int i=1; i<sz; ++i)
    if(array[i] > max)
      max = array[i];
  return max;
}

template <typename T> matrix<T> matrix<T>::max_n_entries(const unsigned int &n) const {
  unsigned int sz = dim_x*dim_y;
  if(n <= sz) {
    matrix<T> M(n,1);
    matrix<T> A(*this);
    T min = min_entrie();
    unsigned int l;

    for(unsigned int j=0; j<n; ++j) {
      l=0;
      M[j][0] = A.array[0];
      for(unsigned int i=1; i<sz; ++i)
        if(A.array[i] > M[j][0]) {
          M[j][0] = A.array[i];
          l = i;
        }
      A.array[l] = min;
    }

    return M;
  }
  else {
    std::cerr << "---------------------------------------------------------------------------------------\n"
              << "ERROR: Attempt to find max "<< n <<" entries of matrix with only " << sz << " elements!\n"
              << "---------------------------------------------------------------------------------------\n";
    return matrix<T>();
  }
}

template <typename T> const matrix<T>& matrix<T>::keep_max_n_entries(const unsigned int &n) {
unsigned int sz = dim_x*dim_y;
  if(n <= sz) {
    matrix<T> A(*this);
    T min = min_entrie();
    
    for(unsigned int i=0; i<sz; ++i)
      if(array[i] != min)
        array[i] = 0;
    
    T max;
    unsigned int l, j;
    for(j=0; j<n; ++j) {
      l=0;
      max = A.array[0];
      for(unsigned int i=1; i<sz; ++i)
        if(A.array[i] > max) {
          max = A.array[i];
          l = i;
        }
      if(max == min)  break;
      A.array[l] = min;
      array[l] = max;
    }

    l=0;
    for(; n-j; l++)  
      if(array[l] == min) 
        j++;
    
    for(; l<sz; l++)
      if(array[l] == min)
        array[l] = 0;
  }
  else {
    std::cerr << "-----------------------------------------------------------------------------------------------\n"
              << "ERROR: Attempt to keep the largest "<< n <<" entries of matrix with only " << sz << " elements!\n"
              << "-----------------------------------------------------------------------------------------------\n";
  }
  return *this;
}

template <typename T> matrix<bool> matrix<T>::flag_max_n_entries(const unsigned int &n) const {
  matrix<T> M(*this);
  M.keep_max_n_entries(n);
  matrix<bool> B(M.dim_y, M.dim_x);
  for(unsigned int i=0; i<dim_y; i++)
    for(unsigned int j=0; j<dim_x; j++)
      if(M[i][j]) B[i][j] = true;
      else B[i][j] = false;

  return B;
}

template <typename T> T matrix<T>::min_entrie() const {
  T min = array[0];
  unsigned int sz = dim_x*dim_y;
  for(unsigned int i=1; i<sz; ++i)
    if(array[i] < min)
      min = array[i];
  return min;
}

template <typename T> matrix<T> matrix<T>::min_n_entries(const unsigned int &n) const {
  unsigned int sz = dim_x*dim_y;
  if(n <= sz) {
    matrix<T> M(n,1);
    matrix<T> A(*this);
    T max = max_entrie();
    unsigned int l;

    for(unsigned int j=0; j<n; ++j) {
      l=0;
      M[j][0] = A.array[0];
      for(unsigned int i=1; i<sz; ++i)
        if(A.array[i] < M[j][0]) {
          M[j][0] = A.array[i];
          l = i;
        }
      A.array[l] = max;
    }

    return M;
  }
  else {
    std::cerr << "---------------------------------------------------------------------------------------\n"
              << "ERROR: Attempt to find min "<< n <<" entries of matrix with only " << sz << " elements!\n"
              << "---------------------------------------------------------------------------------------\n";
    return matrix<T>();
  }
}

template <typename T> const matrix<T>& matrix<T>::keep_min_n_entries(const unsigned int &n) {
unsigned int sz = dim_x*dim_y;
  if(n <= sz) {
    matrix<T> A(*this);
    T max = max_entrie();
    
    for(unsigned int i=0; i<sz; ++i)
      if(array[i] != max)
        array[i] = 0;
    
    T min;
    unsigned int l, j;
    for(j=0; j<n; ++j) {
      l=0;
      min = A.array[0];
      for(unsigned int i=1; i<sz; ++i)
        if(A.array[i] < min) {
          min = A.array[i];
          l = i;
        }
      if(min == max)  break;
      A.array[l] = max;
      array[l] = min;
    }

    l=0;
    for(; n-j; l++)  
      if(array[l] == max) 
        j++;
    
    for(; l<sz; l++)
      if(array[l] == max)
        array[l] = 0;
    
  }
  else {
    std::cerr << "------------------------------------------------------------------------------------------------\n"
              << "ERROR: Attempt to keep the smallest "<< n <<" entries of matrix with only " << sz << " elements!\n"
              << "------------------------------------------------------------------------------------------------\n";
  }
  return *this;
}

template <typename T> matrix<bool> matrix<T>::flag_min_n_entries(const unsigned int &n) const {
  matrix<T> M(*this);
  M.keep_min_n_entries(n);
  matrix<bool> B(M.dim_y, M.dim_x);
  for(unsigned int i=0; i<dim_y; i++)
    for(unsigned int j=0; j<dim_x; j++)
      if(M[i][j]) B[i][j] = true;
      else B[i][j] = false;

  return B;
}

template <typename T> const matrix<T>& matrix<T>::remove_max_n_entries(const unsigned int &n) {
  matrix<bool> B = flag_max_n_entries(n);
  for(unsigned int i=0; i<dim_y; ++i)
    for(unsigned int j=0; j<dim_x; ++j)
      if(B[i][j])
        array[i*dim_x+j] = 0;
  return *this;
}

template <typename T> const matrix<T>& matrix<T>::remove_min_n_entries(const unsigned int &n) {
  matrix<bool> B = flag_min_n_entries(n);
  for(unsigned int i=0; i<dim_y; ++i)
    for(unsigned int j=0; j<dim_x; ++j)
      if(B[i][j])
        array[i*dim_x+j] = 0;
  return *this;
}


template <typename T> unsigned int matrix<T>::count_nonzero_entries() const {
  unsigned int count = 0;
  unsigned int sz = dim_x * dim_y;
  for(unsigned int i=0; i<sz; ++i)
    if(array[i])  count++;
  
  return count;
}

template <typename T> matrix<bool> matrix<T>::flag_nonzero_entries() const {
  matrix<bool> M(dim_y,dim_x);
  for(unsigned int i=0; i<dim_y; i++)
    for(unsigned int j=0; j<dim_y; j++)
      if((*this)[i][j])
        M[i][j] = true;
      else
        M[i][j] = false;
  return M;
}

template <typename T> void matrix<T>::clear() {
  delete[] array;
  array = NULL;
  dim_x = dim_y = 0;
}

template <typename T> void matrix<T>::print(std::ostream& os) const {
  for(unsigned int i=0; i<dim_y; i++) {
    for(unsigned int j=0; j < dim_x-1; j++)
      os << array[i*dim_x + j] << ' ';
    os << array[i*dim_x + dim_x-1] << '\n';
  }
}

template <typename T> void matrix<T>::save(const std::string &file_name) const {
  std::ofstream ofs(file_name.c_str());
  if(!ofs.is_open()) {
    std::cerr << "-----------------------------------------------------------------\n"
              << "ERROR: Cannot open the file " << file_name << " to save the matrix!\n"
              << "-----------------------------------------------------------------\n";
  }
  this -> print(ofs);
  ofs.close();
}

template <typename T>  matrix<T> abs(const matrix<T> &M) {
  matrix<T> abs_M(M);
  return abs_M.abs();
}

template <typename T> double root_mean_square(const matrix<T> &M1, const matrix<T> &M2) {
  unsigned int sz = M1.dim_x * M1.dim_y;
  
  if(sz == M2.dim_x*M2.dim_y) {
    double d = 0;
    for(unsigned int i=0; i<sz; i++)
      d += (M1.array[i] - M2.array[i])*(M1.array[i] - M2.array[i]);
    return sqrt(d/double(sz));
  }
  else {
    std::cerr << "--------------------------------------------------------------------------------------------\n"
              << "ERROR: Attempt to compute Euclidean distance of matrices with defferent numbers of elements!\n"
              << "--------------------------------------------------------------------------------------------\n";
    return 0;
  }
}
template <typename T> double Euclid_distance(const matrix<T> &M1, const matrix<T> &M2) {
  unsigned int sz = M1.dim_x * M1.dim_y;
  
  if(sz == M2.dim_x*M2.dim_y) {
    double d = 0;
    for(unsigned int i=0; i<sz; i++)
      d += (M1.array[i] - M2.array[i])*(M1.array[i] - M2.array[i]);
    return sqrt(d);
  }
  else {
    std::cerr << "--------------------------------------------------------------------------------------------\n"
              << "ERROR: Attempt to compute Euclidean distance of matrices with defferent numbers of elements!\n"
              << "--------------------------------------------------------------------------------------------\n";
    return 0;
  }
}

template <typename T> unsigned int Hamming_distance(const matrix<T> &M1, const matrix<T> &M2) {
  unsigned int sz = M1.dim_x * M1.dim_y;
  
  if(sz == M2.dim_x*M2.dim_y) {
    unsigned int count = 0;
    for(unsigned int i=0; i<sz; i++)
      if(M1.array[i] != M2.array[i])
        ++count;
    return count;
  }
  else {
    std::cerr << "------------------------------------------------------------------------------------------\n"
              << "ERROR: Attempt to compute Hamming distance of matrices with defferent numbers of elements!\n"
              << "------------------------------------------------------------------------------------------\n";
    return 0;
  }
}

template <typename T> std::ostream& operator<<(std::ostream &os, const matrix<T> &M) {
  M.print(os);
  return os;
}

/////////////////////////////////////////////////////////////////////////////////
//////////////////////// COLUMN VECTOR DEFINITIONS //////////////////////////////
/////////////////////////////////////////////////////////////////////////////////
template <typename T> col_vector<T>::col_vector() : matrix<T>() {}
template <typename T> col_vector<T>::col_vector(const unsigned int &n) : matrix<T>(n,1) {}
template <typename T> col_vector<T>::col_vector(const col_vector& v) : matrix<T>(v) {}
template <typename T> col_vector<T>::col_vector(const matrix<T> &M) : matrix<T>(M) {
  matrix<T>::dim_y = matrix<T>::dim_x*matrix<T>::dim_y;
  matrix<T>::dim_x = 1;
}

template <typename T> unsigned int col_vector<T>::size() const {
  return matrix<T>::dim_y;
}

template <typename T> T& col_vector<T>::operator[](const unsigned int &i) const {
  return matrix<T>::array[i];
}

/////////////////////////////////////////////////////////////////////////////////
///////////////////////// ROW VECTOR DEFINITIONS ////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////
template <typename T> row_vector<T>::row_vector() : matrix<T>() {}
template <typename T> row_vector<T>::row_vector(const unsigned int &n) : matrix<T>(1,n) {}
template <typename T> row_vector<T>::row_vector(const row_vector &v) : matrix<T>(v) {}
template <typename T> row_vector<T>::row_vector(const matrix<T> &M) : matrix<T>(M) {
  matrix<T>::dim_x = matrix<T>::dim_x*matrix<T>::dim_y;
  matrix<T>::dim_y = 1;
}

template <typename T> unsigned int row_vector<T>::size() const {
  return matrix<T>::dim_x;
}

template <typename T> T& row_vector<T>::operator[](const unsigned int &i) const { return matrix<T>::array[i]; }

#endif
