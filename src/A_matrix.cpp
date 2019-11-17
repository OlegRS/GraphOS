#include "matrix.hpp"

A_matrix::A_matrix() : symm_matrix() {}
A_matrix::A_matrix(const char* file_name) : symm_matrix<bool>(file_name) {
  //Check that matrix only has entries 0 or 1
  for(unsigned int i=0; i<dim_y*dim_x; i++)
    if(matrix<bool>::array[i]!=0 && array[i]!=1) {
      std::cerr << "----------------------------------------------------------------------------\n";
      std::cerr << "ERROR loading from file: Adjacency matrix should have entries either 0 or 1\n";
      std::cerr << "----------------------------------------------------------------------------\n";
      exit(1);
    }
}

A_matrix::A_matrix(const unsigned int &sz) : symm_matrix<bool>(sz) {
  unsigned int size = dim_x*dim_y;
  for(unsigned int i=0; i<size; ++i)
    array[i] = false;
}

A_matrix::A_matrix(const A_matrix& A)  {
  matrix<bool>::dim_x = dim_y = A.dim_x;
  unsigned int sz = dim_x*dim_y;
  
  array = new bool[sz];
  
  for(unsigned int i=0; i<sz; i++)
    array[i] = A.array[i];
}

const A_matrix& A_matrix::operator=(const A_matrix &A) {
  dim_x = dim_y = A.dim_x;
  unsigned int sz = dim_x*dim_y;

  delete[] array;  
  array = new bool[sz];

  for(unsigned int i=0; i<sz; i++)
      array[i] = A.array[i];

  return *this;
}

bool A_matrix::check_consistency() const {
  //Check that matrix is square
  if(dim_x != dim_y) {
    std::cerr << "------------------------------------------------\n";
    std::cerr << "ERROR: Adjacency matrix should be square matrix\n";
    std::cerr << "------------------------------------------------\n";
    return false;
  }
  
  //Check that matrix is symmetric with all entries either 1 or 0
  for(unsigned int i=0; i<dim_y; i++)
    for(unsigned int j=0; j<dim_x; j++) {
      if(array[i*dim_x+j] != array[j*dim_x+i]) {
        std::cerr << "---------------------------------------------------\n";
        std::cerr << "ERROR: Adjacency matrix should be symmetric matrix\n";
        std::cerr << "---------------------------------------------------\n";
        return false;
      }
      if(i==j and array[i*dim_x+j]!=0) {
        std::cerr << "-------------------------------------------------------------------------\n";
        std::cerr << "ERROR: Diagonal elements of adjacency matrix of simple graph should be 0\n";
        std::cerr << "-------------------------------------------------------------------------\n";
        return false;
      }
    }
  return true;
}

unsigned int A_matrix::num_links() const {
  unsigned int num_links = 0;
  for(unsigned int i=0; i<dim_x; ++i)
    for(unsigned int j=0; j<i; ++j)
      num_links+=(*this)[i][j];
  return num_links;
}

void A_matrix::save(const std::string& file_name) const {
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

std::ostream& operator<<(std::ostream& os, const A_matrix& A) {
  for(unsigned int i=0; i<A.dim_y; i++) {
    for(unsigned int j=0; j < A.dim_x-1; j++)
      os << A.array[i*A.dim_x + j] << ' ';
    os << A.array[i*A.dim_x + A.dim_x-1] << '\n';
  }
  return os;
};

A_matrix& A_matrix::GB_Metropolis_generator(double(&H)(const A_matrix&), const unsigned int &N_iters, const int &seed, const bool &initialize_randomly, const double &temp) {
  srand(seed);
  //// Initialiasing adjacency matrix ////
  A_matrix A_new(dim_x);
  if(initialize_randomly) {
    for(unsigned int i=0; i<dim_x; ++i)
      for(unsigned int j=0; j<i; ++j) {
        (*this)[i][j] = (*this)[j][i] = rand()%2;
      }
    for(unsigned int i=0; i<dim_x; ++i)
      (*this)[i][i] = 0;
  }
  A_new = (*this);
  //// Running Metropolis algorithm ////
  unsigned int i,j;
  for(unsigned int n=0; n<N_iters; ++n) {
    do {
      i = rand()%dim_x;
      j = rand()%dim_x;
    } while(i==j);
    A_new[i][j] = A_new[j][i] = !(*this)[i][j];
    if(rand()/(double)RAND_MAX < exp(1/temp*(H(*this)-H(A_new))))
      (*this)[i][j] = (*this)[j][i] = A_new[i][j];
    else
      A_new[i][j] = A_new[j][i] = (*this)[i][j];
  }
 
  return *this;
}

A_matrix& A_matrix::MF_GB_Metropolis_generator(double(&H)(const unsigned int& N_links), const unsigned int &N_iters, const int &seed, const bool &initialize_randomly, const double &temp) {//Metropolis dynamics optimized for Hamiltonians that depend only on the number of links
  //// Initialiasing adjacency matrix ////
  A_matrix A_new(dim_x);
  unsigned int L = 0; //Number of links
  if(initialize_randomly) {
    for(unsigned int i=0; i<dim_x; ++i)
      for(unsigned int j=0; j<i; ++j) {
        L += (*this)[i][j] = (*this)[j][i] = rand()%2;
      }
    for(unsigned int i=0; i<dim_x; ++i)
      (*this)[i][i] = 0;
  }
  else //Computing number of links in the original adjacency matrix
    for(unsigned int i=0; i<dim_x; ++i)
      for(unsigned int j=0; j<i; ++j) {
        L += (*this)[i][j];
      }
  A_new = (*this);

  //// Running Metropolis algorithm ////
  unsigned int i,j;
  for(unsigned int n=0; n<N_iters; ++n) {
    do {
      i = rand()%dim_x;
      j = rand()%dim_x;
    } while(i==j);
    A_new[i][j] = A_new[j][i] = !(*this)[i][j];
    if((*this)[i][j] == 0)
      if(rand()/(double)RAND_MAX < exp(1/temp*(H(L)-H(L+1)))) {
        (*this)[i][j] = (*this)[j][i] = A_new[i][j];
        L++;
      }
      else
        A_new[i][j] = A_new[j][i] = (*this)[i][j];
    else
      if(rand()/(double)RAND_MAX < exp(1/temp*(H(L)-H(L-1)))) {
        (*this)[i][j] = (*this)[j][i] = A_new[i][j];
        L--;
      }
      else
        A_new[i][j] = A_new[j][i] = (*this)[i][j];

  }
 
  return *this;
}
