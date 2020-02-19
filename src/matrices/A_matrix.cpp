#include "A_matrix.hpp"
#include "col_vector.tpp"
#include "matrix.tpp"
#include "symm_matrix.tpp"
#include <math.h>
#include <fstream>
#include <functional>

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

unsigned int A_matrix::num_links() const {
  unsigned int N_links = 0;
  for(unsigned int i=0; i<dim_x; ++i)
    for(unsigned int j=0; j<i; ++j)
      N_links += array[i*dim_x+j];
  return N_links;
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

A_matrix& A_matrix::set_Erdos_Renyi(const double &p, const unsigned int &seed) {
  std::default_random_engine generator;
  generator.seed(seed);
  std::uniform_int_distribution<> distribution(0, RAND_MAX);
  auto rnd = std::bind(distribution, generator);
  
  for(unsigned int i=0; i<dim_x; ++i)
    for(unsigned int j=0; j<i; ++j)
      rnd() < p*RAND_MAX ? array[dim_x*i+j]=array[dim_x*j+i]=true : array[dim_x * i + j]=array[dim_x*j+i]=false;

  for(unsigned int i=0; i<dim_x; ++i)
    array[dim_x * i + i] = false;
  
  return *this;
}

col_vector<col_vector<unsigned int> > A_matrix::random_node_pairs_col_vector(const unsigned int &N_pairs, const int &seed) const {
  if(dim_x*(dim_x-1)/2 < N_pairs) {
    std::cerr << "------------------------------------------------------------------------------------------\n"
              << "ERROR: Attempt to obtain too many distinct pairs of nodes from a graph which is too small!\n"
              << "------------------------------------------------------------------------------------------\n";
    exit(1);
  }

  std::default_random_engine generator;
  generator.seed(seed);
  std::uniform_int_distribution<> distribution(0, RAND_MAX);
  auto rnd = std::bind(distribution, generator);

  col_vector<col_vector<unsigned int> > pairs(N_pairs);
  unsigned int i, j;
  for(unsigned int l=0; l<N_pairs; ++l) {
    do {
      i = rnd()%dim_x;
      j = rnd()%dim_x;
    } while(i==j);
    pairs[l] = col_vector<unsigned int>(2);
    pairs[l][0] = i; pairs[l][1] = j;
    for(unsigned int k=0; k<l; ++k)
      if((pairs[k][0] == pairs[l][0] && pairs[k][1] == pairs[l][1]) || (pairs[k][0] == pairs[l][1] && pairs[k][1] == pairs[l][0]))
        --l;
  }
  
  return pairs;
}

A_matrix& A_matrix::GB_Metropolis_generator(double(&H)(const A_matrix&), const unsigned int &N_iters, const int &seed, const bool &initialize_randomly, const double &temp) {
  srand(seed);
  //// Initialiasing adjacency matrix ////
  A_matrix A_new(dim_x);
  if(initialize_randomly) {
    // Initialising Erdos-Renyi with random p
    double p = rand();
    for(unsigned int i=0; i<dim_x; ++i)
      for(unsigned int j=0; j<i; ++j)
        if(rand() < p)
          (*this)[i][j] = (*this)[j][i] = 1;
        else
          (*this)[i][j] = (*this)[j][i] = 0;
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
    if(rand() < exp(1/temp*(H(*this)-H(A_new)))*RAND_MAX)
      (*this)[i][j] = (*this)[j][i] = A_new[i][j];
    else
      A_new[i][j] = A_new[j][i] = (*this)[i][j];
  }
 
  return *this;
}

A_matrix& A_matrix::single_link_MF_GB_Metropolis_generator(double(&H)(const unsigned int& N_links), const unsigned int &N_iters, const int &seed, const bool &initialize_randomly, const double &temp) {//Metropolis dynamics optimized for Hamiltonians that depend only on the number of links
  //// Initialiasing adjacency matrix ////
  std::default_random_engine generator;
  generator.seed(seed);
  std::uniform_int_distribution<> distribution(0, RAND_MAX);
  auto rnd = std::bind(distribution, generator);

  A_matrix A_new(dim_x);
  unsigned int L = 0; //Number of links
  if(initialize_randomly) {
    // Initialising Erdos-Renyi with random p
    double p = rnd();
    for(unsigned int i=0; i<dim_x; ++i)
      for(unsigned int j=0; j<i; ++j)
        if(rnd() < p)
          L += (*this)[i][j] = (*this)[j][i] = 1;
        else
          (*this)[i][j] = (*this)[j][i] = 0;
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
      i = rnd()%dim_x;
      j = rnd()%dim_x;
    } while(i==j);
    A_new[i][j] = A_new[j][i] = !(*this)[i][j];
    if((*this)[i][j] == 0)
      if(rnd()/(double)RAND_MAX < exp(1/temp*(H(L)-H(L+1)))) {
        (*this)[i][j] = (*this)[j][i] = A_new[i][j];
        L++;
      }
      else
        A_new[i][j] = A_new[j][i] = (*this)[i][j];
    else
      if(rnd()/(double)RAND_MAX < exp(1/temp*(H(L)-H(L-1)))) {
        (*this)[i][j] = (*this)[j][i] = A_new[i][j];
        L--;
      }
      else
        A_new[i][j] = A_new[j][i] = (*this)[i][j];

  }
 
  return *this;
}

A_matrix& A_matrix::MF_GB_Metropolis_generator(double(&H)(const unsigned int& N_links), const unsigned int &N_iters, const int &seed, const unsigned int &N_pairs_max, const bool &initialize_randomly, const double &temp) {//Metropolis dynamics optimized for Hamiltonians that depend only on the number of links
  //// Initialiasing adjacency matrix ////
  std::default_random_engine generator;
  generator.seed(seed);
  std::uniform_int_distribution<> distribution(0, RAND_MAX);
  auto rnd = std::bind(distribution, generator);

  if(initialize_randomly)
    set_Erdos_Renyi(rnd()/(double)RAND_MAX, rnd());

  //// Running Metropolis algorithm ////
  unsigned int L = num_links();
  unsigned int L_new = L;
  for(unsigned int n=0; n<N_iters; ++n) {
    unsigned int N_pairs = rnd()%N_pairs_max+1;
    col_vector< col_vector<unsigned int> > pairs = random_node_pairs_col_vector(N_pairs, rnd());
    for(unsigned int l=0; l<N_pairs; ++l)
      (*this)[pairs[l][0]][pairs[l][1]] ? --L_new : ++L_new;

    if(rnd() < exp(1/temp*(H(L)-H(L_new)))*RAND_MAX)
       // Accept the proposal and flip the pairs
      for(unsigned int l=0; l<N_pairs; ++l) {
        if((*this)[pairs[l][0]][pairs[l][1]])
          (*this)[pairs[l][0]][pairs[l][1]]=(*this)[pairs[l][1]][pairs[l][0]]=false;
        else
          (*this)[pairs[l][0]][pairs[l][1]]=(*this)[pairs[l][1]][pairs[l][0]]=true;
        L = L_new;
      }
    else
      L_new  = L;
  }
         
  return *this;
}

std::ostream& operator<<(std::ostream& os, const A_matrix& A) {
  for(unsigned int i=0; i<A.dim_y; i++) {
    for(unsigned int j=0; j < A.dim_x-1; j++)
      os << A.array[i*A.dim_x + j] << ' ';
    os << A.array[i*A.dim_x + A.dim_x-1] << '\n';
  }
  return os;
};
