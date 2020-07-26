#include "A_matrix.hpp"
#include "../aux_math.hpp"
#include "col_vector.tpp"
#include "matrix.tpp"
#include "symm_matrix.tpp"

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

unsigned int A_matrix::degree(const unsigned int& i) const {
  unsigned int k=0;
  for(unsigned j=0; j<dim_x; ++j)
    k+=array[dim_x*i + j];
  return k;
}

col_vector<unsigned int> A_matrix::degree_sequence_col_vec() const {
  col_vector<unsigned int> deg_seq(dim_x);
  for(unsigned int i=0; i<dim_x; ++i)
    deg_seq[i] = degree(i);
  return deg_seq;
}

long double A_matrix::num_p_stars(const unsigned int& i, const unsigned int& p) const {
  return aux_math::binom(degree(i), p);
}

col_vector<long double> A_matrix::num_p_stars_sequence_col_vec(const unsigned int& p) const {
  col_vector<long double> p_stars_sequence(dim_x);
  for(unsigned int i=0; i<dim_x; ++i)
    p_stars_sequence[i] = num_p_stars(i, p);
  return p_stars_sequence;
}

long double A_matrix::average_p_stars(const unsigned int& p) const {
  long double s = 0;
  for(unsigned int i=0; i<dim_x; ++i)
    s+=num_p_stars(i, p);
  return s/dim_x;
}

unsigned long long A_matrix::num_triangles(const unsigned int& n) const {
  // Finding adjacent nodes
  unsigned int k = degree(n);
  if(k<2) return 0;
  
  col_vector<unsigned int> neighbours(k); // Adjacent nodes
  unsigned int j=0;
  for(unsigned int i=0; i<dim_x && j<k; ++i)
    if((*this)[i][n])
      neighbours[j++]=i;
  
  // Counting connected adjacent nodes
  unsigned long long Nt = 0; // Number of triangles of the node
  for(j=0; j<k; ++j)
    for(unsigned int l=0; l<j; ++l)
      Nt += (*this)[neighbours[j]][neighbours[l]];
  return Nt;
}

unsigned long long A_matrix::num_triangles() const {
  symm_matrix<uint64_t> S = symm_matrix<uint64_t>(*this);
  return (S.sqr()*S).tr()/6.; //INEFFICIENT!!!
}

double A_matrix::loc_clust_coeff(const unsigned int& n) const {
  long double num_2s = num_p_stars(n,2);
  if(num_2s > 0)
    return num_triangles(n)/num_2s;
  else
    return 0;
}

col_vector<double> A_matrix::loc_clust_coeff_col_vec() const {
  col_vector<double> lcc(dim_x);
  for(unsigned int i=0; i<dim_x; ++i)
    lcc[i] = loc_clust_coeff(i);
  return lcc;
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

A_matrix& A_matrix::set_Erdos_Renyi(const double &p, prng &rnd) {
  unsigned int rand_max = rnd.rand_max();
  for(unsigned int i=0; i<dim_x; ++i)
    for(unsigned int j=0; j<i; ++j)
      rnd() < p*rand_max ? array[dim_x*i+j]=array[dim_x*j+i]=true : array[dim_x*i+j]=array[dim_x*j+i]=false;

  for(unsigned int i=0; i<dim_x; ++i)
    array[dim_x * i + i] = false;
  
  return *this;
}

col_vector<col_vector<unsigned int> > A_matrix::random_node_pairs_col_vector(const unsigned int &N_pairs, prng &rnd) const {
  if(dim_x*(dim_x-1)/2 < N_pairs) {
    std::cerr << "------------------------------------------------------------------------------------------\n"
              << "ERROR: Attempt to obtain too many distinct pairs of nodes from a graph which is too small!\n"
              << "------------------------------------------------------------------------------------------\n";
    exit(1);
  }
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

A_matrix& A_matrix::single_link_GB_Metropolis_generator(double(&H)(const A_matrix&), const unsigned int &N_iters, prng &rnd, const bool &initialize_randomly, const double &temp) {
  //// Initialiasing adjacency matrix ////
  unsigned int rand_max = rnd.rand_max();
  A_matrix A_new(dim_x);
  if(initialize_randomly) {
    // Initialising Erdos-Renyi with random p
    double p = rnd();
    for(unsigned int i=0; i<dim_x; ++i)
      for(unsigned int j=0; j<i; ++j)
        if(rnd() < p)
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
      i = rnd()%dim_x;
      j = rnd()%dim_x;
    } while(i==j);
    A_new[i][j] = A_new[j][i] = !(*this)[i][j];
    if(rnd() < exp(1/temp*(H(*this)-H(A_new)))*rand_max)
      (*this)[i][j] = (*this)[j][i] = A_new[i][j];
    else
      A_new[i][j] = A_new[j][i] = (*this)[i][j];
  }
 
  return *this;
}

A_matrix& A_matrix::single_link_MF_GB_Metropolis_generator(double(&H)(const unsigned int& N_links), const unsigned int &N_iters, prng &rnd, const bool &initialize_randomly, const double &temp) {//Metropolis dynamics optimized for Hamiltonians that depend only on the number of links
  unsigned int rand_max = rnd.rand_max();
  //// Initialiasing adjacency matrix ////
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
      if(rnd()/(double)rand_max < exp(1/temp*(H(L)-H(L+1)))) {
        (*this)[i][j] = (*this)[j][i] = A_new[i][j];
        L++;
      }
      else
        A_new[i][j] = A_new[j][i] = (*this)[i][j];
    else
      if(rnd()/(double)rand_max < exp(1/temp*(H(L)-H(L-1)))) {
        (*this)[i][j] = (*this)[j][i] = A_new[i][j];
        L--;
      }
      else
        A_new[i][j] = A_new[j][i] = (*this)[i][j];

  }
 
  return *this;
}

A_matrix& A_matrix::MF_GB_Metropolis_generator(double(&H)(const unsigned int& N_links), const unsigned int &N_iters, prng &rnd, const unsigned int &N_pairs_max, const bool &initialize_randomly, const double &temp) {//Metropolis dynamics optimized for Hamiltonians that depend only on the number of links
  unsigned int rand_max = rnd.rand_max();
  //// Initialiasing adjacency matrix ////
  if(initialize_randomly)
    set_Erdos_Renyi(rnd()/(double)rand_max, rnd);

  //// Running Metropolis algorithm ////
  unsigned int L = num_links();
  unsigned int L_new = L;
  for(unsigned int n=0; n<N_iters; ++n) {
    unsigned int N_pairs = rnd()%N_pairs_max+1;
    col_vector< col_vector<unsigned int> > pairs = random_node_pairs_col_vector(N_pairs, rnd);
    for(unsigned int l=0; l<N_pairs; ++l)
      (*this)[pairs[l][0]][pairs[l][1]] ? --L_new : ++L_new;

    if(rnd() < exp(1/temp*(H(L)-H(L_new)))*rand_max)
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

A_matrix& A_matrix::sample_p_star_model(const unsigned int &N_iters, prng& rnd, const col_vector<double> &T, const unsigned int &N_pairs_max, const bool &initialize_randomly, const double &temp) {
  unsigned int rand_max = rnd.rand_max();
  if(initialize_randomly)
    set_Erdos_Renyi(rnd()/rand_max, rnd);
  
  unsigned int p = T.size(); //p-star model
  col_vector<double> T_rescaled = T; //Rescaled parameters
  for(unsigned int s=0; s<p; ++s)
    T_rescaled[s] = T[s]/pow(dim_x, s); //(s+1)! is already accounted for in the number of stars

  // Running Metropolis dynamics
  for(unsigned int n=0; n<N_iters; ++n) {
    unsigned int N_pairs = rnd()%N_pairs_max + 1;
    col_vector<col_vector<unsigned int>  > np = random_node_pairs_col_vector(N_pairs, rnd);
    double delta_H = 0;
    for(unsigned int i=0; i<N_pairs; ++i) {
      unsigned int k1=degree(np[i][0]);
      unsigned int k2=degree(np[i][1]);
      if((*this)[np[i][0]][np[i][1]]) { // If there is a link, remove it
        (*this)[np[i][0]][np[i][1]] = (*this)[np[i][1]][np[i][0]] = 0;
        for(unsigned int s=1; s<=p; ++s)
          delta_H += T_rescaled[s-1]*s*(aux_math::binom(k1,s)/k1 + aux_math::binom(k2,s)/k2); // C_n^k - C_(n-1)^k = k/n*C_n^k
      }
      else { // If there is no link, add it
        (*this)[np[i][0]][np[i][1]] = (*this)[np[i][1]][np[i][0]] = 1;
        for(unsigned int s=1; s<=p; ++s)
          delta_H -= T_rescaled[s-1]*s*(aux_math::binom(k1+1,s)/(k1+1) + aux_math::binom(k2+1,s)/(k2+1));
      }   
    }

    if(rnd() > exp(-1/temp*delta_H)*rand_max)
      //Reject the proposal and return everything to how it was by flipping link states again
      for(unsigned int i=0; i<N_pairs; ++i) {
        if((*this)[np[i][0]][np[i][1]]) // If there is a link, remove it
          (*this)[np[i][0]][np[i][1]] = (*this)[np[i][1]][np[i][0]] = 0;
        else
          (*this)[np[i][0]][np[i][1]] = (*this)[np[i][1]][np[i][0]] = 1;
      }
  }
  return *this;
}

A_matrix& A_matrix::sample_p_star_model_with_single_link_Metropolis(const unsigned int &N_iters, prng& rnd, const col_vector<double> &T, const bool &initialize_randomly, const double &temp) {
  unsigned int rand_max = rnd.rand_max();
  
  if(initialize_randomly)
    set_Erdos_Renyi(rnd()/rand_max, rnd);
  
  // Obtaining initial degree sequences
  col_vector<unsigned int> k = degree_sequence_col_vec();
  unsigned int p = T.size();

  col_vector<double> T_rescaled = T; //Rescaled parameters
  for(unsigned int s=0; s<p; ++s)
    T_rescaled[s] = T[s]/pow(dim_x, s); //(s+1)! is already accounted for in the number of stars

  // Running Metropolis dynamics
  unsigned int i,j;
  for(unsigned int n=0; n<N_iters; ++n) {
    do {
      i = rnd()%dim_x;
      j = rnd()%dim_x;
    } while(i==j);
    
    double delta_H=0;
    if((*this)[i][j]) { //If there is a link, propose to remove it
      for(unsigned int s=1; s<=p; ++s)
        delta_H += T_rescaled[s-1]*s*(aux_math::binom(k[i],s)/k[i] + aux_math::binom(k[j],s)/k[j]); // C_n^k - C_(n-1)^k = k/n*C_n^k

      if(rnd() < exp(-1/temp*delta_H)*rand_max) {
        (*this)[i][j] = (*this)[j][i] = 0;
        --k[i]; --k[j];
      }
    } 
    else { //If there is no link, propose to add it
      for(unsigned int s=1; s<=p; ++s)
        delta_H -= T_rescaled[s-1]*s*(aux_math::binom(k[i]+1,s)/(k[i]+1) + aux_math::binom(k[j]+1,s)/(k[j]+1)); // C_n^k - C_(n-1)^k = k/n*C_n^k

      if(rnd() < exp(-1/temp*delta_H)*rand_max) {
        (*this)[i][j] = (*this)[j][i] = 1;
        ++k[i]; ++k[j];
      }
    }
  }
  return *this;
}

A_matrix& A_matrix::sample_triad_model_with_single_link_Metropolis(const unsigned int &N_iters, prng& rnd, const double h, const double sigma, const double tau, const bool &initialize_randomly, const double &temp) {
  unsigned int rand_max = rnd.rand_max();
  
  if(initialize_randomly)
    set_Erdos_Renyi(rnd()/rand_max, rnd);

  double h_rescaled = 2*h;
  double sigma_rescaled = 2*sigma/dim_x;
  double tau_rescaled = 6*tau/dim_x;
  
  // Obtaining initial degree sequences
  col_vector<unsigned int> k = degree_sequence_col_vec();

  // Running Metropolis dynamics
  unsigned int i,j;
  for(unsigned int n=0; n<N_iters; ++n) {
    do {
      i = rnd()%dim_x;
      j = rnd()%dim_x;
    } while(i==j);
    
    double delta_H=0;
    if((*this)[i][j]) { // If there is a link, propose to remove it
      delta_H += h_rescaled; // Field contribution 
      delta_H += sigma_rescaled*2*(aux_math::binom(k[i],2)/k[i] + aux_math::binom(k[j],2)/k[j]); // 2-star contribution (C_n^k - C_(n-1)^k = k/n*C_n^k)
      // Triangles contribution
      unsigned long long N_t = 0; // Number of triangles containing an edge (i,j)
      for(unsigned int k=0; k<dim_x; ++k)
        if((*this)[i][k] && (*this)[j][k])
          ++N_t;
      delta_H += tau_rescaled * N_t;

      if(rnd() < exp(-1/temp*delta_H)*rand_max) {
        (*this)[i][j] = (*this)[j][i] = 0;
        --k[i]; --k[j];
      }
    } 
    else { //If there is no link, propose to add it
        delta_H -= h_rescaled; // Field contribution 
        delta_H -= sigma_rescaled*2*(aux_math::binom(k[i]+1,2)/(k[i]+1) + aux_math::binom(k[j]+1,2)/(k[j]+1)); // 2-star contribution (C_n^k - C_(n-1)^k = k/n*C_n^k)
        // Triangles contribution
        unsigned long long N_t = 0; // Number of triangles containing an edge (i,j)
        for(unsigned int k=0; k<dim_x; ++k)
          if((*this)[i][k] && (*this)[j][k])
            ++N_t;
        delta_H -= tau_rescaled * N_t;

      if(rnd() < exp(-1/temp*delta_H)*rand_max) {
        (*this)[i][j] = (*this)[j][i] = 1;
        ++k[i]; ++k[j];
      }
    }
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
