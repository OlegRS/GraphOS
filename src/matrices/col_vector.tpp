#ifndef __COL_VECTOR_TPP__
#define __COL_VECTOR_TPP__

template <typename T> col_vector<T>::col_vector() : matrix<T>() {}
template <typename T> col_vector<T>::col_vector(const unsigned int &n) : matrix<T>(n,1) {}
template <typename T> col_vector<T>::col_vector(const unsigned int& size, const T* array) : matrix<T>(size, 1) {
  for(unsigned int i=0; i<size; ++i)
    matrix<T>::array[i] = array[i];
}
template <typename T> col_vector<T>::col_vector(const col_vector& v) : matrix<T>(v) {}
template <typename T> col_vector<T>::col_vector(const matrix<T> &M) : matrix<T>(M) {
  matrix<T>::dim_y = matrix<T>::dim_x*matrix<T>::dim_y;
  matrix<T>::dim_x = 1;
}
template <typename T> col_vector<T>::col_vector(const std::list<T>& L) : matrix<T>(L.size(), 1) {
  unsigned int i=0;
  for(typename std::list<T>::const_iterator it_L=L.begin(); it_L!=L.end(); ++it_L) {
    matrix<T>::array[i] = *it_L;
    ++i;
  }
}

template <typename T> unsigned int col_vector<T>::size() const {
  return matrix<T>::dim_y;
}

template <typename T> T col_vector<T>::avrg() const {
  T avrg = 0;
  for(unsigned int i=0; i<matrix<T>::dim_y; ++i)
    avrg += matrix<T>::array[i];
  return avrg/(double)matrix<T>::dim_y;
}

template <typename T> std::list<T> col_vector<T>::unique_elements() const {
  std::list<T> ue;
  for(unsigned int i=0; i<matrix<T>::dim_y; ++i) {
    typename std::list<T>::iterator it_ue = ue.begin();
    while(it_ue!=ue.end()) {
      if(*it_ue==(*this)[i])
        break;
      ++it_ue;
    }
    if(it_ue==ue.end())
      ue.push_back((*this)[i]);
  }
  return ue;
}

template <typename T> T& col_vector<T>::operator[](const unsigned int &i) const {
  return matrix<T>::array[i];
}

#endif
