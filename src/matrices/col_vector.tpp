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

template <typename T> unsigned int col_vector<T>::size() const {
  return matrix<T>::dim_y;
}

template <typename T> T& col_vector<T>::operator[](const unsigned int &i) const {
  return matrix<T>::array[i];
}

#endif
