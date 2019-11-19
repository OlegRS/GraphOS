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
