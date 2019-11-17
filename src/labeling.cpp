#include "labeling.hpp"

#define QUIET_MODE  // disables printing of some common warnings and info
// #define SILENT_MODE // disables printing of most warnings and info
#ifdef SILENT_MODE
#define QUIET_MODE
#endif

labeling::labeling() : size(0), radix(0), array(NULL), class_counters(NULL) {}
labeling::labeling(const unsigned int& array_size, const unsigned short& rdx, const unsigned short* arr) : size(array_size), radix(rdx) {
  array = new unsigned short[size];
  class_counters = new unsigned int[radix];
  if(!arr) {
  for(unsigned int i=0; i<size; ++i)
    array[i] = 0;
  class_counters[0] = size;
  for(unsigned short i=1; i<radix; ++i)
    class_counters[i] = 0;
  } 
  else {
    for(unsigned short i=0; i<radix; ++i)
      class_counters[i] = 0;
    for(unsigned int i=0; i<size; ++i)
      if((array[i] = arr[i]) < radix)
        ++class_counters[array[i]];
      else {
        std::cerr << "\n------------------------------------------------\n"
                  << "ERROR: Array is inconsitent with radix!\n"
                  << "------------------------------------------------\n";
        exit(1);
      }
    // if(check_empty_classes()) {
    //   std::cerr << "\n-----------------------------------------------------\n"
    //             << "ERROR: Attempt to create labeling with empty classes!\n"
    //             << "-----------------------------------------------------\n";
    //   exit(1);
    // }
  }
}
labeling::labeling(const labeling& L) : size(L.size), radix(L.radix) {
  array = new unsigned short[size];
  for(unsigned int i=0; i<size; ++i)
    array[i] = L.array[i];
  
  class_counters = new unsigned int[radix];
  for(unsigned int i=0; i<radix; ++i)
    class_counters[i] = L.class_counters[i];
}

const unsigned short& labeling::operator[](const unsigned int &i) const {
  return array[i];
}

labeling& labeling::set_label_by_index(const unsigned int& i, const unsigned short& label) {
  --class_counters[array[i]];
  array[i] = label;
  ++class_counters[label];
  return *this;
}

bool labeling::check_empty_classes() {
  for(unsigned short i=0; i<radix; ++i)
    if(!class_counters[i])
      return true;
  return false;
}

bool labeling::operator++() {
  do {
    if(size<1) {
      std::cerr << "\n------------------------------------------------\n"
                << "ERROR: Attempt to increment empty labeling!\n"
                << "------------------------------------------------\n";
      exit(1);
    }

    unsigned int i=size;
    while(i > 0 && array[i-1]==radix-1) {
      array[i-1] = 0;
      ++class_counters[0];
      --class_counters[radix-1];
      --i;
    }
    if(i > 0) {
      --class_counters[array[i-1]++];
      ++class_counters[array[i-1]];
    }
    
    else {
#ifndef QUIET_MODE
      std::cerr << "\n------------------------------------------------\n"
                << "WARNING: labeling overflow!\n"
                << "------------------------------------------------\n";
#endif
      return true;
    } 
  } while(check_empty_classes());
    
  return false;
}

bool labeling::operator--() {
  do {
    if(size<1) {
      std::cerr << "\n------------------------------------------------\n"
                << "ERROR: Attempt to decrement empty labeling!\n"
                << "------------------------------------------------\n";
      exit(1);
    }

    unsigned int i=size;
    while(i > 0 && array[i-1]==0) {
      array[i-1] = radix-1;
      --class_counters[0];
      ++class_counters[radix-1];
      --i;
    }

    if(i > 0) {
      --class_counters[array[i-1]--];
      ++class_counters[array[i-1]];
    }
    else {
#ifndef QUIET_MODE
      std::cerr << "\n------------------------------------------------\n"
                << "WARNING: labeling overflow (below zero)!\n"
                << "------------------------------------------------\n";
#endif
      return true;
    }
  } while(check_empty_classes());
  return false;
}

std::ostream& operator<<(std::ostream& os, const labeling& num) {
  for(unsigned int i=0; i<num.size; ++i)
    os << num.array[i] << ' ';
  return os;
}

std::istream& operator>>(std::istream& is, labeling& num) {
  for(unsigned int i=0; i<num.size; ++i)
    is >> num.array[i];
  return is;
}

labeling::~labeling() {
  delete[] array;
  delete[] class_counters;
}
