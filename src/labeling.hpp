#ifndef __LABELING_HPP__
#define __LABELING_HPP__

#include <iostream>

class labeling {
private:
  unsigned int size;
  unsigned short radix;
  unsigned short* array;
  unsigned int* class_counters; //Counts nodes in geven class
public:
  labeling();
  labeling(const unsigned int& array_size, const unsigned short& radix, const unsigned short* array=NULL);
  labeling(const labeling&);

  bool check_empty_classes(); //Checks that there is at least one node in each class

  unsigned int num_nodes() const {return size;}
  unsigned short num_classes() const {return radix;}

  const unsigned short& operator[](const unsigned int&) const;
  labeling& set_label_by_index(const unsigned int& index, const unsigned short& label);
  bool operator++(); //returns 1 if labeling cannot be incremented
  bool operator--(); //returns 1 if labeling cannot be decremented
  bool operator==(const labeling&) const;
  bool operator>(const labeling&) const;
  bool operator<(const labeling&) const;
  labeling& operator=(labeling&);

  friend std::ostream& operator<<(std::ostream&, const labeling&);
  friend std::istream& operator>>(std::istream&, labeling&);
  
  ~labeling();
};

#endif
