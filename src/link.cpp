#include "link.hpp"

link& link::operator=(const link &l) {
  node1 = l.node1;
  node2 = l.node2;
  type = l.type;
  weight = l.weight;
  
  return *this;
}

bool link::operator==(const link &l) {
  if((node1!=l.node1 && node1!=l.node2) || (node2!=l.node2 && node2!=l.node1) || weight!=l.weight || type != l.type)
    return false;
  
  return true;
}

bool link::operator!=(const link &l) {
  return !(*this == l);
}

std::ostream& operator<<(std::ostream &os, const link &l) {
  os << l.node1->name << ' ' << l.type << ' ' << l.node2->name;
  return os;
}
