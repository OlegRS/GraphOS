#include "graph.hpp"

graph::label& graph::label::operator=(const graph::label &l) {
  id = l.id;
  name = l.name;
  return *this;
}

graph::label& graph::label::operator=(const std::string &str) {
  name = str;
  return *this;
}

bool graph::label::operator==(const std::string &str) {
  if(name == str)
    return true;
  return false;
}

bool graph::label::operator==(const graph::label &l) {
  if(name == l.name)
    return true;
  return false;
}

std::ostream& operator<<(std::ostream &os, const graph::label &l) {
  os << l.name;
  return os;
}

bool graph::label::operator<(const graph::label &l) {
  return name < l.name;
}
