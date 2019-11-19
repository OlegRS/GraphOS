#include "node.hpp"

void node::set_spin_to_minus1() {
  if(label.name != "-1") {
    label.name = "-1";
    label.id = !label.id;
  }
}

void node::set_spin_to_plus1() {
 if(label.name != "1") {
    label.name = "1";
    label.id = !label.id;
  } 
}

bool node::operator==(const std::string &name_) {
  if(name != name_)
    return false;
  
  return true;
}

bool node::operator==(const node &nd) {
  if(name != nd.name)
    return false;
  
  return true;
}

node& node::operator=(const node &nd) {
  name = nd.name;
  label = nd.label;
  id = nd.id;
  attached_links = nd.attached_links;

  return *this;
}

std::ostream& operator<<(std::ostream &os, const node &nd) {
  os << nd.name << ' ' << nd.label.name << " degree = " << nd.attached_links.size();
  return os;
}
