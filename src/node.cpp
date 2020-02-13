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

std::list<node*> node::adjacent_nodes_list() const {
  std::list<node*> adjacent_nodes;
  for(std::list<std::list<link>::iterator>::const_iterator it_it_links = attached_links.begin(); it_it_links!=attached_links.end(); ++it_it_links) {
    node *p_node1 = (*it_it_links)->get_node1();
    p_node1 != this ? adjacent_nodes.push_back(p_node1) : adjacent_nodes.push_back((*it_it_links)->get_node2());
  }
  return adjacent_nodes;
}

col_vector<node*> node::adjacent_nodes_col_vector() const {
  col_vector<node*> adjacent_nodes(degree());
  unsigned int i = 0;
  for(std::list<std::list<link>::iterator>::const_iterator it_it_links = attached_links.begin(); it_it_links!=attached_links.end(); ++it_it_links) {
    node *p_node1 = (*it_it_links)->get_node1(); 
    p_node1 != this ? adjacent_nodes[i]=p_node1 : adjacent_nodes[i]=(*it_it_links)->get_node2();
    ++i;
  }
  return adjacent_nodes;
}

bool node::check_adjacency_to(const node& nd) const {
  for(std::list<std::list<link>::iterator>::const_iterator it_it_links = attached_links.begin(); it_it_links!=attached_links.end(); ++it_it_links) {
    node *p_node1 = (*it_it_links)->get_node1();
    if(p_node1 != this) {
      if(p_node1 == &nd)
        return true;
    }
    else
      if((*it_it_links)->get_node2() == &nd)
        return true;
  }
  return false;
}

bool node::check_adjacency_to(const node* p_node) const {
  for(std::list<std::list<link>::iterator>::const_iterator it_it_links = attached_links.begin(); it_it_links!=attached_links.end(); ++it_it_links) {
    node *p_node1 = (*it_it_links)->get_node1();
    if(p_node1 != this) {
      if(p_node1 == p_node)
        return true;
    }
    else
      if((*it_it_links)->get_node2() == p_node)
        return true;
  }
  return false;
}

double node::clustering_coefficient() const {
  col_vector<node*> adjacent_nodes = adjacent_nodes_col_vector();
  unsigned int k = adjacent_nodes.size();
  if(k<2) return 0;
  double cc = 0;
  for(unsigned int i=0; i<k; ++i)
    for(unsigned int j=0; j<i; ++j)
      cc += adjacent_nodes[i]->check_adjacency_to(adjacent_nodes[j]);
  
  return 2*cc/(k*(k-1));
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
