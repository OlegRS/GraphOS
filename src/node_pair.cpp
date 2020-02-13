#include "node_pair.hpp"

bool node_pair::linked() const {return node1->check_adjacency_to(node2);}

node_pair& node_pair::link_nodes(const std::string& type, const double& weight) {
  gr->add_link_no_checks(node1, node2, type, weight);
  return *this;
}

node_pair& node_pair::link_nodes_no_checks(const std::string& type, const double& weight) {
  gr->add_link_no_checks(node1, node2, type, weight);
  return *this;
}

node_pair& node_pair::unlink_nodes(const std::string& type, const double& weight) {
  gr->remove_link(node1, node2, type, weight);
  return *this;
};

node_pair& node_pair::flip_link_state() {
  if(linked())
    unlink_nodes();
  else
    link_nodes();
  return *this;
}

node_pair& node_pair::remove_nodes() {
  gr -> remove_node_no_checks(node1);
  gr -> remove_node_no_checks(node2);
  node1 = node2 = NULL;
  return *this;
}

bool node_pair::operator==(node_pair& np) {
  if((node1 == np.node1 && node2 == np.node2) || (node1 == np.node2 && node2 == np.node1))
    return true;
  else
    return false;
}

bool node_pair::operator!=(node_pair& np) {
  if((node1 != np.node1 || node2 != np.node2) && (node1 != np.node2 || node2 != np.node1))
    return true;
  else
    return false;
}
 
std::ostream& operator<<(std::ostream &os, const node_pair &np) {
  os << "Pair_of_nodes_belongs_to_graph_with_pointer " << np.gr << ":\n"
     << "node1: " << *np.node1 << ";\n"
     << "node2: " << *np.node2 << '.';
  return os;
}
