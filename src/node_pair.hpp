#ifndef __NODE_PAIR_HPP__
#define __NODE_PAIR_HPP__

#include "node.hpp"
#include <string>

class node_pair {
protected:
  node *node1 = NULL;
  node *node2 = NULL;
  graph *gr = NULL; // Graph to which the pair belongs 
public:
  node_pair() : node1(NULL), node2(NULL), gr(NULL) {}
  node_pair(node *p_node1, node *p_node2, graph *p_gr) : node1(p_node1), node2(p_node2), gr(p_gr) {}
  node* get_node1() const {return node1;}
  node* get_node2() const {return node2;}

  bool linked() const;
  node_pair& link_nodes(const std::string& type="pp", const double& weight=0);
  node_pair& link_nodes_no_checks(const std::string& type="pp", const double& weight=0);
  
  node_pair& unlink_nodes(const std::string& type="pp", const double& weight=0);

  node_pair& flip_link_state(); //Changes linked pair to unlinked and vice versa.
  
  node_pair& remove_nodes();

  bool operator==(node_pair&);
  bool operator!=(node_pair&);

  friend std::ostream& operator<<(std::ostream&, const node_pair&);
};

#endif
