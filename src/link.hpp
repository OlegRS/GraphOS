#ifndef __LINK_HPP__
#define __LINK_HPP__

#include "graph.hpp"

class link {
  friend class graph;
protected:
  node* node1;
  node* node2;
  std::string type;
  double weight;
public:
  link(const link &l) : node1(l.node1), node2(l.node2), type(l.type), weight(l.weight) {}
  link(node* p_nd1=NULL, node* p_nd2=NULL, const std::string &type_="pp", double weight_=0) : node1(p_nd1), node2(p_nd2), type(type_), weight(weight_) {}
  // link(node&, node&, const std::string& ="pp", double=0);

  void set_weight(const double &weight_) {weight = weight_;}
  void set_type(const std::string &type_) {type = type_;}
  void set_node1(node *p_node1) {node1 = p_node1;}
  void set_node2(node *p_node2) {node2 = p_node2;}

  double get_weight() const {return weight;}
  const std::string& get_type() const {return type;}
  node* get_node1() const {return node1;}
  node* get_node2() const {return node2;}
  
  link& operator=(const link&);

  bool operator==(const link&);
  bool operator!=(const link&);
  
  ~link() {};

  friend std::ostream& operator<<(std::ostream&, const link&);
};

#endif
