#ifndef __NODE_HPP__
#define __NODE_HPP__

#include "graph.hpp"

class node {
  friend class link;
  friend class graph;
  friend class node_pair;
protected:
  
  std::string name;
  graph::label label;
  
  //Variables below only have meaning in a graph
  std::list<std::list<link>::iterator> attached_links;
  unsigned int id;
public:
  node() {}
  node(const std::string &name_) : name(name_) {}
  node(const std::string &name_, const std::string &label_) :  name(name_), label(label_) {}
  node(const unsigned int &id_, const std::string &name_, const unsigned int &label_id, const std::string &label_name) : name(name_), label(label_id,label_name), id(id_)  {}
  node(const node &nd) : name(nd.name), label(nd.label), id(nd.id) {}
  
  void set_name(const std::string &name_) {name = name_;}
  void set_label_name(const std::string &name_) {label.name = name_;}
  void set_label_id(const unsigned int &id_) {label.id = id_;}

  const std::string& get_name() const {return name;}
  const graph::label& get_label() const {return label;}
  const unsigned int& get_id() const {return id;}
  
  unsigned int degree() const { return attached_links.size(); }
  std::list<node*> adjacent_nodes_list() const;
  col_vector<node*> adjacent_nodes_col_vector() const;
  bool check_adjacency_to(const node& ) const;
  bool check_adjacency_to(const node* ) const;
  double clustering_coefficient() const;

  bool operator==(const node&);
  bool operator==(const std::string&);

  node& operator=(const node&);

  //////////////// FOR SPIN MODELS ///////////////
  void set_spin_to_minus1();
  void set_spin_to_plus1();
  //////////////// FOR SPIN MODELS ///////////////
  
  friend std::ostream& operator<<(std::ostream&, const node&);
  friend std::ostream& operator<<(std::ostream&, const link&);
  
  ~node() {};
};

#endif
