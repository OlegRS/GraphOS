//TO DO
  // 1). Search by name over nodes should be improved by sorting the array of nodes or implementing good hash table!
  // 2). Since graph stores labels it will be better to make nodes store pointers to labels instead of labels themselves node::node.set_label() should be changed to something that checks what are the other labels in the graph, and storing pointers instead of labels themselves will make it natural. Also it will be more efficient to implement node::set_label(unsigned int label_id). But this optimization is questionable since nodes outside the graph will have no labels. In this case node should be a protected class in class graph and node should store the pointer to a graph to which it belongs. This can be implememted differently. Node can erase label_name everytime it is added to a grph, since label has id. Alternatively if we get rid of ids we can store labels as char* and everytime we add node to a graph delete[] this char and redirect the pointer to the corresponding label - this is better since node does not have to know anything about the graph in this case and char* replaces the id.
  // 3). Do the same thing with links. Create class graph::link_type : public graph::label and store link types in the graph while only storing pionters to the link types in the links which are in the graph.
  // 4). Implement your own lists (or, better, deque or que) and instead of saving array of itterators to links in each node, save array of pointers
  // 5). Possibly get rid of node and label ids
  // 6). Change link type type from std::string to char*, since empty std::string weighs 56 bytes.
  // 7). "remove_nodes_with_degree" and "remove_nodes_with_label" functions should be optimised
//TO DO
#ifndef GRAPHLIB_H
#define GRAPHLIB_H

#include <time.h>
#include <math.h>
#include <list>

#include "matrix.h"

#include <stdlib.h>
#include <stdio.h>
#include <iostream>
#include <fstream>

/////////////////////////////
class node;
class link;
class graph;
/////////////////////////////

class graph {
  
public:
  
  struct label; //nested class
  
protected:
  
  bool NODES_SORTED = false;
  unsigned int NODES_ARRAY_SIZE;
  
  node* nodes;
  label* labels;
  std::list<link> links;
  
  unsigned int N_nodes; //Number of nodes (all elements of nodes array after N_nodes are treated as inactive)
  unsigned int N_labels; //Number of different labels
  unsigned int N_links; //Number of links
  
  bool remove_link(std::list<link>::iterator&); //Removes link with certain iterator

  graph& remove_node_no_checks(node *p_node);

  label* find_label(const std::string&);
  
  void move_node(const unsigned int&, const unsigned int&); //Moves node in the array redirecting links (needed for sorting)
  void sort_nodes(); //Places all inactive nodes to the right and the rest sorts by name
  
public:
  
  graph();
  graph(const graph&); //Copy constructor
  graph(const unsigned int& N); //Create graph reserving memory for N nodes
  graph(const std::string&, unsigned int N_additional_nodes=0); //Load labeled graph from efficient format to array, allocating more memory than needed to add nodes later on.
  graph(const std::string&, const std::string&);  //Load graph from Cytoscape .sif file and from labels assignment .attrs file

  void append_nodes_array(const unsigned int &N); //Reallocates memory for array of nodes adding extra space for N nodes
  
  unsigned int num_labels() const {return N_labels;}
  unsigned int num_nodes() const {return N_nodes;}
  unsigned int num_links() const {return N_links;}

  void add_node(const node&, const std::string&); //May be useful for creating intersecting graphs with shared nodes
  graph& add_node(const std::string &name, const std::string &label_name);
  graph& add_node(const std::string &name) { return add_node(name,""); }
  graph& add_node_no_checks(const std::string &name, const std::string &label_name);
  
  graph& remove_node_no_checks(const std::string&); //removes node by name, removes its links, but does not check if any labels should be removed
  graph& remove_node(const std::string&); //remove by name and remove the links
  // graph& remove_node_random(); //NOT IMPLEMENTED
  // graph& remove_random_node_with_label(const std::string&); //NOT IMPLEMENTED
  graph& remove_nodes_with_degree(const unsigned int &k);
  graph& remove_nodes_with_label(const std::string&);

  void add_link_no_checks(const unsigned int&, const unsigned int&, std::string="pp", double weight=0);
  graph& add_link(const unsigned int&, const unsigned int&, std::string="pp", double weight=0);
  void add_link(const std::string&, const std::string&, const std::string& ="pp", double=0);
  void add_link(const std::string&, const std::string&, double&);
  void add_link(node*, node*, const std::string& ="pp", double=0);
    
  void remove_link(node*, node*, const std::string& ="pp", double=0);
  void remove_link(const std::string&, const std::string&, const std::string& ="pp", double=0);
  void remove_link(unsigned int&);

  node* find_node(const std::string &node_name) const;
  node* find_node(const unsigned int &node_id) const; //Returns a node with certain position in the array (id)

  const node& get_node(const std::string &node_name) const;
  const node& get_node(const unsigned int &node_id) const; //Returns a node with certain position in the array (id)

  inline link* get_link(const unsigned int &node1_id, const unsigned int &node2_id) const;

  bool is_simple() const;
  
  std::string* get_labels() const;
  matrix<std::string> get_labels_vec() const;

  void replicate_nodes(const unsigned int &replication_factor); //Proportionally creates many nodes with same labels, but different names

  matrix<bool> adjacency_matrix() const;

  void save_labels(const std::string&) const;
  
  void save(const std::string&) const;
  void save_to_sif_and_attrs(const std::string &name) const;

  graph& load(const std::string&, unsigned int N_additional_nodes=0); //Loads the graph with free space for additional nodes

  void clear_links();
  void clear();

  //////////////////// COMPUTATIONAL FUNCTIONS BEGIN ////////////////////////
  const unsigned int& number_of_nodes() const { return N_nodes;  }
  const unsigned int& number_of_links() const { return N_links;  }
  const unsigned int& number_of_labels() const { return N_labels; }
  
  unsigned int max_degree() const;
  unsigned int min_degree() const;
  unsigned int* degree_sequence() const; //Implicitly allocates memory!!!
  col_vector<unsigned int> degree_sequence_col_vec() const;

  col_vector<node*> nodes_with_max_degree_col_vec() const;
  col_vector<node*> nodes_with_min_degree_col_vec() const;
  col_vector<node*> nodes_with_degree_col_vec(unsigned int &k) const;

  std::string* label_sequence() const; //Implicitly allocates memory!!!
  col_vector<std::string> label_sequence_col_vector() const;
  
  double average_degree() const;
  
  double degree_distribution(const unsigned int &k) const;
  double* degree_distribution() const; //Implicitly allocates memory!!!
  col_vector<double> degree_distribution_col_vector() const;

  double label_distribution(std::string &label_name) const;
  double label_distribution(unsigned int &label_id) const;
  double* label_distribution() const; //Implicitly allocates memory!!!
  matrix<double> label_distribution_vec() const;
  
  double joint_distribution_of_connected_node_labels(const std::string &label_1_name, const std::string &label_2_name) const;
  double joint_distribution_of_connected_node_labels(const unsigned int &label_1_id, const unsigned int &label_2_id) const;
  W_matrix joint_distribution_of_connected_node_labels() const;
  
  double joint_distribution_of_connected_node_degrees(const unsigned int &k1, const unsigned int &k2) const;
  W_matrix joint_distribution_of_connected_node_degrees() const;

  double modularity() const;
  double relative_modularity() const; //Returns ratio between modularity and max modularity with current number of labels.

  double assortativity() const;

  //// Spin model functions begin ////
  const graph& randomize_spins(const int &seed);
  const graph& simulate_Glauber_Dynamics(const int &seed, const double &Temperature, const unsigned int &N_steps);
  
  short* spin_sequence() const; //Implicitly allocates memory!!!
  col_vector<short> spin_sequence_col_vector() const;
  matrix<double> interaction_matrix() const;  
  double average_magnetization() const;

  //Associative memory functions
  graph& set_spin_sequence(const col_vector<short>&);
  graph& set_spin_sequence(const row_vector<short>&);
  graph& set_spin_sequence_2D(const matrix<short>&);
  matrix<short> spin_sequence_2D(const unsigned int &dim_y, const unsigned int &dim_x) const; //Usefull to represent 2-D pictures to play with associative memory (with each node being a pixel)
  //// Spin model functions end ////

  //// Classification functions begin ////
  // This set of functions performs inference of labels in the maximum-entropy ensemble with soft constraints on W_labels.
  double loglikelihood() const; //Computes likelihood of the label assignment in the max entropy ensemble with soft constraints on W_labels.
  bool increment_label_assignment();
  graph& classify_nodes_precisely(const unsigned int &num_classes);
  graph& classify_nodes_greedy_heuristic(const unsigned int &num_classes);
  //// Classification functions end ////
  
  
  ///////////////////// COMPUTATIONAL FUNCTIONS END //////////////////////////
  
  ~graph();
  
  friend std::ostream& operator<<(std::ostream& ,const graph&);
};

/////////////////////////////////////////////////////////////////////////
struct graph::label {
  unsigned int id;
  std::string name;

  label();
  label(const std::string& , unsigned int =-1);
  label(const unsigned int&, const std::string&);
  label(const label&);

  label& operator=(const label&);
  label& operator=(const std::string&);

  bool operator==(const std::string&);
  bool operator==(const label&);
  bool operator<(const label&);
};

/////////////////////////////////////////////////////////////////////////
class node {
  friend class link;
  friend class graph;
protected:
  
  std::string name;
  graph::label label;
  
  //Variables below only have meaning in a graph
  std::list<std::list<link>::iterator> attached_links;
  unsigned int id;
public:
  node();
  node(const std::string&);
  node(const std::string&, const std::string&);
  node(const unsigned int&, const std::string&, const unsigned int&, const std::string&);
  node(const node&);
  
  void set_name(const std::string &name_) {name = name_;}
  void set_label_name(const std::string &name_) {label.name = name_;}
  void set_label_id(const unsigned int &id_) {label.id = id_;}

  const std::string& get_name() const {return name;}
  const graph::label& get_label() const {return label;}
  const unsigned int& get_id() const {return id;}
  
  unsigned int degree() const { return this->attached_links.size(); }

  bool operator==(const node&);
  bool operator==(const std::string&);

  node& operator=(const node&);

  //////////////// FOR SPIN MODELS ///////////////
  inline void set_spin_to_minus1(); 
  inline void set_spin_to_plus1();
  //////////////// FOR SPIN MODELS ///////////////
  
  friend std::ostream& operator<<(std::ostream&, const node&);
  friend std::ostream& operator<<(std::ostream&, const link&);
  
  ~node() {};
};

//////////////////////////////////////////////////////////////////////////
class link {
  friend class graph;
protected:
  node* node1;
  node* node2;
  std::string type;
  double weight;
public:
  link(const link&);
  link(node&, node&, const std::string& ="pp", double=0);
  link(node* = NULL, node* = NULL, const std::string& ="pp", double=0); //creates a link, increases the degrees of the nodes.

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
