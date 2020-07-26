#ifndef __GRAPH_HPP__
#define __GRAPH_HPP__

#include <list>
#include <iostream>
#include <fstream>
#include <functional>

class node;
class link;
class graph;
class node_pair;

#include "labeling.hpp"
#include "matrices/matrix.hpp"
#include "matrices/symm_matrix.hpp"
#include "matrices/A_matrix.hpp"
#include "matrices/W_matrix.hpp"
#include "matrices/row_vector.hpp"
#include "matrices/col_vector.hpp"
#include "aux_math.hpp"
#include "randomisation/prng.hpp"

class graph {
  friend class node_pair;
  
public:
  struct label; // Nested struct
  
protected:
  
  unsigned int NODES_ARRAY_SIZE;
  
  node* nodes;
  label* labels = NULL;
  std::list<link> links;
  
  unsigned int N_nodes; //Number of nodes (all elements of nodes array after N_nodes are treated as inactive)
  unsigned int N_labels; //Number of different labels
  unsigned int N_links; //Number of links
  
  bool remove_link(std::list<link>::iterator&); //Removes link with certain iterator
  // graph& remove_link(link*); //NOT IMPLEMENTED (see TODO 9).

  graph& remove_node_no_checks(node *p_node);

  label* find_label(const std::string&);
  
  // void move_node(const unsigned int&, const unsigned int&); //Moves node in the array redirecting links (needed for sorting)
  // void sort_nodes(); //Places all inactive nodes to the right and the rest sorts by name
public:
  
  graph();
  graph(const graph&);
  graph(const unsigned int& N); //Create graph reserving memory for N nodes
  graph(const std::string&, const unsigned int& N_additional_nodes=0); //Load labeled graph from efficient format to array, allocating more memory than needed to add nodes later on.
  // graph(const std::string& sif_file_name, const std::string& attrs_file_name, unsigned int, const unsigned int& N_additional_nodes=0);  //Load graph from Cytoscape .sif file and from labels assignment .attrs file
  graph(const matrix<bool>& adjacency_matrix); //Create simple non-directed graph from adjacency matrix

  void append_nodes_array(const unsigned int &N); //Reallocates memory for array of nodes adding extra space for N nodes
  
  unsigned int num_labels() const {return N_labels;}
  unsigned int num_nodes() const {return N_nodes;}
  unsigned int num_links() const {return N_links;}

  void add_node(const node&, const std::string&); //May be useful for creating intersecting graphs with shared nodes
  graph& add_node(const std::string &name, const std::string &label_name);
  graph& add_node(const std::string &name) { return add_node(name,"_"); }
  graph& add_node_no_checks(const std::string &name, const std::string &label_name="_");
  
  graph& remove_node_no_checks(const std::string&); //removes node by name, removes its links, but does not check if any labels should be removed
  graph& remove_node(const std::string&); //remove by name and remove the links
  // graph& remove_node_random(); //NOT IMPLEMENTED
  // graph& remove_random_node_with_label(const std::string&); //NOT IMPLEMENTED
  graph& remove_nodes_with_degree(const unsigned int &k);
  graph& remove_nodes_with_label(const std::string&);

  graph& add_link_no_checks(const unsigned int&, const unsigned int&, const std::string& ="pp", const double& weight=0);
  graph& add_link_no_checks(node* p_node1, node* p_node2, const std::string& ="pp", const double& weight=0);
  graph& add_link(const unsigned int&, const unsigned int&, const std::string& ="pp", const double& weight=0);
  void add_link(const std::string&, const std::string&, const std::string& ="pp", const double& weight=0);
  void add_link(const std::string&, const std::string&, const double&);
  graph& add_link(node*, node*, const std::string& ="pp", const double& weight=0);
    
  void remove_link(node*, node*, const std::string& type="pp", const double& weight=0);
  void remove_link(const std::string&, const std::string&, const std::string& ="pp", double=0);

  node* find_node(const std::string &node_name) const;
  node* find_node(const unsigned int &node_id) const; //Returns a node with certain position in the array (id)

  const node& get_node(const std::string &node_name) const;
  const node& get_node(const unsigned int &node_id) const; //Returns a node with certain position in the array (id)

  node_pair get_node_pair(unsigned int &id1, unsigned int &id2);
  node_pair random_node_pair(prng& rnd);
  col_vector<node_pair> random_node_pairs_col_vector(const unsigned int &N_pairs, prng& rnd); //Returns N random pairs of nodes

  link* get_link(const unsigned int &node1_id, const unsigned int &node2_id) const;
  link* get_link(const node *p_node1, const node *p_node2) const;

  bool is_simple() const; // Redundant! We already check to not add self-links!
  
  std::string* get_labels() const;
  matrix<std::string> get_labels_vec() const;

  void replicate_nodes(const unsigned int &replication_factor); //Proportionally creates many nodes with same labels, but different names

  A_matrix adjacency_matrix() const;

  void save_labels(const std::string&) const;
  
  void save(const std::string&) const;
  void save_to_sif_and_attrs(const std::string &name) const; //Saves the graph in Cytoscape recognised format

  graph& load(const std::string&, unsigned int N_additional_nodes=0); //Loads the graph with free space for additional nodes
  graph& load_from_sif(const std::string&, unsigned int N_additional_nodes=0);
  graph& load_from_sif_and_attrs(const std::string& sif_file_name, const std::string& attrs_file_name, const unsigned int& N_additional_nodes=0); // NEEDS TESTING
  
  graph& build_from_adjacency_matrix(const matrix<bool>& A); //Clears the graph and rebuilds with adjacency matrix A

  graph& set_ER(const double& p, prng& rnd);
  graph& set_ER_with_random_p(prng& rnd);

  graph& clear_links();
  graph& clear();

  graph& operator=(const graph&);

  //////////////////// COMPUTATIONAL FUNCTIONS BEGIN ////////////////////////
  const unsigned int& number_of_nodes() const { return N_nodes;  } //REDUNDANT (already have num_nodes)
  const unsigned int& number_of_links() const { return N_links;  } //REDUNDANT (already have num_links)
  const unsigned int& number_of_labels() const { return N_labels; } //REDUNDANT (already have num_labels)
  
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

  col_vector<double> clustering_coefficient_sequence() const;

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

  double degree_assortativity() const;

  //// Spin model functions begin ////
  const graph& randomize_spins(prng& rnd);
  const graph& run_Glauber_dynamics_no_fields(prng& rnd, const double &Temperature, const unsigned int &N_steps);
  const graph& run_Glauber_dynamics_with_fields(prng& rnd, const double &Temperature, const unsigned int &N_steps);
  
  short* spin_sequence() const; //Implicitly allocates memory!!!
  col_vector<short> spin_sequence_col_vector() const;
  matrix<double> interaction_matrix() const;  
  double average_magnetization() const;

  //Associative memory functions
  graph& set_spin_sequence(const col_vector<short>&);
  graph& set_spin_sequence(const row_vector<short>&);
  graph& set_spin_sequence_2D(const matrix<short>&);
  matrix<short> spin_sequence_2D(const unsigned int &dim_y, const unsigned int &dim_x) const; //Usefull to represent 2-D pictures to play with associative memory (with each node being a pixel)

  /* Multi-component Curie-Weiss model functions begin */
  const graph& construct_MCW_model(const col_vector<unsigned int>& N_components, const matrix<double>& Interactions, const col_vector<double>& Fields, prng& rnd);
  const graph& update_MCW_model(const col_vector<unsigned int>& N_components, const matrix<double>& Interactions, const col_vector<double>& Fields);
  bool check_MCW_model(col_vector<unsigned int>& Numbers_of_nodes_in_components) const; //Checks if nodes IDs are assigned in a way that the first N_1 IDs correspond to first component, then N_2 ID's to the second and so on. Also checks if the graph is appropriate MCW model with field node.
  
  col_vector<double> CW_components_magnetizations(col_vector<unsigned int>& Numbers_of_nodes_in_components) const;
  col_vector<double> CW_components_magnetizations(unsigned int& N_components) const;
  col_vector<double> CW_components_magnetizations() const;

  col_vector<double> CW_components_average_magnetizations(col_vector<unsigned int>& Numbers_of_nodes_in_components) const;
  col_vector<double> CW_components_average_magnetizations(unsigned int& N_components) const;
  col_vector<double> CW_components_average_magnetizations() const;
  /*  Multi-component Curie-Weiss model functions end  */

  //// Spin model functions end ////

  //// Random Graph generators begin ////
  graph& Metropolis_generator(double P(const matrix<bool>&), const unsigned int &N_nodes, const unsigned int &N_iters, prng& rnd, const bool &initialize_randomly=true); //Generates random graph from the ensemble defined by distribution P, which is passed in a form of a function of the adjacency matrix.
  graph& MF_Metropolis_generator(double P(const unsigned int&), const unsigned int &N_nodes, const unsigned int &N_iters, prng& rnd, const bool &initialize_randomly=true); //Generates random graph from the ensemble defined by distribution P, which is passed in a form of a function of the number of links.
  graph& GB_Metropolis_generator(double H(const matrix<bool>&), const unsigned int &N_nodes, const unsigned int &N_iters, prng& rnd, const bool &initialize_randomly=true, const double &temp=1); //Generates random graph from the ensemble defined by Gibbs-Boltzmann distribution, which Hamiltonian is passed as a function of the adjacency matrix.
  graph& MF_GB_Metropolis_generator(double H(const unsigned int&), const unsigned int &N_nodes, const unsigned int &N_iters, prng& rnd, const bool &initialize_randomly=true, const double &temp=1); //Generates random graph from the ensemble defined by Gibbs-Boltzmann distribution, which Hamiltonian is passed as a function of the number of links.
  
  graph& sample_p_star_model(const unsigned int &N_iters, prng& rnd, const col_vector<double> &T, const unsigned int &N_pairs_max=1, const bool &initialize_randomly=true, const double &temp=1); //Uses Metropolis Dynamics to sample from p-star model with parameters T. Doesn't construct adjacency matrix, PASS PARAMETERS: T[0] corresponds to t_1 and T[p-1] to t_p ! N_pairs is the number of distinclt pairs of nodes picket at random and being linked (unlinked) simultaneously at a single MD step.
  graph& sample_p_star_model_with_single_spin_Metropolis(const unsigned int &N_iters, prng& rnd, const col_vector<double> &T, const bool &initialize_randomly=true, const double &temp=1); //Uses Metropolis Dynamics to sample from p-star model with parameters T. Doesn't construct adjacency matrix, PASS PARAMETERS: T[0] corresponds to t_1 and T[p-1] to t_p ! N_pairs is the number of distinclt pairs of nodes picket at random and being linked (unlinked) simultaneously at a single MD step.
  
  unsigned int count_p_star_model_iterations_until(bool stopping_condition(const graph&), prng& rnd, const col_vector<double>& T, const unsigned int &N_pairs_max=1, const bool &initialize_randomly=true, const double &temp=1);
  //// Random Graph generators end ////

  //// Classification functions begin ////
  // Function below creates the correct labeling in the graph to start iterations over labelings and returns this labeling
  labeling initialize_labeling_for_classification(const unsigned short &N_classes, const std::string* class_names=NULL);
  graph& load_labeling(const labeling&, const std::string* class_names=NULL);

  matrix<double> link_probabilities() const;
  matrix<double> regularized_link_probabilities() const;
  double loglikelihood() const; //Computes likelihood of the label assignment in the max entropy ensemble with soft constraints on W_labels
  double regularized_loglikelihood() const; //Computes likelihood with small "perturbations" making sure there is no division by zero 
  std::list<labeling> modularity_classifier_precise(const unsigned short &num_classes, const std::string* class_names=NULL);
  std::list<labeling> ML_classifier_precise(const unsigned short &num_classes, const std::string* class_names=NULL);
  
  // graph& classify_nodes_precisely(const unsigned short &num_classes);
  // graph& classify_nodes_greedy_heuristic(const unsigned int &num_classes);
    
  //// Classification functions end ////  
  
  ///////////////////// COMPUTATIONAL FUNCTIONS END //////////////////////////
  
  ~graph();
  
  friend std::ostream& operator<<(std::ostream& ,const graph&);
};

#include "label.hpp"
#include "node.hpp"
#include "node_pair.hpp"
#include "link.hpp"


#include "matrices/matrix.tpp"
#include "matrices/symm_matrix.tpp"
#include "matrices/col_vector.tpp"
#include "matrices/row_vector.tpp"

#endif
