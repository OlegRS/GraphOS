#include "graphlib.h"

//   !!! THESE MACROS ARE NOT YET FULLY IMPLEMENTED !!!
// #define INTERPRETER_MODE //Don't exit() on errors, give ERROR message and continue (use with C++ interpreters like "cling" or "cint")

#define QUIET_MODE  // disables printing of some common warnings and info
// #define SILENT_MODE // disables printing of most warnings and info
#ifdef SILENT_MODE
#define QUIET_MODE
#endif

//////////////////////////////////////////////////////////////////
////////////////////////// graph /////////////////////////////////
//////////////////////////////////////////////////////////////////

graph::graph() : /*NODES_SORTED(false),*/ NODES_ARRAY_SIZE(0), nodes(NULL), labels(NULL), N_nodes(0), N_labels(0), N_links(0) {}

graph::graph(const unsigned int& N) : /*NODES_SORTED(false),*/ NODES_ARRAY_SIZE(N), labels(NULL), N_nodes(0), N_labels(0), N_links(0) {
  nodes = new node[N];
}

graph::graph(const graph& gr) : /*NODES_SORTED(gr.NODES_SORTED),*/ NODES_ARRAY_SIZE(gr.NODES_ARRAY_SIZE), N_nodes(gr.N_nodes), N_labels(gr.N_labels), N_links(gr.N_links) {
  nodes = new node[NODES_ARRAY_SIZE];
  for(unsigned int i=0; i<N_nodes; ++i)
    nodes[i] = gr.nodes[i];

  labels = new label[N_labels];
  for(unsigned int i=0; i<N_labels; ++i)
    labels[i] = gr.labels[i];

  // Copying links, redirecting node pointers
  for(std::list<link>::const_iterator it_links = gr.links.begin(); it_links!=gr.links.end(); ++it_links) {
    links.push_back(*it_links);
    links.back().node1 = nodes + (it_links->node1 - gr.nodes);
    links.back().node2 = nodes + (it_links->node2 - gr.nodes);
  }
}

//////////////////////////////////////////////////////////////////

graph::graph(const std::string &graph_file, unsigned int N_additional_nodes) {
#ifndef SILENT_MODE
  std::cerr << "***Loading the graph from \".graph\" file...\n";
#endif
  std::ifstream ifs(graph_file.c_str());
  
  if( !ifs.is_open() ) {
    std::cerr << "------------------------------------------------------------------------\n"
              << "ERROR: Cannot open the file \"" << graph_file << "\" to load the graph!\n"
              << "                      !!! GRAPH IS NOT LOADED !!!\n"
              << "------------------------------------------------------------------------\n";
#ifdef INTERPRETER_MODE
    return;
#else
    exit(1);
#endif
  }

  ifs >> N_nodes >> N_links >> N_labels;
  
  NODES_ARRAY_SIZE = N_nodes + N_additional_nodes;
  
  nodes = new node[NODES_ARRAY_SIZE];
  labels = new label[N_labels];

  std::string  buf1_str, buf2_str, buf3_str;
  unsigned int buf1_int, buf2_int;
  
  for(ifs >> buf1_int >> buf1_str >> buf2_int >> buf2_str; !ifs.eof(); ifs >> buf1_int >> buf1_str >> buf2_int >> buf2_str)
    if(ifs.get() != '\n') {
      if(nodes[buf1_int].name != buf1_str) {
        nodes[buf1_int] = node(buf1_int, buf1_str, buf2_int, buf2_str);
        labels[buf2_int] = buf2_str;
      }

      links.push_back(link());

      links.back().node1 = nodes + buf1_int;
      nodes[buf1_int].attached_links.push_back(--links.end());
        
      ifs >> links.back().type; //Link type

      ifs >> buf1_int >> buf1_str >> buf2_int >> buf2_str;
      
      if(nodes[buf1_int].name != buf1_str) {
        nodes[buf1_int] = node(buf1_int, buf1_str, buf2_int, buf2_str);
        labels[buf2_int] = buf2_str;
      }
          
      links.back().node2 = nodes + buf1_int;
      nodes[buf1_int].attached_links.push_back(--links.end());
      
      ifs >> links.back().weight; //Weight of the link
    }
  
    else {
      nodes[buf1_int] = node(buf1_int, buf1_str, buf2_int, buf2_str);
      labels[buf2_int] = buf2_str;
    }
  
  ifs.close();
  
#ifndef SILENT_MODE
  std::cerr << "*Graph of "<< N_nodes <<" nodes and "<< N_links <<" links is loaded from \"" << graph_file << "\" file.\n";
#endif
}

//////////////////////// !this function can be reimplemented for optimization of search! //////////////////
node* graph::find_node(const std::string &node_name) const {
  for(unsigned int i = 0; i<N_nodes; i++)
    if(nodes[i].name == node_name)
      return nodes + i;
  return NULL;
}
///////////////////////////////////////////////////////////////////////////////////////////////////////////
node* graph::find_node(const unsigned int &node_id) const {
  if(N_nodes > node_id) return nodes + node_id;
  
  std::cerr << "-----------------------------------------------\n"
            << "WARNING: Attempt to access nonexisting node!\n"
            << "-----------------------------------------------\n";
  return NULL;
}
///////////////////////////////////////////////////////////////////////////////////////////////////////////

const node& graph::get_node(const std::string &node_name) const {
  node* p_node = find_node(node_name);

  if(p_node) return *p_node;
  
  std::cerr << "----------------------------------------------------\n"
            << "FATAL ERROR: Attempt to access nonexisting node!\n"
            << "----------------------------------------------------\n";
  exit(1);
}

const node& graph::get_node(const unsigned int &node_id) const {
  if(N_nodes > node_id)
    return nodes[node_id];
  else {
    std::cerr << "----------------------------------------------------\n"
              << "FATAL ERROR: Attempt to access nonexisting node!\n"
              << "----------------------------------------------------\n";
    exit(1);
  } 
}
///////////////////////////////////////////////////////////////////////////////////////////////////////////

link* graph::get_link(const unsigned int &node1_id, const unsigned int &node2_id) const {
  for(std::list<std::list<link>::iterator>::iterator it_it_links = nodes[node1_id].attached_links.begin(); it_it_links!=nodes[node1_id].attached_links.end(); ++it_it_links)
    if( (*it_it_links)->node2 == nodes + node2_id || (*it_it_links)->node1 == nodes + node2_id )
      return &(*(*it_it_links));
  return NULL;
}
///////////////////////////////////////////////////////////////////////////////////////////////////////////

bool graph::is_simple() const {
  for(std::list<link>::const_iterator it_links = links.begin(); it_links!=links.end(); ++it_links)
    if(it_links->node1 == it_links->node2)
      return false;
  return true;
}
///////////////////////////////////////////////////////////////////////////////////////////////////////////

graph::label* graph::find_label(const std::string &label_name) {
  for(unsigned int i=0; i<N_labels; i++)
    if(labels[i] == label_name)
      return labels + i;
  return NULL;    
}
///////////////////////////////////////////////////////////////////////////////////////////////////////////

void graph::save_labels(const std::string &file_name) const {
  std::ofstream ofs( (file_name).c_str() );
  if(ofs.is_open())
    for(unsigned int i=0; i<N_labels; ++i)
      ofs << labels[i].name << '\n';
  else {
    std::cerr << "-----------------------------------------------\n"
              << "ERROR: Cannot open the file to save the labels!\n"
              << "-----------------------------------------------\n";
  }
}

void graph::save(const std::string &file_name) const {
  if(!N_nodes) {
    std::cerr << "-----------------------------------------------------\n"
              << "WARNING: Attempt to save empty graph is ignored!\n"
              << "-----------------------------------------------------\n";
    return;
  }
  
  std::ofstream ofs( (file_name + ".graph").c_str() );
  if(ofs.is_open()) {
    std::cerr << "***Saving the graph in \".graph\" format...\n";
    ofs << N_nodes << ' ' << N_links << ' ' << N_labels << '\n';

  
    for(std::list<link>::const_iterator it_links = links.begin(); it_links != links.end(); ++it_links)
      ofs << it_links->node1->id << ' ' << it_links->node1->name << ' ' << it_links->node1->label.id << ' ' << it_links->node1->label.name << ' ' << it_links->type << ' ' << it_links->node2->id << ' ' << it_links->node2->name << ' ' << it_links->node2->label.id << ' ' << it_links->node2->label.name << ' ' << it_links->weight << '\n';

    for(unsigned int i=0; i<N_nodes; i++)
      if(nodes[i].attached_links.empty())
        ofs << nodes[i].id << ' ' << nodes[i].name << ' ' << nodes[i].label.id << ' ' << nodes[i].label.name << '\n';

    ofs.close();
    std::cerr << "*Graph is saved to \"" << file_name + ".graph" << "\".\n";
  }
  else {
    std::cerr << "-------------------------------------------------------------------\n"
              << "ERROR: Cannot open the file to save the graph in \".graph\" format!\n"
              << "-------------------------------------------------------------------\n";
  }
}

void graph::save_to_sif_and_attrs(const std::string &file_name) const {
    if(!N_nodes) {
    std::cerr << "-----------------------------------------------------\n"
              << "WARNING: Attempt to save empty graph is ignored!\n"
              << "-----------------------------------------------------\n";
    return;
  }
    
  std::ofstream ofs_sif((file_name + ".sif").c_str());
  std::ofstream ofs_attrs((file_name + ".attrs").c_str());
  if(ofs_sif.is_open() && ofs_attrs.is_open()) {
    std::cerr << "***Saving the graph in \".sif - .attrs\" format...\n";
    
    for(std::list<link>::const_iterator it_links = links.begin(); it_links!=links.end(); ++it_links)
      ofs_sif << it_links->node1->name << ' ' << it_links ->type << ' ' << it_links->node2->name << '\n';

    ofs_attrs << "Labels\n";
    for(unsigned int i=0; i<N_nodes; i++) {
      ofs_attrs << nodes[i].name << " = " << nodes[i].label.name << '\n';
      if(!nodes[i].degree())
        ofs_sif << nodes[i].name << '\n';
    }

    ofs_sif.close();
    ofs_attrs.close();

    std::cerr << "*Graph is saved to \"" << file_name + ".sif" << "\" and \"" << file_name + ".attrs" <<"\" files.\n";
  }
  else {
    std::cerr << "---------------------------------------------------------------------------\n"
              << "ERROR: Cannot open the files to save the graph in \".sif - .attrs\" format!\n"
              << "---------------------------------------------------------------------------\n";
  }
}
/////////////////////////////////////////////////////////////////////////
graph& graph::load(const std::string &graph_file, unsigned int N_additional_nodes) {
  std::cerr << "***Loading the graph from \".graph\" file...\n";
  std::ifstream ifs(graph_file.c_str());
  
  if( !ifs.is_open() ) {
    std::cerr << "------------------------------------------------------------------------\n"
              << "ERROR: Cannot open the file \"" << graph_file << "\" to load the graph!\n"
              << "------------------------------------------------------------------------\n"
              << "*Graph is NOT loaded!!!\n";
    return *this;
  }
  
  if(N_nodes) {
    std::cerr << "--------------------------------------------------------------\n"
              << "WARNING: Graph is not empty, implicitly CLEARING the graph!\n"
              << "--------------------------------------------------------------\n";
    clear();
  }

  ifs >> N_nodes >> N_links >> N_labels;
  
  NODES_ARRAY_SIZE = N_nodes + N_additional_nodes;
  
  nodes = new node[NODES_ARRAY_SIZE];
  labels = new label[N_labels];

  std::string  buf1_str, buf2_str, buf3_str;
  unsigned int buf1_int, buf2_int;
  
  for(ifs >> buf1_int >> buf1_str >> buf2_int >> buf2_str; !ifs.eof(); ifs >> buf1_int >> buf1_str >> buf2_int >> buf2_str)
    if(ifs.get() != '\n') {
      if(nodes[buf1_int].name != buf1_str) {
        nodes[buf1_int] = node(buf1_int, buf1_str, buf2_int, buf2_str);
        labels[buf2_int] = buf2_str;
      }

      links.push_back(link());
          
      links.back().node1 = nodes + buf1_int;
      nodes[buf1_int].attached_links.push_back(--links.end());
        
      ifs >> links.back().type; //Link type

      ifs >> buf1_int >> buf1_str >> buf2_int >> buf2_str;
      
      if(nodes[buf1_int].name != buf1_str) {
        nodes[buf1_int] = node(buf1_int, buf1_str, buf2_int, buf2_str);
        labels[buf2_int] = buf2_str;
      }
          
      links.back().node2 = nodes + buf1_int;
      nodes[buf1_int].attached_links.push_back(--links.end());
      
      ifs >> links.back().weight; //Weight of the link
    }
  
    else {
      nodes[buf1_int] = node(buf1_int, buf1_str, buf2_int, buf2_str);
      labels[buf2_int] = buf2_str;
    }
  
  ifs.close();

  std::cerr << "*Graph of "<< N_nodes <<" nodes and "<< N_links <<" links is loaded from \"" << graph_file << "\" file.\n";
  
  return *this;
}
/////////////////////////////////////////////////////////////////////////
bool graph::remove_link(std::list<link>::iterator &it_link) {
  //Searching the iterator of the iterator of the link in first node
  bool iterator_found = false;
  std::list<std::list<link>::iterator>::iterator it_it_links;
  for(it_it_links = it_link->node1->attached_links.begin(); it_it_links != it_link->node1->attached_links.end(); ++it_it_links)
    if(*it_it_links == it_link) {
      iterator_found = true;
      break;
    }
  if(iterator_found) 
    it_link->node1->attached_links.erase(it_it_links); //Deleting iterator to the link from first node
  else
    return false;

  //Seeking for the poiner to the link in second node
  iterator_found = false;
  for(it_it_links = it_link->node2->attached_links.begin(); it_it_links != it_link->node2->attached_links.end(); ++it_it_links)
    if(*it_it_links == it_link) {
      iterator_found = true;
      break;
    }
  if(iterator_found)
    it_link->node2->attached_links.erase(it_it_links); //Deleting iterator to this link from the second node
  else return false;
  
  links.erase(it_link); //Deleting link
  N_links--;
  return true;
}
/////////////////////////////////////////////////////////////////////////////////////////////////////////

void graph::remove_link(node* p_node1, node* p_node2, const std::string &type, double weight) {
  link l(p_node1, p_node2, type, weight);

  bool link_found = false;
  std::list<std::list<link>::iterator>::iterator it_it_links = p_node1->attached_links.begin();
  for( ; it_it_links != p_node1->attached_links.end(); ++it_it_links)
    if(*(*it_it_links) == l) {
      link_found = true;
      break;
    }
  
  if(link_found) {
    p_node1 -> attached_links.erase(it_it_links); //Deleting iterator to the link from first node

    //Seekink for the poiner to the link in second node
    for(it_it_links = p_node2->attached_links.begin(); it_it_links != p_node2->attached_links.end(); ++it_it_links)
      if(*(*it_it_links) == l)
        break;

    links.erase(*it_it_links); //Deleting link
    
    p_node2 -> attached_links.erase(it_it_links); //Deleting iterator to this link from the second node
    
    N_links--;
  }
  else {
    std::cerr << "--------------------------------------------------------\n"
              << "WARNING: Attempt to remove non-existing link is ignored!\n"
              << "--------------------------------------------------------\n";
    return;
  }
}
////////////////////////////////////////////////////////////////////////
void graph::remove_link(const std::string &node1_name, const std::string &node2_name, const std::string &type, double weight) {
  node *p_node1 = find_node(node1_name);
  node *p_node2 = find_node(node2_name);
  
  if(!p_node1) {
    std::cerr << "--------------------------------------------------------------------------------\n"
              << "ERROR: Attempt to remove a link of non-existing node \"" << node1_name << "\"!\n"
              << "--------------------------------------------------------------------------------\n";
    exit(1);
  }

  if(!p_node2) {
    std::cerr << "--------------------------------------------------------------------------------\n"
              << "ERROR: Attempt to remove a link of non-existing node \"" << node2_name << "\"!\n"
              << "--------------------------------------------------------------------------------\n";
    exit(1);
  }

  remove_link(p_node1, p_node2, type, weight);
}
/////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////

void graph::add_link_no_checks(const unsigned int &node1_id, const unsigned int &node2_id, std::string type, double weight) {
  if(!N_nodes) {
    std::cerr << "------------------------------------------------\n"
              << "ERROR: Attempt to add a link to empty graph!\n"
              << "------------------------------------------------\n";
    exit(1);
  }
  
  links.push_back(link(&nodes[node1_id], &nodes[node2_id], type, weight));
  nodes[node1_id].attached_links.push_back(--links.end());
  nodes[node2_id].attached_links.push_back(--links.end());

  N_links++;
}
/////////////////////////////////////////////////////////////////////////
graph& graph::add_link(const unsigned int &i, const unsigned int &j, std::string type, double weight) {
  if(i == j) {
    std::cerr << "--------------------------------------------------\n"
              << "WARNING: Attempt to create a self-link is ignored!\n"
              << "--------------------------------------------------\n";
    return *this;
  }

  link l(nodes+i, nodes+j, type, weight);
  bool link_found = false;
  for(std::list<std::list<link>::iterator>::iterator it_it_links = nodes[i].attached_links.begin(); it_it_links!=nodes[i].attached_links.end(); ++it_it_links)
    if(*(*it_it_links) == l) {
      link_found = true;
      break;
    }
  if(!link_found) {
    links.push_back(l);
    nodes[i].attached_links.push_back(--links.end());
    nodes[j].attached_links.push_back(--links.end());
      
    ++N_links;
  }
  else {
    std::cerr << "--------------------------------------------------------------------------\n"
              << "WARNING: Attempt to add existing link \""<< l <<"\" is ignored!\n"
              << "--------------------------------------------------------------------------\n";
    return *this;
  }    
  return *this;
}

void graph::add_link(node* node1, node* node2, const std::string &type, double weight) {
  
  if(node1 == node2) {
    std::cerr << "--------------------------------------------------\n"
              << "WARNING: Attempt to create a self-link is ignored!\n"
              << "--------------------------------------------------\n";
    return;
  }
  
  link l(node1, node2, type, weight);
  bool link_found = false;
  for(std::list<std::list<link>::iterator>::iterator it_it_links = node1->attached_links.begin(); it_it_links!=node1->attached_links.end(); ++it_it_links)
    if(*(*it_it_links) == l) {
      link_found = true;
      break;
    }
  if(!link_found) {
    links.push_back(l);
    node1->attached_links.push_back(--links.end());
    node2->attached_links.push_back(--links.end());
      
    ++N_links;
  }
  else {
    std::cerr << "--------------------------------------------------------------------------\n"
              << "WARNING: Attempt to add existing link \""<< l <<"\" is ignored!\n"
              << "--------------------------------------------------------------------------\n";
    return;
  }    
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void graph::add_link(const std::string &node1_name, const std::string &node2_name, const std::string &type, double weight) {
  
  node *p_node1 = find_node(node1_name);
  node *p_node2 = find_node(node2_name);
  
  if(!p_node1) {
    std::cerr << "--------------------------------------------------------------------------------\n"
              << "ERROR: Attempt to attach a link to non-existing node \"" << node1_name << "\"!\n"
              << "--------------------------------------------------------------------------------\n";
    exit(1);
  }
  
  if(!p_node2) {
    std::cerr << "--------------------------------------------------------------------------------\n"
              << "ERROR: Attempt to attach a link to non-existing node \"" << node2_name << "\"!\n"
              << "--------------------------------------------------------------------------------\n";
    exit(1);
  }
  
  add_link(p_node1, p_node2, type, weight);
}

void graph::add_link(const std::string &node1_name, const std::string &node2_name, double &weight) {
  add_link(node1_name, node2_name, "pp", weight);
}
////////////////////////////////////////////////////////////////////////////

graph& graph::add_node_no_checks(const std::string &name, const std::string &label_name) {
  if(N_nodes < NODES_ARRAY_SIZE) {
    nodes[N_nodes].name = name;
    nodes[N_nodes].label = label_name;
    nodes[N_nodes].id = N_nodes; //Placing new node to the end of active nodes

    //Setting the id of the label
    label *p_label = find_label(label_name);
    if(p_label)
      nodes[N_nodes].label.id = p_label->id;
    else { //Reallocate memory for labels
      label* labels_ = new label[N_labels+1];
      for(unsigned int i=0; i<N_labels; i++)
        labels_[i] = labels[i];
      delete[] labels;
      labels = labels_;
      labels[N_labels].name = label_name;
      labels[N_labels].id = nodes[N_nodes].label.id = N_labels;
      N_labels++;
    }
  }
  else {
    std::cerr << "-----------------------------------------------------------\n"
              << "WARNING: Implicitly reallocating memory for 100 more nodes!\n"
              << "-----------------------------------------------------------\n";

    append_nodes_array(100);
    
    nodes[N_nodes].name = name;
    nodes[N_nodes].label = label_name;
    nodes[N_nodes].id = N_nodes; //Placing new node to the end of active nodes

    //Setting the id of the label
    label *p_label = find_label(label_name);
    if(p_label)
      nodes[N_nodes].label.id = p_label->id;
    else { //Reallocate memory for labels
      label* labels_ = new label[N_labels+1];
      for(unsigned int i=0; i<N_labels; i++)
        labels_[i] = labels[i];
      delete[] labels;
      labels = labels_;
      labels[N_labels].name = label_name;
      labels[N_labels].id = nodes[N_nodes].label.id = N_labels;
      N_labels++;
    }
  }

  N_nodes++;
  return *this;
}
////////////////////////////////////////////////////////////////////////////
graph& graph::add_node(const std::string &name, const std::string &label_name) {
  for(unsigned int i = 0; i<N_nodes; i++)
    if(nodes[i].name == name) {
      std::cerr << "---------------------------------------------------------------------\n"
                << "WARNING: Attempt to add existing node: \"" << name << "\" is ignored!\n"
                << "---------------------------------------------------------------------\n";
      return *this;
    }
  add_node_no_checks(name, label_name);
  return *this;
}
////////////////////////////////////////////////////////////////////////////
graph& graph::remove_node_no_checks(node *p_node) {
  if(p_node) {
    while(!p_node->attached_links.empty())
      remove_link( *p_node->attached_links.begin() );

    //Copying the last node in the array in place of deleted node
    *p_node = nodes[--N_nodes];
    p_node->id = p_node - nodes;

    //Redirecting links
    for(std::list<std::list<link>::iterator>::iterator it_it_links = p_node->attached_links.begin(); it_it_links != p_node->attached_links.end(); ++it_it_links)
      if((*it_it_links)->node1 == nodes + N_nodes)
        (*it_it_links)->node1 = p_node;
      else
        (*it_it_links)->node2 = p_node;

    nodes[N_nodes].attached_links.clear();  
  }
  else {
    std::cerr << "---------------------------------------------\n"
              << "WARNING: Attempt to remove non-existing node!\n"
              << "---------------------------------------------\n";
  }
  
  return *this;
}

graph& graph::remove_node_no_checks(const std::string &name) {
  node *p_node = find_node(name);
  return remove_node_no_checks(p_node);
}

graph& graph::remove_node(const std::string &name) {
  node *p_node = find_node(name);
      
  if(p_node) {
    label lbl = p_node->label;
    remove_node_no_checks(p_node);
    
    //Checking if any labels now absent in the graph
    //!!! EXTREMELY INEFFICIENT !!!
    // This will have to be redone when nodes will store only pointers to labels!
    bool label_found = false;
    for(unsigned int i=0; i<N_nodes; ++i)
      if(nodes[i].label == lbl)
        label_found = true;

    if(!label_found) {
      label_found = false;
      label* labels_new = new label[N_labels-1];
      unsigned int j=0;
      for(unsigned int i=0; i<N_labels; i++)
        if( labels[i].name != lbl.name ) {
          labels_new[j].name = labels[i].name;
          labels_new[j].id = j;
          if(label_found) {
            for(unsigned int l=0; l<N_nodes; l++)
              if(nodes[l].label == labels_new[j])
                nodes[l].label.id = j;
          }
          ++j;
        }
        else label_found = true;
      delete[] labels;
      labels = labels_new;
      --N_labels;
    }
  }
  else {
    std::cerr << "---------------------------------------------\n"
              << "WARNING: Attempt to remove non-existing node!\n"
              << "---------------------------------------------\n";
  }

  return *this; 
}

graph& graph::remove_nodes_with_degree(const unsigned int &k) {//!!!NEEDS OPTIMIZATION!!!
  bool node_found = false;
  
  unsigned int i=0;
  for(; i<N_nodes; ++i)
    if(nodes[i].degree() == k) {
      node_found = true;
      break;
    }

  if(node_found) {
    remove_node_no_checks(nodes+i);
    for(; i<N_nodes; i++) 
      if(nodes[i].degree() == k)
        remove_node_no_checks(nodes + (i--));
    
    // Now erasing labels that were deleted (if they were)
    //!!! EXTREMELY INEFFICIENT !!!
    // This will have to be redone when nodes will store only pointers to labels!
    
    std::list<std::string> labels_list;
    
    unsigned int N_labels_new = 0;
    for(unsigned int i=0; i<N_labels; i++)
      for(unsigned int j=0; j<N_nodes; j++)
        if(labels[i].name == nodes[j].label.name) {
          labels_list.push_back(labels[i].name);
          N_labels_new++;
          break;
        }

    if(N_labels_new != N_labels) {
      delete[] labels;
      labels = new label[N_labels_new];
      N_labels = N_labels_new;

      unsigned int i=0;
      for(std::list<std::string>::iterator it_labels = labels_list.begin(); it_labels != labels_list.end(); ++it_labels) {
        labels[i].name = *it_labels;
        labels[i].id = i;
        for(unsigned int j=0; j<N_nodes; ++j)
          if(nodes[j].label.name == labels[i].name)
            nodes[j].label.id = i;
        ++i;
      }
    }
  }
  else {
    std::cerr << "---------------------------------------------------------------------------------------\n"
              << "WARNING: Attempt to remove all nodes with non-existing degree " << k <<" is ignored!\n"
              << "---------------------------------------------------------------------------------------\n";
  }
      
  return *this;
}

graph& graph::remove_nodes_with_label(const std::string &lbl) {//!!!NEEDS OPTIMIZATION!!!
  bool label_found = false;
  
  for(unsigned int i=0; i<N_nodes; ++i)
    if(nodes[i].label == lbl) {
      if(!label_found) label_found = true;
      remove_node_no_checks(nodes + i--);
    }

  if(label_found) {
    // Now removing the label
    //!!! INEFFICIENT !!!
    // This will have to be redone when nodes will store only pointers to labels!
    label_found = false;
    label* labels_new = new label[N_labels-1];
    unsigned int j=0;
    for(unsigned int i=0; i<N_labels; i++)
      if( labels[i].name != lbl ) {
        labels_new[j].name = labels[i].name;
        labels_new[j].id = j;
        if(label_found) {
          for(unsigned int l=0; l<N_nodes; l++)
            if(nodes[l].label == labels_new[j])
              nodes[l].label.id = j;
        }
        ++j;
      }
      else label_found = true;
    
    delete[] labels;
    labels = labels_new;
    --N_labels;
  }
  else {
    std::cerr << "----------------------------------------------------------------------------------------\n"
              << "WARNING: Attempt to remove all nodes with non-existing label " << lbl <<" is ignored!\n"
              << "----------------------------------------------------------------------------------------\n";
  }
  
  return *this;
}

//////////////////////////////////////////////////////////////////////////
std::string* graph::get_labels() const {
  std::string *labels_array= new std::string[N_labels];
  for(unsigned int i=0; i<N_labels; i++)
    labels_array[i] = labels[i].name;
  return labels_array;
}

matrix<std::string> graph::get_labels_vec() const {
  matrix<std::string> M(N_labels,1);
  for(unsigned int i=0; i<N_labels; i++)
    M[i][0] = labels[i].name;
  return M;
}
//////////////////////////////////////////////////////////////////////////
void graph::replicate_nodes(const unsigned int &replication_factor) {
  if(N_nodes * replication_factor > NODES_ARRAY_SIZE)
    append_nodes_array(N_nodes * replication_factor - NODES_ARRAY_SIZE);

  node *p_node;
  for(unsigned int i=1; i<replication_factor; ++i)
    for(unsigned int j=0; j<N_nodes; ++j) {
      p_node = nodes + N_nodes*i + j;
      p_node -> name = std::to_string(i) + "r" + nodes[j].name;
      p_node -> label = nodes[j].label;
      p_node -> id = N_nodes*i + j;      
    }

  N_nodes = N_nodes * replication_factor;
    
}

matrix<bool> graph::adjacency_matrix() const {
  matrix<bool> M(N_nodes, N_nodes);
  
  for(unsigned int i=0; i<N_nodes; ++i)
    for(unsigned int j=0; j<N_nodes; ++j)
      M[i][j] = false;
  
  for(std::list<link>::const_iterator it_links = links.begin(); it_links!=links.end(); ++it_links)
    M[it_links->node1-nodes][it_links->node2-nodes] = M[it_links->node2-nodes][it_links->node1-nodes] = 1;
  
  return M;
}

///////////////////////////////////////////////////////////////////////////
////////////////// GRAPH COMPUTATIONAL FUNCTIONS BEGIN ////////////////////
///////////////////////////////////////////////////////////////////////////
unsigned int graph::max_degree() const {
  if(N_nodes) {
    unsigned int deg_max = nodes[0].degree();
    for(unsigned int i=1; i<N_nodes; i++)
      if(nodes[i].degree() > deg_max)
        deg_max = nodes[i].degree();
    return deg_max;
  }
  else {
    std::cerr << "----------------------------------------------------------------------\n"
              << "WARNING: Call of \"max_degree()\" function for empty graph is ignored!\n"
              << "----------------------------------------------------------------------\n";
    return 0;
  }
}

unsigned int graph::min_degree() const {
  if(N_nodes) {
    unsigned int deg_min = nodes[0].degree();
    for(unsigned int i=1; i<N_nodes; i++)
      if(nodes[i].degree() < deg_min)
        deg_min = nodes[i].degree();
    return deg_min;
  }
  else {
    std::cerr << "----------------------------------------------------------------------\n"
              << "WARNING: Call of \"max_degree()\" function for empty graph is ignored!\n"
              << "----------------------------------------------------------------------\n";
    return 0;
  }
}

unsigned int* graph::degree_sequence() const {
#ifndef QUIET_MODE
  std::cerr << "--------------------------------------------------------------------------------------\n"
            << "WARNING: Implicit memory allocation for degree_sequence! Don't forget to delete[] it!\n"
            << "--------------------------------------------------------------------------------------\n";
#endif
  
  unsigned int* degrees = new unsigned int[N_nodes];
  for(unsigned int i=0; i<N_nodes; ++i)
    degrees[i] = nodes[i].degree();
  return degrees;
}

col_vector<unsigned int> graph::degree_sequence_col_vec() const {
  col_vector<unsigned int> degrees(N_nodes);
  for(unsigned int i=0; i<N_nodes; ++i)
    degrees[i] = nodes[i].degree();
  return degrees;

}

col_vector<node*> graph::nodes_with_max_degree_col_vec() const {
  if(N_nodes) {
    unsigned int* degrees = degree_sequence();
    
    unsigned int k_max = max_degree(), count = 0;   
    for(unsigned int i=0; i<N_nodes; i++)
      if(degrees[i] == k_max)
        ++count;

    col_vector<node*> V_nodes(count);
    for(unsigned int i=0; count != 0; ++i)
      if(degrees[i] == k_max)
        V_nodes[--count] = nodes + i;
    
    delete[] degrees;
    return V_nodes;
  }
  else {
    std::cerr << "-----------------------------------------------------------------------------------------\n"
              << "WARNING: Call of \"node_with_max_degree_col_vec()\" function from empty graph is ignored!\n"
              << "-----------------------------------------------------------------------------------------\n";
    return col_vector<node*>();
  }
}

col_vector<node*> graph::nodes_with_min_degree_col_vec() const {
  if(N_nodes) {
    unsigned int* degrees = degree_sequence();
    
    unsigned int k_min = min_degree(), count = 0;   
    for(unsigned int i=0; i<N_nodes; i++)
      if(degrees[i] == k_min)
        ++count;
    
    col_vector<node*> V_nodes(count);
    for(unsigned int i=0; count != 0; ++i)
      if(degrees[i] == k_min)
        V_nodes[--count] = nodes + i;
    
    delete[] degrees;
    return V_nodes;
  }
  else {
    std::cerr << "-----------------------------------------------------------------------------------------\n"
              << "WARNING: Call of \"node_with_min_degree_col_vec()\" function from empty graph is ignored!\n"
              << "-----------------------------------------------------------------------------------------\n";
    return col_vector<node*>();
  }
}

col_vector<node*> graph::nodes_with_degree_col_vec(unsigned int &k) const {
  if(N_nodes) {
    unsigned int* degrees = degree_sequence();
    
    unsigned int count = 0;   
    for(unsigned int i=0; i<N_nodes; i++)
      if(degrees[i] == k)
        ++count;
    
    col_vector<node*> V_nodes(count);
    for(unsigned int i=0; count != 0; ++i)
      if(degrees[i] == k)
        V_nodes[--count] = nodes + i;
    
    delete[] degrees;
    return V_nodes;
  }
  else {
    std::cerr << "-------------------------------------------------------------------------------------\n"
              << "WARNING: Call of \"node_with_degree_col_vec()\" function from empty graph is ignored!\n"
              << "-------------------------------------------------------------------------------------\n";
    return col_vector<node*>();
  }
}

std::string* graph::label_sequence() const {
#ifndef QUIET_MODE
  std::cerr << "-------------------------------------------------------------------------------------\n"
            << "WARNING: Implicit memory allocation for label_sequence! Don't forget to delete[] it!\n"
            << "-------------------------------------------------------------------------------------\n";
#endif
  
  std::string* label_sequence = new std::string[N_nodes];
  for(unsigned int i=0; i<N_nodes; ++i)
    label_sequence[i] = labels[i].name;
  return label_sequence;
}

col_vector<std::string> graph::label_sequence_col_vector() const {
  col_vector<std::string> V(N_nodes);
  for(unsigned int i=0; i<N_nodes; ++i)
    V[i] = nodes[i].label.name;
  return V;
}

double graph::average_degree() const {
  if(N_nodes)
    return 2*N_links/(double)N_nodes;
  else {
    std::cerr << "--------------------------------------------------------------------------\n"
              << "WARNING: Call of \"average_degree()\" function for empty graph is ignored!\n"
              << "--------------------------------------------------------------------------\n";
    return 0;
  }
}

double graph::degree_distribution(const unsigned int &k) const {
  if(N_nodes) {
    unsigned int count = 0;
    for(unsigned int i=0; i<N_nodes; ++i)
      if(nodes[i].degree() == k)
        count++;
    return count/(double)N_nodes;
  }
  else {
    std::cerr << "----------------------------------------------------------------------------------------------------\n"
              << "WARNING: Call of \"degree_distribution(const unsigned int &k)\" function for empty graph is ignored!\n"
              << "----------------------------------------------------------------------------------------------------\n";
    return 0;
  }
}

double* graph::degree_distribution() const {
  if(N_nodes) {
    unsigned int deg_min = min_degree();
    unsigned int deg_max = max_degree();

    double *deg_distr = new double[deg_max+1];
    for(unsigned int i=0; i<=deg_min; i++)
      deg_distr[i] = 0;
    for(unsigned int i=deg_min; i<=deg_max; i++)
      deg_distr[i] = degree_distribution(i);
    return deg_distr;
  }
  else {
    std::cerr << "----------------------------------------------------------------------------------\n"
              << "WARNING: Call of \"degree_distribution()\" function for empty graph is ignored!\n"
              << "----------------------------------------------------------------------------------\n";
    return NULL;
  }
}

col_vector<double> graph::degree_distribution_col_vector() const {
    if(N_nodes) {
    unsigned int deg_min = min_degree();
    unsigned int deg_max = max_degree();

    col_vector<double> V(deg_max+1);
    for(unsigned int i=0; i<deg_min; i++)
      V[i] = 0;
    for(unsigned int i=deg_min; i<=deg_max; i++)
      V[i] = degree_distribution(i);
    return V;
  }
  else {
    std::cerr << "----------------------------------------------------------------------------------\n"
              << "WARNING: Call of \"degree_distribution()\" function for empty graph is ignored!\n"
              << "----------------------------------------------------------------------------------\n";
    return col_vector<double>();
  }
}


double graph::label_distribution(std::string &label_name) const {
  unsigned int count = 0;
  for(unsigned int i=0; i<N_nodes; i++)
    if(nodes[i].label == label_name)
      count++;
  return count/(double)N_nodes;
}

double graph::label_distribution(unsigned int &label_id) const {
  unsigned int count = 0;
  for(unsigned int i=0; i<N_nodes; i++)
    if(nodes[i].label.id == label_id)
      count++;
  return count/(double)N_nodes;
}

double* graph::label_distribution() const {
  double *lbl_distr = new double[N_labels];
  for(unsigned int i=0; i<N_labels; i++)
    lbl_distr[i] = label_distribution(i);
  return lbl_distr;
}

matrix<double> graph::label_distribution_vec() const {
  matrix<double> M(N_labels,1);
  for(unsigned int i=0; i<N_labels; i++)
    M[i][0] = label_distribution(i);
  return M;
}

double graph::joint_distribution_of_connected_node_labels(const std::string &label_1_name, const std::string &label_2_name) const {
  unsigned int count = 0;
  if(label_1_name != label_2_name) {
    for(std::list<link>::const_iterator it_links = links.begin(); it_links!=links.end(); ++it_links)
      if( (it_links->node1->label == label_1_name && it_links->node2->label == label_2_name) || (it_links->node1->label == label_2_name && it_links->node2->label == label_1_name) )
        ++count;
    return count/((double)N_links*2);
  }
  else {
   for(std::list<link>::const_iterator it_links = links.begin(); it_links!=links.end(); ++it_links)
      if( it_links->node1->label == label_1_name && it_links->node2->label == label_2_name )
        ++count;
    return count/(double)N_links;
  }
}

double graph::joint_distribution_of_connected_node_labels(const unsigned int &label_1_id, const unsigned int &label_2_id) const {
  unsigned int count = 0;
  if(label_1_id != label_2_id) {
    for(std::list<link>::const_iterator it_links = links.begin(); it_links!=links.end(); ++it_links)
      if( (it_links->node1->label.id == label_1_id && it_links->node2->label.id == label_2_id) || (it_links->node1->label.id == label_2_id && it_links->node2->label.id == label_1_id) )
        ++count;
    return count/((double)N_links*2);
  }
  else {
    for(std::list<link>::const_iterator it_links = links.begin(); it_links!=links.end(); ++it_links)
      if( it_links->node1->label.id == label_1_id && it_links->node2->label.id == label_2_id )
        ++count;
    return count/((double)N_links);
  }
}

W_matrix graph::joint_distribution_of_connected_node_labels() const {
  if(!N_nodes) {
#ifndef SILENT_MODE
    std::cerr << "------------------------------------------------------------------------------------------------------\n"
              << "WARNING: Call of \"joint_distribution_of_connected_node_labels\" function for empty graph is ignored!\n"
              << "------------------------------------------------------------------------------------------------------\n";
#endif
    return W_matrix();
  }
#ifndef QUIET_MODE
  std::cerr << "***Computing full matrix of joint distribution of connected nodes labels...\n";
#endif
  W_matrix W(N_labels);
  for(unsigned int i=0; i<N_labels; ++i)
    for(unsigned int j=0; j<i; ++j)
      W[i][j] = W[j][i] = joint_distribution_of_connected_node_labels(i, j);

  for(unsigned int i=0; i<N_labels; ++i)
    W[i][i] = joint_distribution_of_connected_node_labels(i, i);
#ifndef QUIET_MODE
  std::cerr << "*Joint distribution of connecred node labels is computed.\n";
#endif
  return W;
}

double graph::joint_distribution_of_connected_node_degrees(const unsigned int &k1, const unsigned int &k2) const {
  //Storing degree sequence for acceleration
  unsigned int *degree_sequence = new unsigned int[N_nodes];
  for(unsigned int i=0; i<N_nodes; i++)
    degree_sequence[i] = nodes[i].degree();

  unsigned int count = 0;
  
  if(k1!=k2) {
    for(std::list<link>::const_iterator it_links = links.begin(); it_links!=links.end(); ++it_links)
      if( (degree_sequence[it_links->node1-nodes] == k1 && degree_sequence[it_links->node2-nodes] == k2) || (degree_sequence[it_links->node1-nodes] == k2 && degree_sequence[it_links->node2-nodes] == k1) )
        ++count;
    
    delete[] degree_sequence;
    return count/((double)N_links*2);
  }
  else{
    for(std::list<link>::const_iterator it_links = links.begin(); it_links!=links.end(); ++it_links)
      if(degree_sequence[it_links->node1-nodes] == k1 && degree_sequence[it_links->node2-nodes] == k2)
        ++count;

    delete[] degree_sequence;
    return count/(double)N_links;
  }
}

W_matrix graph::joint_distribution_of_connected_node_degrees() const {
  if(!N_nodes) {
    std::cerr << "------------------------------------------------------------------------------------------------------\n"
              << "WARNING: Call of \"joint_distribution_of_connected_node_degrees\" function for empty graph is ignored!\n"
              << "------------------------------------------------------------------------------------------------------\n";
    return W_matrix();
  }

  std::cerr << "***Computing full matrix of joint distribution of connected nodes degrees...\n";
  
  unsigned int k_max = max_degree();
  unsigned int k_min = min_degree();
  W_matrix W(k_max+1);

  //Setting all elements below k_min to zero
  for(unsigned int i=0; i<k_min; ++i)
    for(unsigned int j=0; j<=k_max; ++j)
      W[i][j] = W[j][i] = 0;
  
  //Filling non-zero elements
  for(unsigned int i=k_min; i<=k_max; ++i)
    for(unsigned int j=k_min; j<=i; ++j)
      W[i][j] = W[j][i] = joint_distribution_of_connected_node_degrees(i, j);

  std::cerr << "\n*Joint distribution of connecred node degrees is computed.\n";
  return W;
}

double graph::modularity() const {
#ifndef QUIET_MODE
  std::cerr << "***Computing modularity of the graph...\n";
#endif
  unsigned int *K = degree_sequence();
  
  double count1 = 0;
  for(std::list<link>::const_iterator it_links = links.begin(); it_links != links.end(); ++it_links)
    if(it_links->node1->label.id == it_links->node2->label.id)
      count1 += 2;

  double count2 = 0;
  for(unsigned int i=0; i<N_nodes; ++i)
    for(unsigned int j=0; j<i; ++j)
      if(nodes[i].label.id == nodes[j].label.id)
        count2 += K[i]*K[j];

  delete[] K;
#ifndef QUIET_MODE
  std::cerr << "*Modularity is calculated.\n";
#endif
  return (count1 - count2/N_links) / (N_links*4);
}

double graph::relative_modularity() const {
  return modularity() / (1 - 1 /double(N_labels) + 1/double(N_nodes)) * 2;
}

double graph::assortativity() const {
  return joint_distribution_of_connected_node_degrees().assortativity();
}

///////////////////// SPIN DYNAMICS BEGIN ///////////////////////////////////

short* graph::spin_sequence() const{
#ifndef QUIET_MODE
  std::cerr << "-------------------------------------------------------------------------------------\n"
            << "WARNING: Implicit memory allocation for label_sequence! Don't forget to delete[] it!\n"
            << "-------------------------------------------------------------------------------------\n";
#endif

  short* spin_sequence = new short[N_nodes];
  for(unsigned int i=0; i<N_nodes; ++i)
    spin_sequence[i] = std::stoi(labels[i].name);
  return spin_sequence;
}

col_vector<short> graph::spin_sequence_col_vector() const {
  col_vector<short> V(N_nodes);
  for(unsigned int i=0; i<N_nodes; ++i)
    V[i] = std::stoi(nodes[i].label.name);
  return V;
}

matrix<double> graph::interaction_matrix() const {
  matrix<double> M(N_nodes, N_nodes);
  
  for(unsigned int i=0; i<N_nodes; ++i)
    for(unsigned int j=0; j<N_nodes; ++j)
      M[i][j] = 0;
  
  for(std::list<link>::const_iterator it_links = links.begin(); it_links!=links.end(); ++it_links)
    M[it_links->node1-nodes][it_links->node2-nodes] = M[it_links->node2-nodes][it_links->node1-nodes] = it_links->weight;
  
  return M;
}

const graph& graph::randomize_spins(const int &seed) {
  srand(seed);
  if(N_labels!=2 || labels[0].name!="-1" || labels[0].id!=0 || labels[1].name!="1" || labels[1].id!=1) {
    if(labels != NULL)
      delete[] labels;
    N_labels=2;
    labels = new label[N_labels];
    labels[0].id = 0;
    labels[0].name = "-1";
    labels[1].id = 1;
    labels[1].name = "1";
  }
  for(unsigned int i =0; i<N_nodes; i++) {
    if(nodes[i].name != "FIELD_NODE")
      nodes[i].label.name = std::to_string(2*(short)(nodes[i].label.id = rand()%2)-1);
  }
  return *this;
}

const graph& graph::simulate_Glauber_Dynamics_no_fields(const int &seed, const double &T, const unsigned int &N_steps) {
#ifndef QUIET_MODE
  std::cerr << "*** Simulating Glauber Dynamics...\n";
#endif
  //////////////// CONSISTENCY CHECK //////////////
  if(N_labels > 2) {
    std::cerr << "---------------------------------------------------------------------------------------\n"
              << "ERROR: GRAPH IS INSUFFICIENT FOR GLAUBER DYNAMICS (has more than 2 different labels)!\n"
              << "---------------------------------------------------------------------------------------\n";
    return *this;
  };
  
  for(unsigned int i=0; i<N_nodes; ++i)
    if(nodes[i].label.name != "1" && nodes[i].label.name != "-1") {
      std::cerr << "-------------------------------------------------------------------------------------------------------\n"
                << "ERROR: GRAPH IS INSUFFICIENT FOR GLAUBER DYNAMICS (has labels that are different from \"1\" and \"-1\")!\n"
                << "DEFECTED LABEL OF NODE " << nodes[i].name << " (id="<< i <<") is " << nodes[i].label.name << '\n'
                << "-------------------------------------------------------------------------------------------------------\n";
      return *this;
  };

  ////////////// SIMULATING DYNAMICS ////////////////
  srand(seed);
  double RAND_MAX_DOUBLE = double(RAND_MAX);
  
  for(unsigned int l=0; l<N_steps; ++l) {
    unsigned int i = rand()%N_nodes;
    
    double sum = 0; //!!!!! should be set to h_i !!!!!
    for(std::list<std::list<link>::iterator>::iterator it_it_links = nodes[i].attached_links.begin(); it_it_links != find_node(i)->attached_links.end(); ++it_it_links)
      if((*it_it_links)->node1 == nodes + i)
        sum += (*it_it_links)->weight * std::stoi((*it_it_links)->node2->label.name);
      else
        sum += (*it_it_links)->weight * std::stoi((*it_it_links)->node1->label.name);

    if( /*(1+tanh(sum/T))/2*/ 1/(1+exp(-2*sum/T)) > (rand()/RAND_MAX_DOUBLE) )
      nodes[i].set_spin_to_plus1();
    else
      nodes[i].set_spin_to_minus1();
  }
#ifndef QUIET_MODE
  std::cerr <<"* " << N_steps <<" steps of Glauber Dynamics were performed.\n";
#endif
  return *this;
}

const graph& graph::simulate_Glauber_Dynamics_with_fields(const int &seed, const double &T, const unsigned int &N_steps) {
#ifndef QUIET_MODE
  std::cerr << "*** Simulating Glauber Dynamics...\n";
#endif
  //////////////// CONSISTENCY CHECK //////////////
  if(N_labels > 2) {
    std::cerr << "---------------------------------------------------------------------------------------\n"
              << "ERROR: GRAPH IS INSUFFICIENT FOR GLAUBER DYNAMICS (has more than 2 different labels)!\n"
              << "---------------------------------------------------------------------------------------\n";
    return *this;
  };
  
  for(unsigned int i=0; i<N_nodes; ++i)
    if(nodes[i].label.name != "1" && nodes[i].label.name != "-1") {
      std::cerr << "-------------------------------------------------------------------------------------------------------\n"
                << "ERROR: GRAPH IS INSUFFICIENT FOR GLAUBER DYNAMICS (has labels that are different from \"1\" and \"-1\")!\n"
                << "DEFECTED LABEL OF NODE " << nodes[i].name << " (id="<< i <<") is " << nodes[i].label.name << '\n'
                << "-------------------------------------------------------------------------------------------------------\n";
      return *this;
  };

  ////////////// SIMULATING DYNAMICS ////////////////
  srand(seed);
  double RAND_MAX_DOUBLE = double(RAND_MAX);

  for(unsigned int l=0; l<N_steps; ++l) {
    unsigned int i = rand()%N_nodes; // Picking a node at random
    if(nodes[i].name != "FIELD_NODE") {
      double sum = 0;
      for(std::list<std::list<link>::iterator>::iterator it_it_links = nodes[i].attached_links.begin(); it_it_links != find_node(i)->attached_links.end(); ++it_it_links)
        if((*it_it_links)->node1 == nodes + i)
          sum += (*it_it_links)->weight * std::stoi((*it_it_links)->node2->label.name);
        else
          sum += (*it_it_links)->weight * std::stoi((*it_it_links)->node1->label.name);

      if( /*(1+tanh(sum/T))/2*/ 1/(1+exp(-2*sum/T)) > (rand()/RAND_MAX_DOUBLE) )
        nodes[i].set_spin_to_plus1();
      else
        nodes[i].set_spin_to_minus1();
    }
    else
      --l;
  }
#ifndef QUIET_MODE
  std::cerr <<"* " << N_steps <<" steps of Glauber Dynamics were performed.\n";
#endif
  
  return *this;
}

double graph::average_magnetization() const {
  ////////////// CONSISTENCY CHECK //////////////
  if(N_labels > 2) {
    std::cerr << "----------------------------------------------------------------------------------\n"
              << "ERROR: GRAPH IS INAPPROPRIATE SPIN MODEL (has more than 2 different labels)!\n"
              << "----------------------------------------------------------------------------------\n";
#ifdef INTERPRETER_MODE
    return -1;
#else
    exit(1);
#endif
  };
  
  for(unsigned int i=0; i<N_nodes; ++i)
    if(nodes[i].label.name != "1" && nodes[i].label.name != "-1") {
      std::cerr << "----------------------------------------------------------------------------------------------\n"
                << "ERROR: GRAPH IS INAPPROPRIATE SPIN MODEL (has labels that are different from \"1\" and \"-1\")!\n"
                << "DEFECTED LABEL OF NODE " << nodes[i].name << " (id="<< i <<") is " << nodes[i].label.name << '\n'
                << "----------------------------------------------------------------------------------------------\n";
#ifdef INTERPRETER_MODE
    return -1;
#else
    exit(2);
#endif
  };

  bool FIELD_NODE_DETECTED = false;
  double m_avrg = 0;
  for(unsigned int i=0; i<N_nodes; i++)
    if(nodes[i].name != "FIELD_NODE")
      m_avrg += std::stoi(nodes[i].label.name);
    else FIELD_NODE_DETECTED = true;

  if(FIELD_NODE_DETECTED)
    return m_avrg / (N_nodes-1);
  else
    return m_avrg / N_nodes;
}

graph& graph::set_spin_sequence(const col_vector<short> &V) {
  // std::cerr << "V.get_dim_y = "<< V.get_dim_y() << '\n';
  if(V.get_dim_y() == N_nodes) {
    for(unsigned int i=0; i<N_nodes; ++i)
        if(V[i] == 1)
          nodes[i].set_spin_to_plus1();
        else
          nodes[i].set_spin_to_minus1();
    return *this;
  }
  else {
    std::cerr << "------------------------------------------------------------------------------------------------\n"
              << "ERROR: Attempt to set 2-D spin sequence with inconsistent dimensions (dim_x * dim_y != N_nodes)!\n"
              << "------------------------------------------------------------------------------------------------\n";
#ifdef INTERPRETER_MODE
    return *this;
#else
    exit(1);
#endif
  }
}

graph& graph::set_spin_sequence(const row_vector<short> &V) {
  if(V.get_dim_x() == N_nodes) {
    for(unsigned int i=0; i<N_nodes; ++i)
        if(V[i] == 1)
          nodes[i].set_spin_to_plus1();
        else
          nodes[i].set_spin_to_minus1();
    return *this;
  }
  else {
    std::cerr << "------------------------------------------------------------------------------------------------\n"
              << "ERROR: Attempt to set 2-D spin sequence with inconsistent dimensions (dim_x * dim_y != N_nodes)!\n"
              << "------------------------------------------------------------------------------------------------\n";
#ifdef INTERPRETER_MODE
    return *this;
#else
    exit(1);
#endif
  }
}

graph& graph::set_spin_sequence_2D(const matrix<short> &M) {
  if(M.get_dim_x() * M.get_dim_y() == N_nodes) {
    for(unsigned int i=0; i<M.get_dim_y(); ++i)
      for(unsigned int j=0; j<M.get_dim_x(); ++j)
        if(M[i][j] == 1)
          nodes[i*M.get_dim_x() + j].set_spin_to_plus1();
        else
          nodes[i*M.get_dim_x() + j].set_spin_to_minus1();
    
    return *this;
  }
  else {
    std::cerr << "N_nodes = " << N_nodes << '\n'
              << "M.get_dim_x()= " << M.get_dim_x() << '\n'
              << "M.get_dim_y()= " << M.get_dim_y() << '\n';
    
    std::cerr << "------------------------------------------------------------------------------------------------\n"
              << "ERROR: Attempt to set 2-D spin sequence with inconsistent dimensions (dim_x * dim_y != N_nodes)!\n"
              << "------------------------------------------------------------------------------------------------\n";
#ifdef INTERPRETER_MODE
    return *this;
#else
    exit(1);
#endif
  }
}

matrix<short> graph::spin_sequence_2D(const unsigned int &dim_y, const unsigned int &dim_x) const {
  if(dim_x * dim_y == N_nodes) {
    matrix<short> M(dim_y, dim_x);
    
    for(unsigned int i=0; i<dim_y; ++i)
      for(unsigned int j=0; j<dim_x; ++j)
        M[i][j] = std::stoi(nodes[i*dim_x + j].label.name);
    return M;
  }
  else {
    std::cerr << "-----------------------------------------------------------------------------------------------------\n"
              << "ERROR: Attempt to arrange spins in 2-D array with inconsistent dimensions (dim_x * dim_y != N_nodes)!\n"
              << "-----------------------------------------------------------------------------------------------------\n";
#ifdef INTERPRETER_MODE
    return matrix<short>();
#else
    exit(1);
#endif
  }
}

////////// MULTIPARTITES BEGIN ////////////
const graph& graph::construct_multipartite_spin_model(const col_vector<unsigned int>& Np, const matrix<double>& J, const col_vector<double>& F) {
  if(J.get_dim_x() != J.get_dim_y()) {
    std::cerr << "----------------------------------------------------------------------------------\n"
              << "ERROR: Multipartite cannot be created. Reason:\n"
              << "       Interaction matrix is not square!\n"
              << "----------------------------------------------------------------------------------\n";
#ifdef INTERPRETER_MODE
    return -1;
#else
    exit(1);
#endif
  }

  for(unsigned int i=0; i<J.get_dim_x(); ++i)
    for(unsigned int j=0; j<i; ++j)
      if(J[i][j] != J[j][i]) {
        std::cerr << "----------------------------------------------------------------------------------\n"
                  << "ERROR: Multipartite cannot be created. Reason:\n"
                  << "       Interaction matrix is not symmetric!\n"
                  << "----------------------------------------------------------------------------------\n";
#ifdef INTERPRETER_MODE
        return -1;
#else
        exit(1);
#endif
      }
  
  if(J.get_dim_x() != F.size()) {
    std::cerr << "----------------------------------------------------------------------------------\n"
              << "ERROR: Multipartite cannot be created. Reason:\n"
              << "       Dimensions of Interaction matrix " <<  J.get_dim_x() << " are inconsistent with dim of Fields vector " << F.size() << '\n'
              << "       They should be equal!\n"
              << "----------------------------------------------------------------------------------\n";
#ifdef INTERPRETER_MODE
    return -1;
#else
    exit(1);
#endif
  }

  ////////////// ADDING NODES TO THE GRAPH /////////////
  srand(time(NULL));
  for(unsigned int i=0; i<Np.size(); ++i)
    for(unsigned int j=0; j<Np[i]; ++j)
      add_node("PART"+std::to_string(i)+'_'+std::to_string(j), std::to_string(1-2*(rand()%2)) );;
  
  ///////////////// LINKING THE GRAPH //////////////////
  col_vector<unsigned int> Bounds(Np.size()+1);
  for(unsigned int i=0; i<Bounds.size(); ++i)
    Bounds[i] = 0;
  for(unsigned int i=1; i<Bounds.size(); ++i)
    for(unsigned int j=0; j<i; ++j)
      Bounds[i] += Np[j]; // Bounds == [0, N1, N1+N2, N1+N2+N3, N1+N2+N3+N4, ...]

  //Linking diagonal
  for(unsigned int l=1; l<Bounds.size(); ++l)
    if(J[l-1][l-1]!=0)
      for(unsigned int i=Bounds[l-1]; i<Bounds[l]; ++i)
        for(unsigned int j=Bounds[l-1]; j<i; ++j)
          add_link_no_checks(i, j, 'P'+std::to_string(l)+"_P"+std::to_string(l), J[l-1][l-1]);
  // Linking off diagonal
  for(unsigned int l=1; l<Bounds.size(); ++l)
    for(unsigned int m=1; m<l; ++m)
      if(J[l-1][m-1]!=0)
        for(unsigned int i=Bounds[l-1]; i<Bounds[l]; ++i)
          for(unsigned int j=Bounds[m-1]; j<Bounds[m]; ++j)
            add_link_no_checks(i, j, 'P'+std::to_string(m)+"_P"+std::to_string(l), J[l-1][m-1]);

  /////// ADDING FIELDS AS FICTIVE NODE AND LINKS //////
  add_node("FIELD_NODE", "1"); //This node is added last and has ID == N_nodes-1
  for(unsigned int l=1; l<Bounds.size(); l++)
    if(F[l-1]!=0)
      for(unsigned int i=Bounds[l-1]; i<Bounds[l]; ++i)
        add_link(N_nodes-1, i, "FIELD" + std::to_string(l), F[l-1]);
  
  return *this;
}

const graph& graph::update_multipartite_spin_model(const col_vector<unsigned int>& Np, const matrix<double>& J, const col_vector<double>& F) {
  clear_links();
  ///////////////// LINKING THE GRAPH //////////////////
  col_vector<unsigned int> Bounds(Np.size()+1);
  for(unsigned int i=0; i<Bounds.size(); ++i)
    Bounds[i] = 0;
  for(unsigned int i=1; i<Bounds.size(); ++i)
    for(unsigned int j=0; j<i; ++j)
      Bounds[i] += Np[j]; // Bounds == [0, N1, N1+N2, N1+N2+N3, N1+N2+N3+N4, ...]

  //Linking diagonal
  for(unsigned int l=1; l<Bounds.size(); ++l)
    if(J[l-1][l-1]!=0)
      for(unsigned int i=Bounds[l-1]; i<Bounds[l]; ++i)
        for(unsigned int j=Bounds[l-1]; j<i; ++j)
          add_link_no_checks(i, j, 'P'+std::to_string(l)+"_P"+std::to_string(l), J[l-1][l-1]);
  // Linking off diagonal
  for(unsigned int l=1; l<Bounds.size(); ++l)
    for(unsigned int m=1; m<l; ++m)
      if(J[l-1][m-1]!=0)
        for(unsigned int i=Bounds[l-1]; i<Bounds[l]; ++i)
          for(unsigned int j=Bounds[m-1]; j<Bounds[m]; ++j)
            add_link_no_checks(i, j, 'P'+std::to_string(m)+"_P"+std::to_string(l), J[l-1][m-1]);

  /////// ADDING FIELDS AS FICTIVE LINKS TO BIAS NODE ///////
  for(unsigned int l=1; l<Bounds.size(); l++)
    if(F[l-1]!=0)
      for(unsigned int i=Bounds[l-1]; i<Bounds[l]; ++i)
        add_link(N_nodes-1, i, "FIELD" + std::to_string(l), F[l-1]);

  return *this;
}

bool graph::check_multipartite(col_vector<unsigned int>& Np) const {
  ////////////// CONSISTENCY CHECK //////////////
  if(N_labels > 2) {
    std::cerr << "----------------------------------------------------------------------------------\n"
              << "ERROR: GRAPH IS INAPPROPRIATE SPIN MODEL (has more than 2 different labels)!\n"
              << "----------------------------------------------------------------------------------\n";
#ifdef INTERPRETER_MODE
    return -1;
#else
    exit(1);
#endif
  };
  
  for(unsigned int i=0; i<N_nodes; ++i)
    if(nodes[i].label.name != "1" && nodes[i].label.name != "-1") {
      std::cerr << "----------------------------------------------------------------------------------------------\n"
                << "ERROR: GRAPH IS INAPPROPRIATE SPIN MODEL (has labels that are different from \"1\" and \"-1\")!\n"
                << "DEFECTED LABEL OF NODE " << nodes[i].name << " (id="<< i <<") is " << nodes[i].label.name << '\n'
                << "----------------------------------------------------------------------------------------------\n";
#ifdef INTERPRETER_MODE
    return -1;
#else
    exit(2);
#endif
  };

  ////////////// NODES IDS CHECK //////////////
  unsigned int count = 0;
  for(unsigned int l=0; l<Np.size(); ++l) {
    std::size_t pos_ = nodes[count].name.find('_') + 1;
    std::string partite_name = nodes[count].name.substr(0, pos_);
    for(unsigned int i=count; i<count+Np[l]; ++i)
      if(nodes[i].name.substr(0, pos_) != partite_name) {
        std::cerr << "--------------------------------------------------------------\n"
                  << "ERROR: Incorrect assignment of nodes IDs for " << partite_name << '\n'
                  << "The problem with the node[" << i << "], which is named " << nodes[i].name << '\n'
                  << "--------------------------------------------------------------\n";
#ifdef INTERPRETER_MODE
        return -1;
#else
        exit(2);
#endif
      }
    count += Np[l];
  }

  ////////////// NUMBER ON NODES CHECK //////////////
  unsigned int N = 0; //Number of nodes based on the vector provided
  for(unsigned int i=0; i<Np.size(); ++i)
    N += Np[i];
  if(N+1 != N_nodes) {
    if(N == N_nodes)
      std::cerr << "----------------------------------------------------------------------------------\n"
                << "WARNING: There is no feeld node in the graph!\n"
                << "----------------------------------------------------------------------------------\n";
    else {
      std::cerr << "----------------------------------------------------------------------------------\n"
                << "ERROR: Number of nodes in the graph is " << N_nodes << " instead of " << N+1 << '\n'
                << "----------------------------------------------------------------------------------\n";
#ifdef INTERPRETER_MODE
      return -1;
#else
      exit(1);
#endif
    }
  }

  ////////////// FIELD NODE CHECK //////////////
  if(nodes[N].name != "FIELD_NODE") {
    std::cerr << "----------------------------------------------------------------------------------\n"
              << "ERROR: The field node is named " << nodes[N].name << " instead of FIELD_NODE\n"
              << "----------------------------------------------------------------------------------\n";
#ifdef INTERPRETER_MODE
    return -1;
#else
    exit(1);
#endif
  }
  return 1;
}


col_vector<double> graph::partites_magnetizations(col_vector<unsigned int>& Np) const {
  col_vector<double> M(Np.size());

  double m;
  unsigned int count = 0;
  for(unsigned int l=0; l<Np.size(); ++l) {
    m = 0;
    for(unsigned int i=count; i<count+Np[l]; i++) {
      m += std::stoi(nodes[i].label.name);
    }
    M[l] = m;
    count += Np[l];
  }

  return M;
}

col_vector<double> graph::partites_average_magnetizations(col_vector<unsigned int>& Np) const {
  col_vector<double> M(Np.size());

  double m;
  unsigned int count = 0;
  for(unsigned int l=0; l<Np.size(); ++l) {
    m = 0;
    for(unsigned int i=count; i<count+Np[l]; i++) {
      m += std::stoi(nodes[i].label.name);
    }
    M[l] = m/Np[l];
    count += Np[l];
  }
  
  return M;
}

////////// MULTIPARTITES END ////////////

///////////////////// SPIN DYNAMICS END ///////////////////////////////////

///////////////////// CLASSIFICATION BEGIN ////////////////////////////////
labeling graph::initialize_labeling_for_classification(const unsigned short &N_classes, const std::string* class_names) {
  if(N_nodes < N_classes) {
    std::cerr << "--------------------------------------------------------------------------\n"
              << "ERROR: Number of nodes in the graph is smaller than the number of classes!\n"
              << "--------------------------------------------------------------------------\n";
#ifdef INTERPRETER_MODE
    return -1;
#else
    exit(1);
#endif
  }

  labeling L(N_nodes, N_classes);
  unsigned int i=0;
  do L.set_label_by_index(i, i);
  while(++i < N_classes);
  do L.set_label_by_index(i, N_classes-1);
  while(++i < N_nodes);
  
  load_labeling(L, class_names);
  
  return L;
}

graph& graph::load_labeling(const labeling &L, const std::string* class_names) {
  if(N_labels != L.num_classes()) {
    delete[] labels;
    N_labels = L.num_classes();
    labels = new label[N_labels];
  }
  if(class_names) {
    for(unsigned short i=0; i<N_labels; ++i) {
      labels[i].name = class_names[i];
      labels[i].id = i;
    }
  }
  else {
    for(unsigned short i=0; i<N_labels; ++i) {
      // labels[i].name = std::to_string(i);
      labels[i].id = i;
    }
  }
  for(unsigned int i=0; i<N_nodes; ++i) {
    nodes[i].label.id = L[i]; //Label IDs only create overheads even here
    nodes[i].label.name = labels[L[i]].name;
  }
  return *this;
}

matrix<double> graph::link_probabilities() const {//THIS FUNCTION USES LABELS IDs!
  matrix<double> L(N_labels, N_labels);
  for(unsigned int i=0; i<N_labels; ++i)
    L[i][i] = 0;
  for(unsigned int i=0; i<N_labels; ++i)
    for(unsigned int j=0; j<i; j++)
      L[i][j] = L[j][i] = 0;
  
  for(std::list<link>::const_iterator it_links=links.begin(); it_links!=links.end(); ++it_links)
    ++L[it_links->node1->label.id][it_links->node2->label.id];

  for(unsigned int i=0; i<N_labels; ++i)
    for(unsigned int j=0; j<i; j++)
      L[j][i] = (L[i][j] += L[j][i]);
  // Now L[i][j] contains the number of links connecting nodes with labels i and j

  unsigned int *node_counters = new unsigned int[N_labels]; //Count nodes with certain labels
  for(unsigned int l=0; l<N_labels; ++l) {
    node_counters[l] = 0;
    for(unsigned int i=0; i<N_nodes; ++i)
      if(nodes[i].label.id == l)
        ++node_counters[l];
  }
  // Now we know the counts of nodes of each type.
  // We will divide numbers of links in L by total possible numbers of links

  for(unsigned int i=0; i<N_labels; ++i)
    for(unsigned int j=0; j<i; j++)
      L[i][j] = L[j][i] /= node_counters[i]*node_counters[j];
  for(unsigned int i=0; i<N_labels; ++i)
    L[i][i] /= node_counters[i]*(node_counters[i]-1)/2;
  
  delete[] node_counters;
  return L;
}

matrix<double> graph::regularized_link_probabilities() const {
  matrix<double> L(N_labels, N_labels);
  for(unsigned int i=0; i<N_labels; ++i)
    L[i][i] = 0;
  for(unsigned int i=0; i<N_labels; ++i)
    for(unsigned int j=0; j<i; j++)
      L[i][j] = L[j][i] = 0;
  
  for(std::list<link>::const_iterator it_links=links.begin(); it_links!=links.end(); ++it_links)
    ++L[it_links->node1->label.id][it_links->node2->label.id];

  for(unsigned int i=0; i<N_labels; ++i)
    for(unsigned int j=0; j<i; j++)
      L[j][i] = (L[i][j] += L[j][i] + 0.0001);
  for(unsigned int i=0; i<N_labels; ++i)
    L[i][i] += 0.0001;
  // Now L[i][j] contains the number of links connecting nodes with labels i and j

  unsigned int *node_counters = new unsigned int[N_labels]; //Count nodes with certain labels
  for(unsigned int l=0; l<N_labels; ++l) {
    node_counters[l] = 0;
    for(unsigned int i=0; i<N_nodes; ++i)
      if(nodes[i].label.id == l)
        ++node_counters[l];
  }
  // Now we know the counts of nodes of each type.
  // We will divide numbers of links in L by total possible numbers of links

  for(unsigned int i=0; i<N_labels; ++i)
    for(unsigned int j=0; j<i; j++)
      L[i][j] = L[j][i] /= node_counters[i]*node_counters[j] + 0.01;
  for(unsigned int i=0; i<N_labels; ++i)
    L[i][i] /= node_counters[i]*(node_counters[i]-1)/2 + 0.01;
  
  delete[] node_counters;
  return L;
}

double graph::loglikelihood() const {
// This function deliberately does not use regularized_link_probabilities() funaction for optimisation
  matrix<unsigned int> L(N_labels, N_labels);
  for(unsigned int i=0; i<N_labels; ++i)
    L[i][i] = 0;
  for(unsigned int i=0; i<N_labels; ++i)
    for(unsigned int j=0; j<i; j++)
      L[i][j] = L[j][i] = 0;
  
  for(std::list<link>::const_iterator it_links=links.begin(); it_links!=links.end(); ++it_links)
    ++L[it_links->node1->label.id][it_links->node2->label.id];

  for(unsigned int i=0; i<N_labels; ++i)
    for(unsigned int j=0; j<i; j++)
      L[j][i] = (L[i][j] += L[j][i]);
  // Now L[i][j] contains the number of links connecting nodes with labels i and j
  unsigned int *node_counters = new unsigned int[N_labels]; //Count nodes with certain labels
  for(unsigned int l=0; l<N_labels; ++l) {
    node_counters[l] = 0;
    for(unsigned int i=0; i<N_nodes; ++i)
      if(nodes[i].label.id == l)
        ++node_counters[l];
  }
  // Now we know the counts of nodes of each type and can compute maximim numbers of links
  
  double log_likelihood = 0;
  for(unsigned int i=0; i<N_labels; ++i)
    for(unsigned int j=0; j<i; ++j) {
      double L_max = node_counters[i]*node_counters[j];
      double link_prob = L[i][j]/L_max;
      log_likelihood += L[i][j]*log(link_prob) + (L_max-L[i][j])*log(1-link_prob);
    }
  for(unsigned int i=0; i<N_labels; ++i) {
    double L_max = node_counters[i]*(node_counters[i]-1)/2;
    double link_prob = L[i][i]/L_max;
    log_likelihood += L[i][i]*log(link_prob) + (L_max-L[i][i])*log(1-link_prob);
  }
  delete[] node_counters;
  return log_likelihood;
}

double graph::regularized_loglikelihood() const {
  // This function deliberately does not use regularized_link_probabilities() funaction for optimisation
  matrix<unsigned int> L(N_labels, N_labels);
  for(unsigned int i=0; i<N_labels; ++i)
    L[i][i] = 0;
  for(unsigned int i=0; i<N_labels; ++i)
    for(unsigned int j=0; j<i; j++)
      L[i][j] = L[j][i] = 0;
  
  for(std::list<link>::const_iterator it_links=links.begin(); it_links!=links.end(); ++it_links)
    ++L[it_links->node1->label.id][it_links->node2->label.id];

  for(unsigned int i=0; i<N_labels; ++i)
    for(unsigned int j=0; j<i; j++)
      L[j][i] = (L[i][j] += L[j][i]);
  // Now L[i][j] contains the number of links connecting nodes with labels i and j
  unsigned int *node_counters = new unsigned int[N_labels]; //Count nodes with certain labels
  for(unsigned int l=0; l<N_labels; ++l) {
    node_counters[l] = 0;
    for(unsigned int i=0; i<N_nodes; ++i)
      if(nodes[i].label.id == l)
        ++node_counters[l];
  }
  // Now we know the counts of nodes of each type and can compute maximim numbers of links
  
  double log_likelihood = 0;
  for(unsigned int i=0; i<N_labels; ++i)
    for(unsigned int j=0; j<i; ++j) {
      double L_max = node_counters[i]*node_counters[j];
      double link_prob = (L[i][j] + 0.0001)/(L_max + 0.01);
      log_likelihood += L[i][j]*log(link_prob) + (L_max-L[i][j])*log(1-link_prob);
    }
  for(unsigned int i=0; i<N_labels; ++i) {
    double L_max = node_counters[i]*(node_counters[i]-1)/2;
    double link_prob = (L[i][i] + 0.0001)/(L_max + 0.01);
    log_likelihood += L[i][i]*log(link_prob) + (L_max-L[i][i])*log(1-link_prob);
  }
  delete[] node_counters;
  return log_likelihood;
}

std::list<labeling> graph::modularity_classifier_precise(const unsigned short &num_classes, const std::string* class_names) {
  std::list<labeling> labelings;
  labeling L = initialize_labeling_for_classification(num_classes, class_names);
  labelings.push_back(L);
  double modularity, modularity_max = -1.79769e+308; //Initializing with min(double)
  do {
    load_labeling(L);
    modularity = this->modularity();
    if(modularity > modularity_max) {
      labelings.clear();
      labelings.push_back(L);
      modularity_max = modularity;
    }
    else 
      if(modularity == modularity_max)
        labelings.push_back(L);
  } while(!--L);
#ifndef QUIET_MODE
  std::cerr << "MAX_MODULARITY = " <<  modularity_max << '\n';
  std::cerr << "MAX_RELATIVE_MODULARITY = " << this->load_labeling(labelings.back()).relative_modularity()) << '\n';
#endif
  return labelings;
}

std::list<labeling> graph::ML_classifier_precise(const unsigned short &num_classes, const std::string* class_names) {
  std::list<labeling> labelings;
  labeling L = initialize_labeling_for_classification(num_classes, class_names);
  double loglikelihood, loglikelihood_max = -1.79769e+308; //Initializing with min(double)
  do {
    load_labeling(L);
    loglikelihood = regularized_loglikelihood();
    if(loglikelihood > loglikelihood_max) {
      labelings.clear();
      labelings.push_back(L);
      loglikelihood_max = loglikelihood;
    }
    else 
      if(loglikelihood == loglikelihood_max)
        labelings.push_back(L);  
  } while(!--L);
#ifdef QUIET_MODE
  std::cerr << "MAX_LOGLIKELIHOOD = " << loglikelihood_max << '\n';
#endif
  return labelings;
}
////////////////////// CLASSIFICATION END /////////////////////////////////

///////////////////////////////////////////////////////////////////////////
///////////////// GRAPH COMPUTATIONAL FUNCTIONS END ///////////////////////
///////////////////////////////////////////////////////////////////////////

graph& graph::clear_links() {
  std::list<link>::iterator it_links;
  while(!links.empty()) {
    it_links = links.begin();
    remove_link(it_links);
  }  
  return *this;
}

graph& graph::clear() {
  delete[] nodes;
  nodes = NULL;
  delete[] labels;
  labels = NULL;
  links.clear();
    
  N_nodes = 0;
  N_links = 0;
  N_labels = 0;

  NODES_ARRAY_SIZE = 0;

  return *this;
}

graph& graph::operator=(const graph& gr) {
  clear();

  /*NODES_SORTED = gr.NODES_SORTED;*/
  NODES_ARRAY_SIZE = gr.NODES_ARRAY_SIZE;
  nodes = new node[NODES_ARRAY_SIZE];

  N_nodes = gr.N_nodes;
  N_labels = gr.N_labels;
  N_links = gr.N_links;
  
  for(unsigned int i=0; i<N_nodes; ++i)
    nodes[i] = gr.nodes[i];

  labels = new label[N_labels];
  for(unsigned int i=0; i<N_labels; ++i)
    labels[i] = gr.labels[i];

  // Copying links, redirecting node pointers
  for(std::list<link>::const_iterator it_links = gr.links.begin(); it_links!=gr.links.end(); ++it_links) {
    links.push_back(*it_links);
    links.back().node1 = nodes + (it_links->node1 - gr.nodes);
    links.back().node2 = nodes + (it_links->node2 - gr.nodes);
  }

  return *this;
}
/////////////////////////////////////////////////////////////////////////
void graph::append_nodes_array(const unsigned int& n) {
  NODES_ARRAY_SIZE+=n;    
  node* nodes_ = new node[NODES_ARRAY_SIZE];
  for(unsigned int i=0; i<N_nodes; i++) {
    nodes_[i] = nodes[i];
    for(std::list<std::list<link>::iterator>::iterator it_it_links = nodes_[i].attached_links.begin(); it_it_links!=nodes_[i].attached_links.end(); ++it_it_links)
      if((*it_it_links)->node1 == nodes+i)
        (*it_it_links)->node1 = nodes_ + i;
      else
        (*it_it_links)->node2 = nodes_ + i;
  }
  delete[] nodes;
  nodes = nodes_;
}
/////////////////////////////////////////////////////////////////////////
graph::~graph() {
  delete[] nodes;
  delete[] labels;
}
/////////////////////////////////////////////////////////////////////////
std::ostream& operator<<(std::ostream &os, const graph &gr) {
  os << "-------------------------GRAPH BEGINNING------------------------------\n";
  os << "PARAMETERS:\n";
  os << "* N_nodes = " << gr.N_nodes << '\n';
  os << "* N_links = " << gr.N_links << '\n';
  os << "* N_labels = " << gr.N_labels << '\n';

  os << "NODES:\n";
  for(unsigned int i=0; i<gr.N_nodes; i++) {
    os << "* " << gr.nodes[i] << '\n';
  }

  os << "LINKS:\n";
  for(std::list<link>::const_iterator it_links = gr.links.begin(); it_links != gr.links.end() ; ++it_links)
    os << "* " << *it_links << '\n';

  os << "LABELS:\n";
  for(unsigned int i=0; i<gr.N_labels; i++)
    os << "* " << gr.labels[i].name << '\n';
  
  os << "----------------------------GRAPH END---------------------------------\n";
    
  return os;
}

//////////////////////////////////////////////////////////////////
/////////////////////// graph::label /////////////////////////////
//////////////////////////////////////////////////////////////////

graph::label::label() : id(-1) {}

graph::label::label(const std::string &name_, unsigned int id_) : id(id_), name(name_) {}

graph::label::label(const unsigned int &id_, const std::string &name_) : id(id_), name(name_) {}

graph::label::label(const label &l) : id(l.id), name(l.name) {}

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

//////////////////////////////////////////////////////////////////
//////////////////////////// link ////////////////////////////////
//////////////////////////////////////////////////////////////////

link::link(const link &l) : node1(l.node1), node2(l.node2), type(l.type), weight(l.weight) {}

link::link(node* p_nd1, node* p_nd2, const std::string &type_, double weight_) : node1(p_nd1), node2(p_nd2), type(type_), weight(weight_) {}
///////////////////////////////////////////////////////////////
link& link::operator=(const link &l) {
  node1 = l.node1;
  node2 = l.node2;
  type = l.type;
  weight = l.weight;
  
  return *this;
}
///////////////////////////////////////////////////////////////
bool link::operator==(const link &l) {
  if((node1!=l.node1 && node1!=l.node2) || (node2!=l.node2 && node2!=l.node1) || weight!=l.weight || type != l.type)
    return false;
  
  return true;
}

bool link::operator!=(const link &l) {
  return !(*this == l);
}
///////////////////////////////////////////////////////////////
std::ostream& operator<<(std::ostream &os, const link &l) {
  os << l.node1->name << ' ' << l.type << ' ' << l.node2->name;
  return os;
}

//////////////////////////////////////////////////////////////////
//////////////////////////// node ////////////////////////////////
//////////////////////////////////////////////////////////////////

node::node() {}

node::node(const std::string &name_) : name(name_) {}

node::node(const unsigned int &id_, const std::string &name_, const unsigned int &label_id, const std::string &label_name) : name(name_), label(label_id,label_name), id(id_)  {}

node::node(const std::string &name_, const std::string &label_) :  name(name_), label(label_) {}

node::node(const node &nd) : name(nd.name), label(nd.label), id(nd.id) {}

//////////////////////////////////////////////////////////////////
inline void node::set_spin_to_minus1() {
  if(label.name != "-1") {
    label.name = "-1";
    label.id = !label.id;
  }
}

inline void node::set_spin_to_plus1() {
 if(label.name != "1") {
    label.name = "1";
    label.id = !label.id;
  } 
}
//////////////////////////////////////////////////////////////////

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
//////////////////////////////////////////////////////////////////
node& node::operator=(const node &nd) {
  name = nd.name;
  label = nd.label;
  id = nd.id;
  attached_links = nd.attached_links;

  return *this;
}
//////////////////////////////////////////////////////////////////

std::ostream& operator<<(std::ostream &os, const node &nd) {
  os << nd.name << ' ' << nd.label.name << " degree = " << nd.attached_links.size();
  return os;
}
