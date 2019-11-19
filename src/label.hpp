#ifndef __LABEL_HPP__
#define __LABEL_HPP__

#include "graph.hpp"

struct graph::label {
  unsigned int id;
  std::string name;

  label() : id(-1) {}
  label(const std::string &name_, const unsigned int &id_=-1) : id(id_), name(name_) {}
  label(const unsigned int &id_, const std::string &name_) : id(id_), name(name_) {}
  label(const label &l) : id(l.id), name(l.name) {};

  label& operator=(const label&);
  label& operator=(const std::string&);

  bool operator==(const std::string&);
  bool operator==(const label&);
  bool operator<(const label&);
};

#endif
