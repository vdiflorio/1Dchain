#ifndef PARAM_H
#define PARAM_H

#include <unordered_map>
#include <iostream>
#include <string>

struct
params_t {

  std::unordered_map<std::string, int> iparams;
  std::unordered_map<std::string, double> dparams;
  std::unordered_map<std::string, std::string> sparams;
  std::unordered_map<std::string, bool> bparams;

  void
  read (const std::string &fname);

  void
  write (const std::string &fname);
};

std::ostream&
operator<< (std::ostream& os, const params_t& p);

std::istream&
operator>> (std::istream& is, params_t& p);

extern params_t p;

#endif
