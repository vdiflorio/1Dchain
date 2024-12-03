#include <fstream>
#include <iomanip>
#include <nlohmann/json.hpp>

#include "param.h"

params_t p{
  {
    {"dim", 1},
    {"N", 10},
    {"step", 50000},
    {"catene_CI", 10},
    {"catene_scelte", 100}
  },
  {
    {"dt", 5e-4},
    {"m", 1.0},
    {"a", 1.0},
    {"chi", 1.0},
    {"alpha", 0.0},
    {"beta", 0.5},
    {"Tl", 1.0},
    {"Tr", 5.0},
    {"thetaL", 1.0},
    {"thetaR", 1.0},
    {"Kb", 1.0}
  },
  {
    {"ddata", "densita.dat"},
    {"tddata", "temp_dens.dat"}
  },
  {
    {"save_conditions", false}
  }
}; 

std::ostream&
operator<< (std::ostream& os, const params_t& p) {
  nlohmann::json my_json{{"iparams", p.iparams}, 
                         {"dparams", p.dparams}, 
                         {"sparams", p.sparams},
                         {"bparams", p.bparams}};
  os << std::setw(4) << my_json;
  return os;
};

std::istream&
operator>> (std::istream& is, params_t& p) {
  nlohmann::json my_json;
  is >> my_json;

  p.iparams = my_json["iparams"].get<std::unordered_map<std::string, int>> ();
  p.dparams = my_json["dparams"].get<std::unordered_map<std::string, double>> ();
  p.sparams = my_json["sparams"].get<std::unordered_map<std::string, std::string>> ();
  p.bparams = my_json["bparams"].get<std::unordered_map<std::string, bool>> ();
  return is;
}

void
params_t::read (const std::string &fname) {
  std::ifstream fparametri(fname.c_str ());
  fparametri >> (*this);
};

void
params_t::write (const std::string &fname) {
  std::ofstream fparametri(fname.c_str ());
  fparametri << (*this);
  fparametri.close ();
};
