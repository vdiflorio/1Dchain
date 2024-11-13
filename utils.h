#include <iostream>
#include <cmath>


double omega_0(std::vector<double> &Y, double T);

double TTCF(double (*obs)(std::vector<double> &), double omega ,std::vector<double> &Y,double T);

double observable(std::vector<double> &Y);

void read_conditions(std::vector<std::vector<double>>& condizioni, int num_condizioni, int neq);

int save_condizioni_iniziali(int num_catene);

