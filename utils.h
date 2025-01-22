#ifndef UTILS_H
#define UTILS_H
#include <iostream>
#include <cmath>


double omega_0 (std::vector<double> &Y, double T);

double TTCF (std::function<double (std::vector<double> &)> obs, double omega,std::vector<double> &Y,double T);

double observable (std::vector<double> &Y);

double dumb_observable (std::vector<double> &Y);

double observable_tot (std::vector<double> &Y);

double observable_bulk (std::vector<double> &Y);

void read_conditions (std::vector<double>& condizioni, int num_condizioni, int neq);

int save_condizioni_iniziali (int num_catene);


void compute_mean (int num_catene);

#endif
