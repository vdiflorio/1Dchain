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

void generate_condition (const std::vector<double>& base_cond,
                         std::vector<double>& new_cond,
                         int neq);

void read_conditions_subset (std::vector<double>& condizioni, int neq, const int max_catene);

void compute_mean ( );

void timing_RK (int neq);

#endif
