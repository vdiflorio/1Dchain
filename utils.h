/*
 * This file is part of 1Dchain.
 *
 * Copyright (C) 2026
 * Vincenzo Di Florio, Davide Carbone
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 */

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

double observable_bulk_pinning (std::vector<double> &Y);

void read_conditions (std::vector<double>& condizioni, int num_condizioni, int neq);

int save_condizioni_iniziali (int num_catene);

void generate_condition (const std::vector<double>& base_cond,
                         std::vector<double>& new_cond,
                         int neq);

void read_conditions_subset (std::vector<double>& condizioni, int neq, const int max_catene, int job_id);

void compute_mean ( );

void timing_RK (int neq);

#endif
