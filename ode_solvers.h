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

#ifndef ODE_SOLVERS_H
#define ODE_SOLVERS_H

#include <iostream>
#include <cmath>
#include <iomanip>
#include <fstream>
#include <vector>
#include <functional>

void EulerStep (double t, std::vector<double> &Y,
                std::function<void (double, std::vector<double> &, std::vector<double> &)> dYdt,
                double dt, int neq);

void RK2Step (double t, std::vector<double> &Y,
              std::function<void (double, std::vector<double> &, std::vector<double> &)> dYdt,
              double dt, int neq);

void RK4Step (double t, std::vector<double> &Y,
              std::function<void (double, std::vector<double> &, std::vector<double> &)> dYdt,
              double dt, int neq);

void RK4Step_fast (double t, std::vector<double> &Y,
                   std::function<void (double, std::vector<double> &, std::vector<double> &)> dYdt,
                   double dt, int neq);

void RK5Step (double t, std::vector<double> &Y,
              std::function<void (double, std::vector<double> &, std::vector<double> &)> dYdt,
              double dt, int neq);

void PositionVerletStep (double *x, double *v, void (*acc) (double *,double *),
                         double dt, int neq);

#endif
