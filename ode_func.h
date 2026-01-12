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

#ifndef CONSTANTS_H
#define CONSTANTS_H

#include <iostream>
#include <cmath>
#include <iomanip>
#include <fstream>
#include <vector>
#include <algorithm>


void betaFPUT (double t, std::vector<double> &Y, std::vector<double> &R); //real FPUT

void LJ (double t, std::vector<double> &Y, std::vector<double> &R);

void AlfaBetaFPUT_initial (double t, std::vector<double> &Y, std::vector<double> &R);

void AlfaBetaFPUT (double t, std::vector<double> &Y, std::vector<double> &R);

void LepriChain (double t, std::vector<double> &Y, std::vector<double> &R);

void LepriChain_initial (double t, std::vector<double> &Y, std::vector<double> &R);

#endif // CONSTANTS_H



