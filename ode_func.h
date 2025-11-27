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



