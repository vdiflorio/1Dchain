#ifndef CONSTANTS_H
#define CONSTANTS_H

#include <iostream>
#include <cmath>
#include <iomanip>
#include <fstream>
#include <vector>
#include <algorithm>

#define  chi     1.0  //"intensità" potenziale armonico 

#define  Alpha   0.0 //"intensità" potenziale anarmonico cubico
#define  dim     1    //spazio dimensioni di oscillazione
#define  thetaL  1.0  //tempo di rilassamento termost sx
#define  thetaR  1.0  //tempo di rilassamento termost dx
#define  Kb      1.0  //Costante di Boltzmann
#define  res     1.0







// Declare global variables
extern int N;       // Number of particles
extern double Tl;   // Temperature of left thermostat
extern double Tr;   // Temperature of right thermostat
const int  a     =  1.0;  //distanza di eq. tra le part.

const int  m      = 1.0;
const int bet   = 1.0;  //"intensità" potenziale anarmonico  quartico

#endif // CONSTANTS_H
void Chain (double t, std::vector<double> &Y, std::vector<double> &R);   //our FPUT

void Chain1 (double t, std::vector<double> &Y, std::vector<double> &R);  //real FPUT

void SPC (double t, std::vector<double> &Y, std::vector<double> &R);

void LJ (double t, std::vector<double> &Y, std::vector<double> &R);

void AlfaBeta (double t, std::vector<double> &Y, std::vector<double> &R);

void AlfaBeta_corrected (double t, std::vector<double> &Y, std::vector<double> &R);

// void generateRandomIntegers(int min_val, int max_val, std::vector<int>* random_numbers);



