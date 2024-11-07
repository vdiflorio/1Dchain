#include <iostream>
#include <cmath>
#include <iomanip>
#include <fstream>
using namespace std;

#define  m       1.0  //massa delle particelle
#define  a       1.0  //distanza di eq. tra le part.
#define  chi     1.0  //"intensità" potenziale armonico 
#define  beta    1.0  //"intensità" potenziale anarmonico  quartico
#define  Alpha   1.0 //"intensità" potenziale anarmonico cubico
#define  N       20   //# di masse mobili
#define  dim     1    //spazio dimensioni di oscillazione
#define  thetaL  1.0  //tempo di rilassamento termost sx
#define  thetaR  1.0  //tempo di rilassamento termost dx
#define  Tl      5.0  //Temp termostato sx
#define  Tr      1.0  //Temp termostato dx
#define  Kb      1.0  //Costante di Boltzmann
#define  res     1.0

void Chain (double t, double *Y, double *R);   //our FPUT

void Chain1 (double t, double *Y, double *R);  //real FPUT

void SPC (double t, double *Y, double *R);

void LJ (double t, double *Y, double *R);

void AlfaBeta (double t, double *Y, double *R);

