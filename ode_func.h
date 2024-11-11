#include <iostream>
#include <cmath>
#include <iomanip>
#include <fstream>
#include <vector>
#include <algorithm>



using namespace std;



#define  chi     1.0  //"intensità" potenziale armonico 

#define  Alpha   1.0 //"intensità" potenziale anarmonico cubico
#define  N       20   //# di masse mobili
#define  dim     1    //spazio dimensioni di oscillazione
#define  thetaL  1.0  //tempo di rilassamento termost sx
#define  thetaR  1.0  //tempo di rilassamento termost dx
#define  Tl      1.0  //Temp termostato sx
#define  Tr      1.0  //Temp termostato dx
#define  Kb      1.0  //Costante di Boltzmann
#define  res     1.0

const int  a     =  1.0;  //distanza di eq. tra le part.

const int  m      = 1.0;
const int bet   = 1.0;  //"intensità" potenziale anarmonico  quartico
void Chain (double t, double *Y, double *R);   //our FPUT

void Chain1 (double t, double *Y, double *R);  //real FPUT

void SPC (double t, double *Y, double *R);

void LJ (double t, double *Y, double *R);

void AlfaBeta (double t, double *Y, double *R);

// void generateRandomIntegers(int min_val, int max_val, std::vector<int>* random_numbers);



