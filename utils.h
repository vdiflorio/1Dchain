#include <iostream>
#include <cmath>


double omega_0(double *Y, double T);

double TTCF(double (*obs)(double *), double omega ,double *Y,double T);

double observable(double *Y);