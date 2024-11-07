#include <iostream>
#include <cmath>
#include <iomanip>
#include <fstream>
using namespace std;

void EulerStep (double t, double *Y, void (*dYdt)(double, double *, double *),
                double dt, int neq);

void RK2Step   (double t, double *Y, void (*dYdt)(double, double *, double *),
                double dt, int neq);

void RK4Step   (double t, double *Y, void (*dYdt)(double, double *, double *),
                double dt, int neq);

void RK5Step   (double t, double *Y, void (*dYdt)(double, double *, double *),
                double dt, int neq);

void PositionVerletStep (double *x, double *v, void (*acc)(double *,double *),
                         double dt, int neq);

