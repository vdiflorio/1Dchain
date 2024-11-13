#include <iostream>
#include <cmath>
#include <iomanip>
#include <fstream>
#include <vector>

void EulerStep (double t, std::vector<double>  &Y, void (*dYdt)(double, std::vector<double> &, std::vector<double> &),
                double dt, int neq);

void RK2Step   (double t, std::vector<double>  &Y, void (*dYdt)(double, std::vector<double> &, std::vector<double> &),
                double dt, int neq);

void RK4Step   (double t, std::vector<double>  &Y, void (*dYdt)(double, std::vector<double> &, std::vector<double> &),
                double dt, int neq);

void RK5Step   (double t, std::vector<double>  &Y, void (*dYdt)(double, std::vector<double> &, std::vector<double> &),
                double dt, int neq);

void PositionVerletStep (double *x, double *v, void (*acc)(double *,double *),
                         double dt, int neq);

