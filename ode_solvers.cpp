#include "ode_solvers.h"
#include "param.h"


void EulerStep (double t, std::vector<double> &Y,
                std::function<void (double, std::vector<double> &, std::vector<double> &)> dYdt,
                double dt, int neq)
// Take one step dt using Euler method for the solution of dY/dt = rhs.
// Here neq is the number of ODE (the dimensionality of Y[]) and *dYdt() is
// a pointer to the function that calculates the right hand side
// of the system of equations.
{
  int n;
  static std::vector<double> rhs (neq);
  rhs.assign (neq, 0.0);

  dYdt (t, Y, rhs);

  for (n = 0; n < neq; n++) {
    Y[n] += dt*rhs[n];
  }

}



void RK2Step (double t, std::vector<double> &Y,
              std::function<void (double, std::vector<double> &, std::vector<double> &)> dYdt,
              double dt, int neq)
// Take one step dt using Euler method for the solution of dY/dt = rhs.
// Here neq is the number of ODE (the dimensionality of Y[]) and *dYdt() is
// a pointer to the function that calculates the right hand side
// of the system of equations.
{
  int n;
  static std::vector<double> y1, k1, k2;
  y1.assign (neq, 0.0);
  k1.assign (neq, 0.0);
  k2.assign (neq, 0.0);

  dYdt (t, Y, k1);

  for (n = 0; n < neq; n++) {
    y1[n] = Y[n]+0.5*dt*k1[n];
  }

  dYdt (t+0.5*dt, y1, k2);

  for (n = 0; n < neq; n++) {
    Y[n] += dt*k2[n];
  }
}



void RK4Step (double t, std::vector<double> &Y,
              std::function<void (double, std::vector<double> &, std::vector<double> &)> dYdt,
              double dt, int neq)
// Take one step dt using Euler method for the solution of dY/dt = rhs.
// Here neq is the number of ODE (the dimensionality of Y[]) and *dYdt() is
// a pointer to the function that calculates the right hand side
// of the system of equations.
{
  int n;
  //static std::vector<double> y1, y2, y3, k1, k2, k3, k4;
  static std::vector<double> y1 (neq), y2 (neq), y3 (neq), k1 (neq), k2 (neq), k3 (neq), k4 (neq);
// Clear or reset vectors only if needed
  std::fill (y1.begin(), y1.end(), 0.0);
  std::fill (y2.begin(), y2.end(), 0.0);
  std::fill (y3.begin(), y3.end(), 0.0);

  dYdt (t, Y, k1);


  for (n = 0; n < neq; n++) {
    y1[n] = Y[n]+ 0.5*dt*k1[n];
  }

  dYdt (t+0.5*dt, y1, k2);

  for (n = 0; n < neq; n++) {
    y2[n] = Y[n]+ 0.5*dt*k2[n];
  }

  dYdt (t+0.5*dt, y2, k3);

  for (n = 0; n < neq; n++) {
    y3[n] = Y[n]+ dt*k3[n];
  }

  dYdt (t+dt, y3, k4);

  for (n = 0; n < neq; n++) {
    Y[n] += (dt/6.0)* (k1[n] + 2.0*k2[n] + 2.0*k3[n] +k4[n]);
  }
}


void RK5Step (double t, std::vector<double> &Y,
              std::function<void (double, std::vector<double> &, std::vector<double> &)> dYdt,
              double dt, int neq)
// Take one step dt using Euler method for the solution of dY/dt = rhs.
// Here neq is the number of ODE (the dimensionality of Y[]) and *dYdt() is
// a pointer to the function that calculates the right hand side
// of the system of equations.
{
  int n;
  static std::vector<double> y1 (neq), y2 (neq), y3 (neq), y4 (neq), y5 (neq);
  static std::vector<double> k1 (neq), k2 (neq), k3 (neq), k4 (neq), k5 (neq), k6 (neq);


  dYdt (t, Y, k1);

  for (n = 0; n < neq; n++) {
    y1[n] = Y[n] + 0.25*dt*k1[n];
  }

  dYdt (t+0.25*dt, y1, k2);

  for (n = 0; n < neq; n++) {
    y2[n] = Y[n] + 0.125*dt*k1[n] + 0.125*dt*k2[n];
  }

  dYdt (t+0.25*dt, y2, k3);

  for (n = 0; n < neq; n++) {
    y3[n] = Y[n] - 0.5*dt*k2[n] + dt*k3[n];
  }

  dYdt (t+0.5*dt, y3, k4);

  for (n = 0; n < neq; n++) {
    y4[n] = Y[n] + (3./16.)*dt*k1[n] + (9./16.)*dt*k4[n];
  }

  dYdt (t+ (3./4.)*dt, y4, k5);

  for (n = 0; n < neq; n++) {
    y5[n] = Y[n] - (3./7.)*dt*k1[n] + (2./7.)*dt*k2[n] +
            (12./7.)*dt*k3[n] - (12./7.)*dt*k4[n] + (8./7.)*dt*k5[n];
  }

  dYdt (t+dt, y5, k6);

  for (n = 0; n < neq; n++) {
    Y[n] += (dt/90.)* (7.*k1[n] + 32.*k3[n] + 12.*k4[n] +32.*k5[n] + 7.*k6[n]);
  }
}

void PositionVerletStep (double *x, double *v, void (*acc) (double *,double *),
                         double dt, int neq)
{
  double a[neq];
  int i;

  for (i=0; i<neq; i++) {
    x[i] += 0.5*dt*v[i];
  }

  acc (x,a);

  for (i=0; i<neq; i++) {
    v[i] += dt*a[i];
  }

  for (i=0; i<neq; i++) {
    x[i] += 0.5*dt*v[i];
  }

}

