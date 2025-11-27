
#include "ode_func.h"
#include "param.h"



void betaFPUT (double t, std::vector<double> &Y, std::vector<double> &R)
//
//RSH of motion of system
///////////////////////////////////////////
{
  double Tl = p.dparams["Tl"];
  double Tr = p.dparams["Tr"];
  double m = p.dparams["m"];
  double a = p.dparams["a"];
  double thetaL = p.dparams["thetaL"];
  double thetaR = p.dparams["thetaR"];
  double chi = p.dparams["chi"];
  double bet = p.dparams["beta"];
  int dim = p.iparams["dim"];
  int N = p.iparams["N"];

  int i,j,k,n;
  double r1 =0.0, r2 =0.0;
  k = 2*dim;

  double Psi1, Psi2; // DEFINISCO LE VARIABILI DEL TERMOSTATO
  double kin1 =0.0, kinN=0.0; // energia kin della prima e ultima
  // particlella mobile
  Psi1 = Y[k* (N+2)];
  Psi2 = Y[k* (N+2)+1];

  for (i=0; i<dim; i++) {
    kin1 += Y[k+dim+i]*Y[k+dim+i]/m;
    kinN += Y[k*N+dim+i]*Y[k*N+dim+i]/m;
  }

  //eq. for termostat
  R[k* (N+2)] = (kin1/Tl -1.0)/ (thetaL*thetaL); //left
  R[k* (N+2)+1] = (kinN/Tr -1.0)/ (thetaR*thetaR); //right


  ////////////////////////////////////
  i = 0; //prima particella fissa

  R[0] = 0.0; //la prima particella rimane ferma
  R[dim] = 0.0;

  ////////////////////////////////////

  ////////////////////////////////////

  i = 1; //la prima particella che si può muovere
  n = k*i;

  r1 = (Y[n+k] - Y[n] - a);
  r2 = (Y[n] - Y[n-k] - a);

  // eq. for position
  R[n] = Y[n+dim]/m;

  // eq. for momentum
  R[n+dim] = chi* (Y[n+k] + Y[n-k] - 2.0*Y[n] ) +
             bet* (r1*r1*r1 - r2*r2*r2) - Psi1*R[n];

  ////////////////////////////////////


  // le altre N-2 masse

  for (i = 2; i < N; i++) {
    n = k*i;

    r2 =r1;
    r1 = (Y[n+k] - Y[n] - a);

    // eq. for position
    R[n] = Y[n+dim]/m;

    // eq. for momentum
    R[n+dim] = chi* (Y[n+k] + Y[n-k] - 2.0*Y[n] ) +
               bet* (r1*r1*r1 - r2*r2*r2);
  }

  //l' ultima massa mobile
  i = N;
  n = k*i;
  r2 = r1;
  r1 = (Y[n+k] - Y[n] - a);

  // eq. for position
  R[n] = Y[n+dim]/m;

  // eq. for momentum
  R[n+dim] = chi* (Y[n+k] + Y[n-k] - 2.0*Y[n] ) +
             bet* (r1*r1*r1 - r2*r2*r2) - Psi2*R[n] ;


  ////////////////////////////////////
  i = N+1; //la massa finale fissa
  n = k*i;

  R[n] = 0.0; //è sempre ferma
  R[n+dim] = 0.0;

  ////////////////////////////////////
}

void AlfaBetaFPUT_initial (double t, std::vector<double> &Y, std::vector<double> &R)
//
//RSH of motion of system
///////////////////////////////////////////
{
  const double Tl = p.dparams["Tl"];
  const double Tr = Tl; // inizializzo a temperatura uguale al termostato di sinistra
  const double m = p.dparams["m"];
  const double a = p.dparams["a"];
  const double thetaL = p.dparams["thetaL"];
  const double thetaR = p.dparams["thetaR"];
  const double chi = p.dparams["chi"];
  const double bet = p.dparams["beta"];
  const double Alpha = p.dparams["alpha"];
  const int dim = p.iparams["dim"];
  const int N = p.iparams["N"];

  int i,j,k,n;
  double r1 =0.0, r2 =0.0;
  k = 2*dim;

  double Psi1, Psi2; // DEFINISCO LE VARIABILI DEL TERMOSTATO
  double kin1 =0.0, kinN=0.0; // energia kin della prima e ultima
  // particlella mobile
  Psi1 = Y[k* (N+2)];
  Psi2 = Y[k* (N+2)+1];

  for (i=0; i<dim; i++) {
    kin1 += Y[k+dim+i]*Y[k+dim+i]/m;
    kinN += Y[k*N+dim+i]*Y[k*N+dim+i]/m;
  }

  //eq. for termostat
  R[k* (N+2)] = (kin1/Tl -1.0)/ (thetaL*thetaL); //left
  R[k* (N+2)+1] = (kinN/Tr -1.0)/ (thetaR*thetaR); //right


  ////////////////////////////////////
  i = 0; //prima particella fissa

  R[0] = 0.0; //la prima particella rimane ferma
  R[dim] = 0.0;

  ////////////////////////////////////

  i = 1; //la prima particella che si può muovere
  n = k*i;

  r1 = (Y[n+k] - Y[n] - a);
  r2 = (Y[n] - Y[n-k] - a);


  // eq. for position
  R[n] = Y[n+dim]/m;

  // eq. for momentum
  R[n+dim] = chi* (Y[n+k] + Y[n-k] - 2.0*Y[n] ) +
             Alpha* (r1*r1 - r2*r2) +
             bet* (r1*r1*r1 - r2*r2*r2) - Psi1*R[n];

  ////////////////////////////////////

  // le altre N-2 masse

  for (i = 2; i < N; i++) {
    n = k*i;

    r2 =r1;
    r1 = (Y[n+k] - Y[n] - a);

    // eq. for position
    R[n] = Y[n+dim]/m;

    // eq. for momentum
    R[n+dim] = chi* (Y[n+k] + Y[n-k] - 2.0*Y[n] ) +
               Alpha* (r1*r1 - r2*r2)+ bet* (r1*r1*r1 - r2*r2*r2);
  }

  //l' ultima massa mobile
  i = N;
  n = k*i;
  r2 = r1;
  r1 = (Y[n+k] - Y[n] - a);

  // eq. for position
  R[n] = Y[n+dim]/m;

  // eq. for momentum
  R[n+dim] = chi* (Y[n+k] + Y[n-k] - 2.0*Y[n] ) +
             Alpha* (r1*r1 - r2*r2)+
             bet* (r1*r1*r1 - r2*r2*r2) - Psi2*R[n] ;


  ////////////////////////////////////
  i = N+1; //la massa finale fissa
  n = k*i;

  R[n] = 0.0; //è sempre ferma
  R[n+dim] = 0.0;

  ////////////////////////////////////
}


void AlfaBetaFPUT (double t, std::vector<double> &Y, std::vector<double> &R)
//
//RSH of motion of system
///////////////////////////////////////////
{
  const double Tl = p.dparams["Tl"];
  const double Tr = p.dparams["Tr"];
  const double m = p.dparams["m"];
  const double a = p.dparams["a"];
  const double thetaL = p.dparams["thetaL"];
  const double thetaR = p.dparams["thetaR"];
  const double chi = p.dparams["chi"];
  const double bet = p.dparams["beta"];
  const double Alpha = p.dparams["alpha"];
  const int dim = p.iparams["dim"];
  const int N = p.iparams["N"];

  int i,j,k,n;
  double r1 =0.0, r2 =0.0;
  k = 2*dim;

  double Psi1, Psi2; // DEFINISCO LE VARIABILI DEL TERMOSTATO
  double kin1 =0.0, kinN=0.0; // energia kin della prima e ultima
  // particlella mobile
  Psi1 = Y[k* (N+2)];
  Psi2 = Y[k* (N+2)+1];

  for (i=0; i<dim; i++) {
    kin1 += Y[k+dim+i]*Y[k+dim+i]/m;
    kinN += Y[k*N+dim+i]*Y[k*N+dim+i]/m;
  }

  //eq. for termostat
  R[k* (N+2)] = (kin1/Tl -1.0)/ (thetaL*thetaL); //left
  R[k* (N+2)+1] = (kinN/Tr -1.0)/ (thetaR*thetaR); //right


  ////////////////////////////////////
  i = 0; //prima particella fissa

  R[0] = 0.0; //la prima particella rimane ferma
  R[dim] = 0.0;

  ////////////////////////////////////

  i = 1; //la prima particella che si può muovere
  n = k*i;

  r1 = (Y[n+k] - Y[n] - a);
  r2 = (Y[n] - Y[n-k] - a);


  // eq. for position
  R[n] = Y[n+dim]/m;

  // eq. for momentum
  R[n+dim] = chi* (Y[n+k] + Y[n-k] - 2.0*Y[n] ) +
             Alpha* (r1*r1 - r2*r2) +
             bet* (r1*r1*r1 - r2*r2*r2) - Psi1*R[n];

  ////////////////////////////////////

  // le altre N-2 masse

  for (i = 2; i < N; i++) {
    n = k*i;

    r2 =r1;
    r1 = (Y[n+k] - Y[n] - a);

    // eq. for position
    R[n] = Y[n+dim]/m;

    // eq. for momentum
    R[n+dim] = chi* (Y[n+k] + Y[n-k] - 2.0*Y[n] ) +
               Alpha* (r1*r1 - r2*r2)+ bet* (r1*r1*r1 - r2*r2*r2);
  }

  //l' ultima massa mobile
  i = N;
  n = k*i;
  r2 = r1;
  r1 = (Y[n+k] - Y[n] - a);

  // eq. for position
  R[n] = Y[n+dim]/m;

  // eq. for momentum
  R[n+dim] = chi* (Y[n+k] + Y[n-k] - 2.0*Y[n] ) +
             Alpha* (r1*r1 - r2*r2)+
             bet* (r1*r1*r1 - r2*r2*r2) - Psi2*R[n] ;


  ////////////////////////////////////
  i = N+1; //la massa finale fissa
  n = k*i;

  R[n] = 0.0; //è sempre ferma
  R[n+dim] = 0.0;

  ////////////////////////////////////
}



void LJ (double t, std::vector<double> &Y, std::vector<double> &R)
{

  double Tl = p.dparams["Tl"];
  double Tr = p.dparams["Tr"];
  double m = p.dparams["m"];
  double thetaL = p.dparams["thetaL"];
  double thetaR = p.dparams["thetaR"];
  int dim = p.iparams["dim"];
  int N = p.iparams["N"];

  int i,j,k,n;
  double r1 =0.0, r2 =0.0;
  k = 2*dim;

  double Psi1, Psi2; // DEFINISCO LE VARIABILI DEL TERMOSTATO
  double kin1 =0.0, kinN=0.0; // energia kin della prima e ultima
  // particlella mobile
  Psi1 = Y[k* (N+2)];
  Psi2 = Y[k* (N+2)+1];

  for (i=0; i<dim; i++) {
    kin1 += Y[k+dim+i]*Y[k+dim+i]/m;
    kinN += Y[k*N+dim+i]*Y[k*N+dim+i]/m;
  }

  //eq. for termostat
  R[k* (N+2)] = (kin1/Tl -1.0)/ (thetaL*thetaL); //left
  R[k* (N+2)+1] = (kinN/Tr -1.0)/ (thetaR*thetaR); //right


  ////////////////////////////////////
  i = 0; //prima particella fissa

  R[0] = 0.0; //la prima particella rimane ferma
  R[dim] = 0.0;

  ////////////////////////////////////

  ////////////////////////////////////

  i = 1; //la prima particella che si può muovere
  n = k*i;

  r1 = Y[n+k] - Y[n];
  r2 = Y[n] - Y[n-k];


  // eq. for position
  R[n] = Y[n+dim]/m;

  // eq. for momentum
  R[n+dim] = -12.0* (1.0/pow (r1,13) - 1.0/pow (r1,7)
                     - 1.0/pow (r2,13) + 1.0/pow (r2,7)) - Psi1*R[n];



  ////////////////////////////////////


  // le altre N-2 masse

  for (i = 2; i < N; i++) {
    n = k*i;

    r2 =r1;
    r1 = Y[n+k] - Y[n];

    // eq. for position
    R[n] = Y[n+dim]/m;

    // eq. for momentum
    R[n+dim] = -12.0* (1.0/pow (r1,13) - 1.0/pow (r1,7)
                       - 1.0/pow (r2,13) + 1.0/pow (r2,7));
  }

  //l' ultima massa mobile
  i = N;
  n = k*i;
  r2 = r1;
  r1 = Y[n+k] - Y[n];

  // eq. for position
  R[n] = Y[n+dim]/m;

  // eq. for momentum
  R[n+dim] = -12.0* (1.0/pow (r1,13) - 1.0/pow (r1,7)
                     - 1.0/pow (r2,13) + 1.0/pow (r2,7)) - Psi2*R[n] ;


  ////////////////////////////////////
  i = N+1; //la massa finale fissa
  n = k*i;

  R[n] = 0.0; //è sempre ferma
  R[n+dim] = 0.0;

  ////////////////////////////////////
}


void LepriChain (double t, std::vector<double> &Y, std::vector<double> &R)
//
// LepriChain
// ----------
// Calcola il lato destro (RHS) del sistema di equazioni del moto per una
// catena mono-dimensionale di masse accoppiate secondo un modello di tipo
// FPU-β (potenziale quartico), con condizioni al contorno fisse e
// termostati di Nosé–Hoover applicati alla prima e all’ultima massa mobile.
//
// Il potenziale totale implementato è:
//
//     V = Σ_i [ 1/2 (x_{i+1} - x_i - a)^2
//              + 1/2 (x_i - i a)^2
//              + 1/4 (x_i - i a)^4 ]
//
// dove x_i sono le coordinate delle masse e a è la distanza di equilibrio.
// Nel codice la forza derivante dal potenziale locale è:
//        F_local = - (r + r^3),  con r = (x_i - i a)
//
// Le equazioni integrate sono:
//   - ẋ_i = p_i / m
//   - ṗ_i = (x_{i+1} + x_{i-1} - 2 x_i) + F_local
//
// Per le prime e ultime masse mobili sono presenti i termini dei termostati
// di Nosé–Hoover:
//
//   ṗ_1  ← ṗ_1  - Ψ_L ẋ_1
//   ṗ_N  ← ṗ_N  - Ψ_R ẋ_N
//
// Le variabili Ψ_L e Ψ_R evolvono secondo:
//   Ψ̇ = (K/T - 1) / θ²
//
// La funzione richiede il vettore di stato Y contenente, per ogni massa,
// prima le posizioni e poi i momenti, ed inserisce in R le loro derivate.
//
///////////////////////////////////////////
{
  double Tl = p.dparams["Tl"];
  double Tr = p.dparams["Tr"];
  double m = p.dparams["m"];
  double a = p.dparams["a"];
  double thetaL = p.dparams["thetaL"];
  double thetaR = p.dparams["thetaR"];
  double chi = p.dparams["chi"];
  double bet = p.dparams["beta"];
  int dim = p.iparams["dim"];
  int N = p.iparams["N"];
  double spring = 1.0;

  int i,j,k,n;
  double r1 =0.0, r2 =0.0;
  k = 2*dim;

  double Psi1, Psi2; // DEFINISCO LE VARIABILI DEL TERMOSTATO
  double kin1 =0.0, kinN=0.0; // energia kin della prima e ultima
  // particlella mobile
  Psi1 = Y[k* (N+2)];
  Psi2 = Y[k* (N+2)+1];

  for (i=0; i<dim; i++) {
    kin1 += Y[k+dim+i]*Y[k+dim+i]/m;
    kinN += Y[k*N+dim+i]*Y[k*N+dim+i]/m;
  }

  //eq. for termostat
  R[k* (N+2)] = (kin1/Tl -1.0)/ (thetaL*thetaL); //left
  R[k* (N+2)+1] = (kinN/Tr -1.0)/ (thetaR*thetaR); //right


  ////////////////////////////////////
  i = 0; //prima particella fissa

  R[0] = 0.0; //la prima particella rimane ferma
  R[dim] = 0.0;

  ////////////////////////////////////

  ////////////////////////////////////

  i = 1; //la prima particella che si può muovere
  n = k*i;

  r1 = (Y[n] - a);


  // eq. for position
  R[n] = Y[n+dim]/m;

  // eq. for momentum
  R[n+dim] = Y[n+k] + Y[n-k] - 2.0*Y[n] +
             spring* (- r1 - r1*r1*r1 ) - Psi1*R[n];

  ////////////////////////////////////


  // le altre N-2 masse

  for (i = 2; i < N; i++) {
    n = k*i;

    r1 = (Y[n] - i*a);

    // eq. for position
    R[n] = Y[n+dim]/m;

    // eq. for momentum
    R[n+dim] = Y[n+k] + Y[n-k] - 2.0*Y[n] +
               spring* (- r1 - r1*r1*r1);
  }

  //l' ultima massa mobile
  i = N;
  n = k*i;
  r1 = (Y[n] - i*a);

  // eq. for position
  R[n] = Y[n+dim]/m;

  // eq. for momentum
  R[n+dim] = Y[n+k] + Y[n-k] - 2.0*Y[n] +
             spring* (- r1 - r1*r1*r1) - Psi2*R[n] ;


  ////////////////////////////////////
  i = N+1; //la massa finale fissa
  n = k*i;

  R[n] = 0.0; //è sempre ferma
  R[n+dim] = 0.0;

  ////////////////////////////////////
}

void LepriChain_initial (double t, std::vector<double> &Y, std::vector<double> &R)
//
// LepriChain
// ----------
// Calcola il lato destro (RHS) del sistema di equazioni del moto per una
// catena mono-dimensionale di masse accoppiate secondo un modello di tipo
// FPU-β (potenziale quartico), con condizioni al contorno fisse e
// termostati di Nosé–Hoover applicati alla prima e all’ultima massa mobile.
//
// Il potenziale totale implementato è:
//
//     V = Σ_i [ 1/2 (x_{i+1} - x_i - a)^2
//              + 1/2 (x_i - i a)^2
//              + 1/4 (x_i - i a)^4 ]
//
// dove x_i sono le coordinate delle masse e a è la distanza di equilibrio.
// Nel codice la forza derivante dal potenziale locale è:
//        F_local = - (r + r^3),  con r = (x_i - i a)
//
// Le equazioni integrate sono:
//   - ẋ_i = p_i / m
//   - ṗ_i = (x_{i+1} + x_{i-1} - 2 x_i) + F_local
//
// Per le prime e ultime masse mobili sono presenti i termini dei termostati
// di Nosé–Hoover:
//
//   ṗ_1  ← ṗ_1  - Ψ_L ẋ_1
//   ṗ_N  ← ṗ_N  - Ψ_R ẋ_N
//
// Le variabili Ψ_L e Ψ_R evolvono secondo:
//   Ψ̇ = (K/T - 1) / θ²
//
// La funzione richiede il vettore di stato Y contenente, per ogni massa,
// prima le posizioni e poi i momenti, ed inserisce in R le loro derivate.
//
///////////////////////////////////////////
{
  const double Tl = p.dparams["Tl"];
  const double Tr = Tl; // inizializzo a temperatura uguale al termostato di sinistra
  double m = p.dparams["m"];
  double a = p.dparams["a"];
  double thetaL = p.dparams["thetaL"];
  double thetaR = p.dparams["thetaR"];
  double chi = p.dparams["chi"];
  double bet = p.dparams["beta"];
  int dim = p.iparams["dim"];
  int N = p.iparams["N"];
  double spring = 1.0;

  int i,j,k,n;
  double r1 =0.0, r2 =0.0;
  k = 2*dim;

  double Psi1, Psi2; // DEFINISCO LE VARIABILI DEL TERMOSTATO
  double kin1 =0.0, kinN=0.0; // energia kin della prima e ultima
  // particlella mobile
  Psi1 = Y[k* (N+2)];
  Psi2 = Y[k* (N+2)+1];

  for (i=0; i<dim; i++) {
    kin1 += Y[k+dim+i]*Y[k+dim+i]/m;
    kinN += Y[k*N+dim+i]*Y[k*N+dim+i]/m;
  }

  //eq. for termostat
  R[k* (N+2)] = (kin1/Tl -1.0)/ (thetaL*thetaL); //left
  R[k* (N+2)+1] = (kinN/Tr -1.0)/ (thetaR*thetaR); //right


  ////////////////////////////////////
  i = 0; //prima particella fissa

  R[0] = 0.0; //la prima particella rimane ferma
  R[dim] = 0.0;

  ////////////////////////////////////

  ////////////////////////////////////

  i = 1; //la prima particella che si può muovere
  n = k*i;

  r1 = (Y[n] - a);


  // eq. for position
  R[n] = Y[n+dim]/m;

  // eq. for momentum
  R[n+dim] = Y[n+k] + Y[n-k] - 2.0*Y[n] +
             spring* (- r1 - r1*r1*r1 ) - Psi1*R[n];

  ////////////////////////////////////


  // le altre N-2 masse

  for (i = 2; i < N; i++) {
    n = k*i;

    r1 = (Y[n] - i*a);

    // eq. for position
    R[n] = Y[n+dim]/m;

    // eq. for momentum
    R[n+dim] = Y[n+k] + Y[n-k] - 2.0*Y[n] +
               spring* (- r1 - r1*r1*r1);
  }

  //l' ultima massa mobile
  i = N;
  n = k*i;
  r1 = (Y[n] - i*a);

  // eq. for position
  R[n] = Y[n+dim]/m;

  // eq. for momentum
  R[n+dim] = Y[n+k] + Y[n-k] - 2.0*Y[n] +
             spring* (- r1 - r1*r1*r1) - Psi2*R[n] ;


  ////////////////////////////////////
  i = N+1; //la massa finale fissa
  n = k*i;

  R[n] = 0.0; //è sempre ferma
  R[n+dim] = 0.0;

  ////////////////////////////////////
}
