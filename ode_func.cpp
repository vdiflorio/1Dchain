
#include "ode_func.h"
#include "param.h"



void  betaFPUT (double t, std::vector<double> &Y, std::vector<double> &R)
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

  int     i,j,k,n;
  double  r1 =0.0, r2 =0.0;
  k = 2*dim;
  
  double Psi1, Psi2;          // DEFINISCO LE VARIABILI DEL TERMOSTATO
  double kin1 =0.0, kinN=0.0; // energia kin della prima e ultima 
                              // particlella mobile
  Psi1 = Y[k*(N+2)];
  Psi2 = Y[k*(N+2)+1];
  
  for(i=0; i<dim; i++){
    kin1 += Y[k+dim+i]*Y[k+dim+i]/m;
    kinN += Y[k*N+dim+i]*Y[k*N+dim+i]/m;
  }

  //eq. for termostat
  R[k*(N+2)]   = (kin1/Tl -1.0)/(thetaL*thetaL); //left
  R[k*(N+2)+1] = (kinN/Tr -1.0)/(thetaR*thetaR); //right
  
  
  //////////////////////////////////// 
  i = 0; //prima particella fissa
  
  R[0]   = 0.0;     //la prima particella rimane ferma
  R[dim] = 0.0;
  
  ////////////////////////////////////
  
  ////////////////////////////////////
  
  i = 1;  //la prima particella che si può muovere   
  n = k*i;
  
  r1 = (Y[n+k] - Y[n] - a);
  r2 = (Y[n] - Y[n-k] - a);   
  
  // eq. for position
  R[n] = Y[n+dim]/m;
   
  // eq. for momentum
  R[n+dim] = chi*(Y[n+k] + Y[n-k] - 2.0*Y[n] ) + 
             bet*(r1*r1*r1 - r2*r2*r2) - Psi1*R[n];
   
  ////////////////////////////////////
  
  
  // le altre N-2 masse
  
  for(i = 2; i < N; i++){
    n = k*i;
       
    r2 =r1;    
    r1 = (Y[n+k] - Y[n] - a);
    
    // eq. for position
    R[n] = Y[n+dim]/m;
   
    // eq. for momentum
    R[n+dim] = chi*(Y[n+k] + Y[n-k] - 2.0*Y[n] ) + 
               bet*(r1*r1*r1 - r2*r2*r2);
  }

  //l' ultima massa mobile 
  i = N;
  n = k*i;
  r2 = r1;
  r1 = (Y[n+k] - Y[n] - a);
  
  // eq. for position
  R[n] = Y[n+dim]/m;
   
  // eq. for momentum
  R[n+dim] = chi*(Y[n+k] + Y[n-k] - 2.0*Y[n] ) + 
             bet*(r1*r1*r1 - r2*r2*r2) - Psi2*R[n] ;   
  
    
  ////////////////////////////////////
  i = N+1; //la massa finale fissa
  n = k*i;
  
  R[n] = 0.0;        //è sempre ferma
  R[n+dim] = 0.0; 
  
  ////////////////////////////////////
}


//////////////////////////////////////////////////////////////////////


void   AlfaBetaFPUT (double t, std::vector<double> &Y, std::vector<double> &R)
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
  double Alpha = p.dparams["alpha"];
  int dim = p.iparams["dim"];
  int N = p.iparams["N"];

  int     i,j,k,n;
  double  r1 =0.0, r2 =0.0;
  k = 2*dim;
  
  double Psi1, Psi2;          // DEFINISCO LE VARIABILI DEL TERMOSTATO
  double kin1 =0.0, kinN=0.0; // energia kin della prima e ultima 
                              // particlella mobile
  Psi1 = Y[k*(N+2)];
  Psi2 = Y[k*(N+2)+1];
  
  for(i=0; i<dim; i++){
    kin1 += Y[k+dim+i]*Y[k+dim+i]/m;
    kinN += Y[k*N+dim+i]*Y[k*N+dim+i]/m;
  }

  //eq. for termostat
  R[k*(N+2)]   = (kin1/Tl -1.0)/(thetaL*thetaL); //left
  R[k*(N+2)+1] = (kinN/Tr -1.0)/(thetaR*thetaR); //right
  
  
  //////////////////////////////////// 
  i = 0; //prima particella fissa
  
  R[0]   = 0.0;     //la prima particella rimane ferma
  R[dim] = 0.0;
  
  ////////////////////////////////////
    
  i = 1;  //la prima particella che si può muovere   
  n = k*i;
  
  r1 = (Y[n+k] - Y[n] - a);
  r2 = (Y[n] - Y[n-k] - a);   
  

  // eq. for position
  R[n] = Y[n+dim]/m;
   
  // eq. for momentum
  R[n+dim] = chi*(Y[n+k] + Y[n-k] - 2.0*Y[n] ) + 
             Alpha*(r1*r1 - r2*r2) + 
             bet*(r1*r1*r1 - r2*r2*r2) - Psi1*R[n];
    
  ////////////////////////////////////
  
  // le altre N-2 masse
  
  for(i = 2; i < N; i++){
    n = k*i;
       
    r2 =r1;    
    r1 = (Y[n+k] - Y[n] - a);
    
    // eq. for position
    R[n] = Y[n+dim]/m;
   
    // eq. for momentum
    R[n+dim] = chi*(Y[n+k] + Y[n-k] - 2.0*Y[n] ) + 
               Alpha*(r1*r1 - r2*r2)+  bet*(r1*r1*r1 - r2*r2*r2);
  }

  //l' ultima massa mobile 
  i = N;
  n = k*i;
  r2 = r1;
  r1 = (Y[n+k] - Y[n] - a);
  
  // eq. for position
  R[n] = Y[n+dim]/m;
   
  // eq. for momentum
  R[n+dim] = chi*(Y[n+k] + Y[n-k] - 2.0*Y[n] ) + 
             Alpha*(r1*r1 - r2*r2)+ 
             bet*(r1*r1*r1 - r2*r2*r2) - Psi2*R[n] ;   
  
    
  ////////////////////////////////////
  i = N+1; //la massa finale fissa
  n = k*i;
  
  R[n] = 0.0;        //è sempre ferma
  R[n+dim] = 0.0; 
  
  ////////////////////////////////////
}



void LJ (double t, std::vector<double> &Y, std::vector<double> &R){
	
  double Tl = p.dparams["Tl"];
  double Tr = p.dparams["Tr"];
  double m = p.dparams["m"];
  double thetaL = p.dparams["thetaL"];
  double thetaR = p.dparams["thetaR"];
  int dim = p.iparams["dim"];
  int N = p.iparams["N"];

  int     i,j,k,n;
  double  r1 =0.0, r2 =0.0;
  k = 2*dim;
  
  double Psi1, Psi2;          // DEFINISCO LE VARIABILI DEL TERMOSTATO
  double kin1 =0.0, kinN=0.0; // energia kin della prima e ultima 
                              // particlella mobile
  Psi1 = Y[k*(N+2)];
  Psi2 = Y[k*(N+2)+1];
  
  for(i=0; i<dim; i++){
    kin1 += Y[k+dim+i]*Y[k+dim+i]/m;
    kinN += Y[k*N+dim+i]*Y[k*N+dim+i]/m;
  }

  //eq. for termostat
  R[k*(N+2)]   = (kin1/Tl -1.0)/(thetaL*thetaL); //left
  R[k*(N+2)+1] = (kinN/Tr -1.0)/(thetaR*thetaR); //right
  
  
  //////////////////////////////////// 
  i = 0; //prima particella fissa
  
  R[0]   = 0.0;     //la prima particella rimane ferma
  R[dim] = 0.0;
  
  ////////////////////////////////////
  
  ////////////////////////////////////
  
  i = 1;  //la prima particella che si può muovere   
  n = k*i;
  
  r1 = Y[n+k] - Y[n];
  r2 = Y[n] - Y[n-k];   
  

  // eq. for position
  R[n] = Y[n+dim]/m;
   
  // eq. for momentum
  R[n+dim] = -12.0*(1.0/pow(r1,13) - 1.0/pow(r1,7) 
  								- 1.0/pow(r2,13) + 1.0/pow(r2,7)) - Psi1*R[n];
   
    
  
  ////////////////////////////////////
  
  
  // le altre N-2 masse
  
  for(i = 2; i < N; i++){
    n = k*i;
       
    r2 =r1;    
    r1 = Y[n+k] - Y[n];
    
    // eq. for position
    R[n] = Y[n+dim]/m;
   
    // eq. for momentum
    R[n+dim] = -12.0*(1.0/pow(r1,13) - 1.0/pow(r1,7) 
    								- 1.0/pow(r2,13) + 1.0/pow(r2,7));
  }

  //l' ultima massa mobile 
  i = N;
  n = k*i;
  r2 = r1;
  r1 = Y[n+k] - Y[n];
  
  // eq. for position
  R[n] = Y[n+dim]/m;
   
  // eq. for momentum
  R[n+dim] = -12.0*(1.0/pow(r1,13) - 1.0/pow(r1,7) 
  								- 1.0/pow(r2,13) + 1.0/pow(r2,7)) - Psi2*R[n] ;   
  
    
  ////////////////////////////////////
  i = N+1; //la massa finale fissa
  n = k*i;
  
  R[n] = 0.0;        //è sempre ferma
  R[n+dim] = 0.0; 
  
  ////////////////////////////////////
}

	

