//
//
////////////////////////////////////////////////////////////////////////////

#include "ode_solvers.h"
#include "ode_func.h"


void linear_fit (int ndata, double *X, double *Y, double& beta1, double& q);


int main()
{
  ofstream ddata;
  ofstream tdata;
  ofstream pdata;   
  ddata.open("densita.dat");
  tdata.open("temperatura.dat"); 
  pdata.open("dati_iniziali.dat");

  ddata   << setiosflags(ios::scientific); 
  ddata   << setprecision(8);
  tdata  << setiosflags(ios::scientific); 
  tdata  << setprecision(8);
  pdata   << setiosflags(ios::scientific); 
  pdata   << setprecision(8);

  int      neq = (N+2)*2*dim + 2;
  double   X[neq];
  long int step = 8000000;
  long int h;
  long int no_step  = 5000;
  long int progress = 0;
  int      k,i,j,l,n;
  double   t, dt;
  double   X0[dim];  //unità di spaziatura
  double   X_eq[dim*(N+2)]; //vettore contenete posizioni di equilibrio meccanico (potrebbe essere cancellato)
  double   T[N];
  double   R[dim*(N+2)];
  double   R_m[dim*(N+2)];
  double   S[N];
  long int count = step - no_step;  //# di steps per le medie temporali
  double   Xm[(N+2)*dim];
  double   Var[(N+2)*dim], Dev[(N+2)*dim]; //varianza nello spostamento
  
  k  = 2*dim;   
  t  = 0.0;    //tempo zero
  dt = 1.e-3;  //intervallo di integrazione

  X0[0] = a;   // nodes only along x-axis
  #if dim > 1
   for(i = 1; i<dim; i++) {X0[i] = 0.0;}
  #endif
  
  
  // INIZILIZZAZIONE VETTORI
  for(i=0; i < neq; i++) { X[i] = 0.0; } 
  for(i=0;i<N; i++) { 
    T[i]=0.0; 
    S[i]  = 0.0;
  }
  for(i=0;i<(N+2)*dim; i++){
    R_m[i] = 0.0;
    R[i]   = 0.0;  
    Xm[i]  = 0.0;
    Var[i] = 0.0;
    Dev[i] = 0.0;
  }
  
  X[k*(N+2)]   = 1.0; //inizializzazione termostato di sinistra
  X[k*(N+2)+1] = 1.0; //inizializzazione termostato di destra
  ////////////////////////////////////////////////////////////////
  
  // Nodes position at mechanical equilibrium (potrebbe non essere utile per il codice)
  for(j=0;j<=N+1; j++){
    l = dim*j;
    for(i=0;i<dim; i++){
      X_eq[l+i] = j*X0[i];
    }
  }
  ////////////////////////////////////////////////////////////////
  
  double alfa = 1.0;  // fattore di stretching
  
  // CONDIZIONE INIZIALE per posizione e velocità
  srand48(time(NULL));          // Initialize the sequence
  for (j= 1; j<=N; j++){
    n=k*j;
    l=dim*j;
    for (i= 0; i<dim; ++i){
      X[n+i]= X_eq[l+i]*alfa;// + drand48()*0.8 - 0.4;  //posizione
      X[n+i+dim]= 0.1;  //velocità
      pdata <<X[n+i]<< "  " << X[n+i+dim]<<endl;
    }
  }
  X[(N+1)*k] = X_eq[(N+1)*dim]*alfa; //l'ultima particella è fissa
  
  
  pdata <<endl<<endl;

  //EVOLUZIONE SISTEMA
  for(h=1; h<=step; h++){
    if(h%(int(step)/100) == 0){
      cout << progress << "% executed" << endl; // Show the progress of simulation 
      progress++;
    }
    RK4Step(t, X, AlfaBeta, dt,neq);   // integration of the function
    t += dt; 
    
    if(h>no_step){  //no prendo dati per un numero di passi uguali a no_step  
      for(j=0; j<=N+1; j++){
        n = k*j;
        l = dim*j;
        for(i=0; i<dim; i++){
          Xm[l+i] += X[n+i]; 
          Var[l+i] += X[n+i]*X[n+i];
        }
      }
      for(j=0; j<N;j++){
        n = k*(j+1);
        for(i=0; i<dim; i++){          
          T[j] += X[n+dim+i]*X[n+dim+i]/m;
        } 
      }  
    }   
  }
  
  
  for(i=0;i<(N+2)*dim; i++){
    Xm[i]   = Xm[i]/(double)count;
    Var[i]  = Var[i]/(double)count - Xm[i]*Xm[i] ;
    Dev[i]  = sqrt(Var[i]);
  }
  
  for(j=0; j<=N+1; j++){
    n = dim*j;
    for(i=0; i<dim; i++){
      R_m[n+i] = Xm[n+i+dim] - Xm[n+i];
    }
  }
   
      
  for(j=0; j<N; j++){
    l = dim*(j+1);
    S[j] =0.0;
    for(i=0; i<dim; i++){
      S[j] += R_m[l+i]*R_m[l+i];  
    }
    S[j] = sqrt(S[j]);
    T[j] = T[j]/(double)count;
  }   
  
  for(j=0; j<N; j++){        
    ddata  << 1.0*j/N << " " << S[j] << endl;     
  }
  ddata << endl << endl;
  
  
  //linear fit of linear relation and density
  double c1, c2;
  
  linear_fit (N, S, T, c1, c2);
  cout <<"c1 = "<< c1 << "   c2 = " << c2<<endl;
  
  for(j=0; j<N; j++){        
    ddata  << 1.0*j/N << " " << c1*S[j] + c2 << endl;   
    tdata << 1.0*j/N << " " << T[j] << endl;  
  }

  ddata <<endl<<endl;
  ddata << c1 << " " << c2 <<endl;
  
  ////////////////////////////////////////////////////////////////  
  
  ddata.close();
  tdata.close();  
  pdata.close();  
  return 0;
}




//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////

void linear_fit (int ndata, double *X, double *Y, double& beta1, double& q){
 
  int  i;
  double sumX=0, sumX2=0, sumY=0, sumXY=0;
  
  /* Calculating Required Sum */
  for(i=0;i<ndata;i++){
    sumX  +=  X[i];
    sumX2 +=  X[i]*X[i];
    sumY  +=  Y[i];
    sumXY +=  X[i]*Y[i];
  }
  /* Calculating m and q */
  beta1 = (ndata*sumXY-sumX*sumY)/(ndata*sumX2-sumX*sumX);
  q = (sumY - beta1*sumX)/ndata;
}
