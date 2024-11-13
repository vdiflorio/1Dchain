//
//
////////////////////////////////////////////////////////////////////////////

#include "ode_solvers.h"
#include "ode_func.h"
#include "utils.h"


int main()
{

  int      neq = (N+2)*2*dim + 2;
  std::vector<double> X(neq);
  long int step = 8000000;
  long int h;
  long int no_step  = 3000000;
  long int progress = 0;
  int      k,i,j,l,n;
  double   t, dt;
  double   X0[dim];  //unità di spaziatura
  double   X_eq[dim*(N+2)]; //vettore contenete posizioni di equilibrio meccanico (potrebbe essere cancellato)

  std::vector<std::vector<double>> condizioni;

  int N_catene = 1;

  
  k  = 2*dim;   
  t  = 0.0;    //tempo zero
  dt = 1.e-3;  //intervallo di integrazione

  X0[0] = a;   // nodes only along x-axis
  // Nodes position at mechanical equilibrium (potrebbe non essere utile per il codice)
  for(j=0;j<=N+1; j++){
    l = dim*j;
    for(i=0;i<dim; i++){
      X_eq[l+i] = j*X0[i];
    }
  }
  
  for (int ii = 0; ii<N_catene; ++ii){
    std::cout << "\n\n CATENA NUMERO: " << ii+1 << std::endl << std::endl;

    // INIZILIZZAZIONE VETTORI
    for(i=0; i < neq; i++) { X[i] = 0.0; } 
    X[k*(N+2)]   = 1.0; //inizializzazione termostato di sinistra
    X[k*(N+2)+1] = 1.0; //inizializzazione termostato di destra
    ////////////////////////////////////////////////////////////////
        
    double alfa = 1.0;  // fattore di stretching
    
    // CONDIZIONE INIZIALE per posizione e velocità
    srand48(time(NULL));          // Initialize the sequence
    for (j= 1; j<=N; j++){
      n=k*j;
      l=dim*j;
      for (i= 0; i<dim; ++i){
        X[n+i]= X_eq[l+i]*alfa + drand48()*0.8 - 0.4;  //posizione
        X[n+i+dim]= drand48()*0.2 - 0.1;  //velocità
      }
    }
    X[(N+1)*k] = X_eq[(N+1)*dim]*alfa; //l'ultima particella è fissa
    
    

    //EVOLUZIONE SISTEMA
    for(h=1; h<=step - no_step; h++){
      if(h%(int(step)/100) == 0){
        std::cout << progress << "% executed" << std::endl; // Show the progress of simulation 
        progress++;
      }
      RK4Step(t, X, Chain1, dt,neq);   // integration of the function
      // RK4Step(t, X, AlfaBeta_corrected, dt,neq);   // integration of the function
      t += dt; 
    }
    for(h=step - no_step; h<=step; h++){  
      
      if(h%(int(step)/100) == 0){
        std::cout << progress << "% executed" << std::endl; // Show the progress of simulation 
        progress++;
      }  

      if (drand48()<0.02)
        condizioni.push_back(X);
   
    }
  }
  std::ofstream outFile("condizioni.bin", std::ios::binary);
  if (!outFile) {
    std::cerr << "Errore nell'apertura del file per scrittura!" << std::endl;
    return 1;
  }

  // Scrive tutti i dati in un'unica operazione per massima efficienza
  outFile.write(reinterpret_cast<char*>(condizioni.data()), condizioni.size() * neq * sizeof(double));

  outFile.close();
  if (!outFile.good()) {
    std::cerr << "Errore durante la scrittura dei dati nel file!" << std::endl;
    return 1;
  }

  std::cout << "Condizioni scritte con successo in formato binario." << std::endl;


  return 0;
}


