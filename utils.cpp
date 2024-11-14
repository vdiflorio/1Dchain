// File for auxiliary functions (ex. Omega)
#include "ode_func.h"
#include "ode_solvers.h"

#include "utils.h"

#include <random>


double omega_0(std::vector<double> &Y, double T){
    // Function to compute the dissipation function wrt the initial density
    int     j,k;   
    k = 2*dim;
    double  omega=0.0;
    
    for(j =0; j<dim; j++){
        omega+=Y[k*(N+2)+1]*Y[k*N+dim+j]*Y[k*N+dim+j]*(1.0/Tr-1.0/T)/m +
               Y[k*(N+2)]*Y[k+dim+j]*Y[k+dim+j]*(1.0/Tl-1.0/T)/m;
    } 

    return omega;
}


double TTCF(double (*obs)(std::vector<double> &), double omega,std::vector<double> &Y,double T){
    return omega*obs(Y);
double TTCF(double (*obs)(double *), double omega ,double *Y,double T){
    // Function to compute the TTCF



    return obs(Y)*omega;
}

double observable( double *Y){
    // Function to compute any observable, that is function of the phase space state


    return 1;
};


void read_conditions(std::vector<double>& condizioni, int num_condizioni, int neq) {
  
  std::ifstream inFile("condizioni.bin", std::ios::binary);
  if (!inFile) {
    std::cerr << "Errore nell'apertura del file per lettura!" << std::endl;
    return;
  }
  // Dimensione di una condizione in byte
  int dimensione_condizione = neq * sizeof(double); 
  
  // Ottieni la dimensione totale del file
  // Posizionarsi alla fine del file per ottenere la sua dimensione
  inFile.seekg(0, std::ios::end);
  std::streampos file_size = inFile.tellg();
  inFile.seekg(0, std::ios::beg);  // Torna all'inizio del file


  // Calcola il numero di condizioni nel file
  int numero_condizioni_tot = file_size / dimensione_condizione;

  if (num_condizioni > numero_condizioni_tot) {
    std::cerr << "Errore: il numero di condizioni richiesto eccede il numero di condizioni nel file." << std::endl;
    inFile.close();
    return;
  }


  // Generatore di numeri casuali
  std::mt19937 gen(std::random_device{}());  // Generatore basato su random_device
  std::uniform_int_distribution<int> dist(0, numero_condizioni_tot - 1); // Seleziona indici casuali

  // Riservare spazio per le condizioni
  condizioni.reserve(num_condizioni * neq);

  for (int i = 0; i < num_condizioni; ++i) {
    int indice_casuale = dist(gen); // Seleziona un indice casuale

    // Calcolare l'offset dell'indice
    std::streampos offset = indice_casuale * dimensione_condizione;

    // Posizionarsi nel file all'offset desiderato
    inFile.seekg(offset);

    // Leggere la condizione direttamente nel vettore 1D
    inFile.read(reinterpret_cast<char*>(&condizioni[i * neq]), neq * sizeof(double));
    // Verifica se la lettura è riuscita
    if (!inFile) {
      std::cerr << "Errore durante la lettura del file!" << std::endl;
      break;
    }

  }
  

  inFile.close();


  // // Stampare le condizioni casuali per verificarne il contenuto
  // for (int ii = 0; ii < num_condizioni; ++ii) {
  //   for (int jj= 0; jj < neq; ++jj) {
  //       std::cout << std::setprecision(3) << condizioni[ii * neq + jj] << " ";
  //   }
  //   std::cout << std::endl;
  // }
}


int save_condizioni_iniziali(int num_catene)
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
  int pd = 0;
  for (int ii = 0; ii<num_catene; ++ii){
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
      // if(h%(int(step)/100) == 0){
      //   std::cout << progress << "% executed" << std::endl; // Show the progress of simulation 
      //   progress++;
      // }
      RK4Step(t, X, Chain1, dt,neq);   // integration of the function
      // RK4Step(t, X, AlfaBeta_corrected, dt,neq);   // integration of the function
      t += dt; 
    }
    for(h=step - no_step; h<=step; h++){  
      
      // if(h%(int(step)/100) == 0){
      //   std::cout << progress << "% executed" << std::endl; // Show the progress of simulation 
      //   progress++;
      // }  
      RK4Step(t, X, Chain1, dt,neq);   // integration of the function
      // RK4Step(t, X, AlfaBeta_corrected, dt,neq);   // integration of the function
      t += dt; 
      if (drand48()<0.0001){
        condizioni.push_back(X);
        pd++;
      }
   
    }
  }

  int dimension_condizioni;
  dimension_condizioni = condizioni.size();

  std::ofstream outFile("condizioni.bin", std::ios::binary);
  if (!outFile) {
    std::cerr << "Errore nell'apertura del file per scrittura!" << std::endl;
    return 1;
  }

  // Scrivere le condizioni nel file, una per volta
  std::cout << "\nNumero condizioni salvate: " << pd << std::endl;
  for (const auto& condizione : condizioni) {
      outFile.write(reinterpret_cast<const char*>(condizione.data()), neq * sizeof(double));
  }
  outFile.close();
  if (!outFile.good()) {
    std::cerr << "Errore durante la scrittura dei dati nel file!" << std::endl;
    return 1;
  }

  std::cout << "Condizioni scritte con successo in formato binario." << std::endl;


  return 0;
}