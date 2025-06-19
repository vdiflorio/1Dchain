// File for auxiliary functions (ex. Omega)
#include "ode_func.h"
#include "ode_solvers.h"
#include "param.h"

#include "utils.h"
#include <cstdint>
#include <sstream>
#include <random>
#include <chrono>
#include <sys/mman.h>
#include <fcntl.h>
#include <unistd.h>

double omega_0 (std::vector<double> &Y, double T)
{
  // Function to compute the dissipation function wrt the initial density

  double Tl = p.dparams["Tl"];
  double Tr = p.dparams["Tr"];
  double m = p.dparams["m"];
  int dim = p.iparams["dim"];
  int N = p.iparams["N"];

  int j,k;
  k = 2*dim;
  double omega=0.0;

  for (j =0; j<dim; j++) {
    omega+=Y[k* (N+2)+1]*Y[k*N+dim+j]*Y[k*N+dim+j]* (1.0/Tr-1.0/T)/m +
           Y[k* (N+2)]*Y[k+dim+j]*Y[k+dim+j]* (1.0/Tl-1.0/T)/m;
  }

  return omega;
}



double TTCF (std::function<double (std::vector<double>&)> obs, double omega,std::vector<double> &Y,double T)
{
  // Function to compute the TTCF
  return obs (Y)*omega;
}

double observable_tot (std::vector<double> &Y)
{
  // Function to compute any observable, that is function of the phase space state
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

  int k = 2*dim;
  int i,j,n;
  double r1 =0.0, r2 =0.0;
  // static std::vector <double> F;
  // F.assign (N, 0.0);
  double flux = 0;
  double F_l, F_r;
  double Psi1, Psi2; // DEFINISCO LE VARIABILI DEL TERMOSTATO
  double kin1 =0.0, kinN=0.0; // energia kin della prima e ultima
  // particlella mobile
  Psi1 = Y[k* (N+2)];
  Psi2 = Y[k* (N+2)+1];

  //eq. for termostat
  F_l = (kin1/Tl -1.0)/ (thetaL*thetaL); //left
  F_r = (kinN/Tr -1.0)/ (thetaR*thetaR); //right


  ////////////////////////////////////

  ////////////////////////////////////

  i = 1; //la prima particella che si può muovere
  n = k*i;
  r1 = (Y[n+k] - Y[n] - a);
  r2 = (Y[n] - Y[n-k] - a);
  // eq. for momentum
  flux = (chi* (Y[n+k] + Y[n-k] - 2.0*Y[n] ) +
          Alpha* (r1*r1 - r2*r2) +
          bet* (r1*r1*r1 - r2*r2*r2) - Psi1*F_l)*Y[n+dim]/m;

  // le altre N-2 masse

  for (i = 2; i < N; i++) {
    n = k*i;
    r2 =r1;
    r1 = (Y[n+k] - Y[n] - a);
    // eq. for momentum
    flux += (chi* (Y[n+k] + Y[n-k] - 2.0*Y[n] ) +
             Alpha* (r1*r1 - r2*r2) +
             bet* (r1*r1*r1 - r2*r2*r2))*Y[n+dim]/m;
  }

  //l' ultima massa mobile
  i = N;
  n = k*i;
  r2 = r1;
  r1 = (Y[n+k] - Y[n] - a);

  // eq. for momentum
  flux += (chi* (Y[n+k] + Y[n-k] - 2.0*Y[n] ) +
           Alpha* (r1*r1 - r2*r2) +
           bet* (r1*r1*r1 - r2*r2*r2) - Psi2*F_r)* Y[n+dim]/m;



  return flux/N;
}




double observable_bulk (std::vector<double> &Y)
{
  // Function to compute any observable, that is function of the phase space state
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

  int k = 2*dim;
  int i,j,n;
  double r1 =0.0;
  double flux = 0;
  int bd_paticle = N*0.15;



  for (i = bd_paticle; i <= N-bd_paticle; i++) {
    n = k*i;

    r1 = (Y[n+k] - Y[n] - a);
    // eq. for momentum
    flux += (chi* (r1) +
             Alpha* (r1*r1) +
             bet* (r1*r1*r1))*Y[n+dim]/m;

  }

  //std::cout<<(N-bd_paticle*2.0)<<"  "<<flux<<" "<<flux/(N-bd_paticle*2.0)<<std::endl;
  return flux/ (N-bd_paticle*2.0);
}

double observable (std::vector<double> &Y)
{
  // Function to compute any observable, that is function of the phase space state
  int dim = p.iparams["dim"];
  int N = p.iparams["N"];

  int k = 2*dim;
  return -Y[k* (N+2)]*Y[k+dim]*Y[k+dim]+Y[k* (N+2)+1]*Y[k*N+dim]*Y[k*N+dim];
}


double dumb_observable (std::vector<double> &Y)
{
  double m = p.dparams["m"];
  double a = p.dparams["a"];
  double chi = p.dparams["chi"];
  double bet = p.dparams["beta"];
  double Alpha = p.dparams["alpha"];
  int dim = p.iparams["dim"];
  int N = p.iparams["N"];

  int k = 2*dim;
  int i,j,n;
  double r1 =0.0;
  double flux = 0;

  int bd_paticle = N*0.5;
  n = k* (bd_paticle);
  r1 = (Y[n+k] - Y[n] - a);
  // eq. for momentum
  flux = (chi* (r1) +
          Alpha* (r1*r1) +
          bet* (r1*r1*r1))*Y[n+dim]/m;


  return flux;
}




void read_conditions (std::vector<double>& condizioni, int num_condizioni, int neq)
{
  // Apri il file binario
  // Creazione del nome del file dinamico
  int dim = p.iparams["dim"];
  int N = p.iparams["N"];

  std::ostringstream name;
  name << "condizioni_" << N << ".bin";
  std::string filename = name.str();

  

  // Apertura del file con il nome dinamico
  int fd = open (filename.c_str(), O_RDONLY);

  if (fd == -1) {
    std::cerr << "Errore nell'apertura del file per lettura!" << std::endl;
    return;
  }

  // Ottieni la dimensione del file
  off_t file_size = lseek (fd, 0, SEEK_END);

  if (file_size == -1) {
    std::cerr << "Errore nell'ottenere la dimensione del file!" << std::endl;
    return;
  }

  // Mappa il file in memoria
  void* file_memory = mmap (NULL, file_size, PROT_READ, MAP_PRIVATE, fd, 0);

  if (file_memory == MAP_FAILED) {
    std::cerr << "Errore nel mappare il file in memoria!" << std::endl;
    return;
  }

  // Dimensione di una condizione in byte
  int dimensione_condizione = neq * sizeof (double);

  // Calcola il numero di condizioni nel file
  int numero_condizioni_tot = file_size / dimensione_condizione;

  if (num_condizioni > numero_condizioni_tot) {
    std::cerr << "Errore: il numero di condizioni richiesto eccede il numero di condizioni nel file." << std::endl;
    munmap (file_memory, file_size);
    close (fd);
    return;
  }

  // Prepara il vettore per le condizioni
  condizioni.reserve (num_condizioni * neq);

  int64_t num_selezioni = num_condizioni/ 2; // Numero di indici da selezionare

  // Genera un vettore con tutti gli indici
  std::vector<int64_t> all_indices (num_condizioni);
  std::iota (all_indices.begin(), all_indices.end(), 0); // Riempie con {0, 1, 2, ..., num_condizioni - 1}

  // Mescola il vettore
  std::random_device rd;
  std::mt19937 gen (rd());
  std::shuffle (all_indices.begin(), all_indices.end(), gen);

  // Seleziona i primi num_selezioni indici
  std::vector<int64_t> indices (all_indices.begin(), all_indices.begin() + num_selezioni);


  // Timer per calcolare la velocità
  auto start_time = std::chrono::high_resolution_clock::now();
  int read_count = 0;

  // Lettura dei dati
  for (int i = 0; i < num_condizioni/2; ++i) {
    // Calcolare l'offset per l'indice casuale
    std::streampos offset = indices[i] * dimensione_condizione;

    // Leggere direttamente dalla memoria mappata
    double* condizione_ptr = reinterpret_cast<double*> (static_cast<char*> (file_memory) + offset);
    std::copy (condizione_ptr, condizione_ptr + neq, condizioni.begin() + i * neq);

    read_count++;

    // Ogni 1000 letture, calcola e stampa la velocità
    if (read_count == 100000) {
      auto current_time = std::chrono::high_resolution_clock::now();
      std::chrono::duration<double> elapsed_time = current_time - start_time;
      double speed = 100000.0 / elapsed_time.count(); // letture al secondo
      std::cout << "Velocità di lettura: " << speed << " condizioni/secondo" << std::endl;

      // Reset del timer per il prossimo intervallo di 1000 letture
      start_time = current_time;

    }
  }

  // // Modifica delle condizioni per il secondo ciclo
  for (int k = num_condizioni / 2; k < num_condizioni; ++k) {
    for (int j = 0; j < neq; ++j) {
      if (j < neq - 2) {
        // Per gli indici 0 a neq-2, prendi i componenti dispari e moltiplica per -1
        if (j % 2 != 0) {
          condizioni[k * neq + j] = -condizioni[ (k - num_condizioni / 2) * neq + j];
        } else {
          // Per i componenti pari, copia semplicemente
          condizioni[k * neq + j] = condizioni[ (k - num_condizioni / 2) * neq + j];
        }
      } else {
        // Per gli ultimi due componenti, moltiplica per -1
        condizioni[k * neq + j] = -condizioni[ (k - num_condizioni / 2) * neq + j];
      }
    }
  }

  // Libera la memoria mappata e chiudi il file
  munmap (file_memory, file_size);
  close (fd);
}


int save_condizioni_iniziali (int num_catene)
{
  int dim = p.iparams["dim"];
  int N = p.iparams["N"];
  double a = p.dparams["a"];

  int neq = (N+2)*2*dim + 2;
  std::vector<double> X (neq);
  long int step = 50000000;
  long int h;
  long int no_step = 8000000;
  long int progress = 0;
  int k,i,j,l,n;
  double t, dt;
  double X0[dim]; //unità di spaziatura
  double X_eq[dim* (N+2)]; //vettore contenete posizioni di equilibrio meccanico (potrebbe essere cancellato)

  std::vector<std::vector<double>> condizioni;

  k = 2*dim;
  t = 0.0; //tempo zero
  dt = 1.e-3; //intervallo di integrazione

  X0[0] = a; // nodes only along x-axis

  // Nodes position at mechanical equilibrium (potrebbe non essere utile per il codice)
  for (j=0; j<=N+1; j++) {
    l = dim*j;

    for (i=0; i<dim; i++) {
      X_eq[l+i] = j*X0[i];
    }
  }

  int pd = 0;

  for (int ii = 0; ii<num_catene; ++ii) {
    std::cout << "\n\n CATENA NUMERO: " << ii+1 << std::endl << std::endl;

    // INIZILIZZAZIONE VETTORI
    for (i=0; i < neq; i++) {
      X[i] = 0.0;
    }

    X[k* (N+2)] = 1.0; //inizializzazione termostato di sinistra
    X[k* (N+2)+1] = 1.0; //inizializzazione termostato di destra
    ////////////////////////////////////////////////////////////////

    double alfa = 1.0; // fattore di stretching

    // CONDIZIONE INIZIALE per posizione e velocità
    srand48 (time (NULL)); // Initialize the sequence

    for (j= 1; j<=N; j++) {
      n=k*j;
      l=dim*j;

      for (i= 0; i<dim; ++i) {
        X[n+i]= X_eq[l+i]*alfa + drand48()*0.8 - 0.4; //posizione
        X[n+i+dim]= drand48()*0.2 - 0.1; //velocità
      }
    }

    X[ (N+1)*k] = X_eq[ (N+1)*dim]*alfa; //l'ultima particella è fissa



    //EVOLUZIONE SISTEMA
    for (h=1; h<=step - no_step; h++) {
      //RK4Step(t, X, betaFPUT, dt,neq);   // integration of the function
      RK4Step (t, X, LepriChain, dt,neq); // integration of the function
      t += dt;
    }

    for (h=step - no_step; h<=step; h++) {
      //RK4Step(t, X, betaFPUT, dt,neq);   // integration of the function
      RK4Step (t, X, LepriChain, dt,neq); // integration of the function
      t += dt;

      if (drand48()<0.005) {
        condizioni.push_back (X);
        pd++;
      }

    }
  }

  int dimension_condizioni;
  dimension_condizioni = condizioni.size();
  std::ostringstream name;
  name << p.sparams["dir_CI"] <<"/condizioni_" << N << ".bin";

  // Apertura del file con il nome dinamico
  std::ofstream outFile (name.str(), std::ios::binary);

  if (!outFile) {
    std::cerr << "Errore nell'apertura del file per scrittura!" << std::endl;
    return 1;
  }

  // Scrivere le condizioni nel file, una per volta
  std::cout << "\nNumero condizioni salvate: " << pd << std::endl;

  for (const auto& condizione : condizioni) {
    outFile.write (reinterpret_cast<const char*> (condizione.data()), neq * sizeof (double));
  }

  outFile.close();

  if (!outFile.good()) {
    std::cerr << "Errore durante la scrittura dei dati nel file!" << std::endl;
    return 1;
  }

  std::cout << "Condizioni scritte con successo in formato binario." << std::endl;


  return 0;
}


void compute_mean ()
{
  int dim = p.iparams["dim"];
  int N = p.iparams["N"];

  double m = p.dparams["m"];
  double a = p.dparams["a"];
  double chi = p.dparams["chi"];
  double bet = p.dparams["beta"];
  double Alpha = p.dparams["alpha"];

  int neq = (N+2)*2*dim + 2;
  std::vector<double> X (neq);
  long int step = 1000000000;
  long int h;
  long int no_step = 1000000000-100000000;
  long int progress = 0;
  int k,i,j,l,n;
  double t, dt;
  double X0[dim]; //unità di spaziatura
  double X_eq[dim* (N+2)]; //vettore contenete posizioni di equilibrio meccanico (potrebbe essere cancellato)
  double r1;
  std::vector<std::vector<double>> condizioni;

  k = 2*dim;
  t = 0.0; //tempo zero
  dt = 1.e-2; //intervallo di integrazione

  X0[0] = a; // nodes only along x-axis

  // Nodes position at mechanical equilibrium (potrebbe non essere utile per il codice)
  for (j=0; j<=N+1; j++) {
    l = dim*j;

    for (i=0; i<dim; i++) {
      X_eq[l+i] = j*X0[i];
    }
  }

  int bd_paticle = N*0.15;
  std::vector<double> cum_mean (N-bd_paticle*2);
  std::vector<double> Temperature (N);

  

  // INIZILIZZAZIONE VETTORI
  for (i=0; i < neq; i++) {
    X[i] = 0.0;
  }

  X[k* (N+2)] = 1.0; //inizializzazione termostato di sinistra
  X[k* (N+2)+1] = 1.0; //inizializzazione termostato di destra
  ////////////////////////////////////////////////////////////////

  double alfa = 1.0; // fattore di stretching

  // CONDIZIONE INIZIALE per posizione e velocità
  srand48 (time (NULL)); // Initialize the sequence

  for (j= 1; j<=N; j++) {
    n=k*j;
    l=dim*j;

    for (i= 0; i<dim; ++i) {
      X[n+i]= X_eq[l+i]*alfa + drand48()*0.8 - 0.4; //posizione
      X[n+i+dim]= drand48() + 0.5; //velocità
    }
  }

  X[ (N+1)*k] = X_eq[ (N+1)*dim]*alfa; //l'ultima particella è fissa

  std::ostringstream convergence_mean;
  convergence_mean << p.sparams["dir"] << "/ergodic_conv_N_" << N << "_Tr_" << p.dparams["Tr"] << ".dat";
  // Open file and save ergodic_mean
  std::ofstream outfile_conv (convergence_mean.str());
  int n_conv = 1000;
  std::vector<double> vect_conv;
  double cum_mean_tmp=0.0;
  int ni;
  //EVOLUZIONE SISTEMA
  for (h=1; h<=step - no_step; h++) {
    // RK4Step(t, X, betaFPUT, dt,neq);   // integration of the function
    RK4Step (t, X, AlfaBetaFPUT, dt,neq); // integration of the function
    t += dt;
    ni = k*int(N*0.5);
    r1 = (X[ni+k] - X[ni] - a);
    cum_mean_tmp += (chi* (r1) +
                    Alpha* (r1*r1) +
                    bet* (r1*r1*r1))*X[ni+dim]/m;
    if (h%n_conv == 0){
      vect_conv.push_back (cum_mean_tmp/h);
    }
  }

  for (h=step - no_step; h<=step; h++) {
    // RK4Step(t, X, betaFPUT, dt,neq);   // integration of the function
    RK4Step (t, X, AlfaBetaFPUT, dt,neq); // integration of the function
    t += dt;

    ni = k*int(N*0.5);
    r1 = (X[ni+k] - X[ni] - a);
    cum_mean_tmp += (chi* (r1) +
                    Alpha* (r1*r1) +
                    bet* (r1*r1*r1))*X[ni+dim]/m;
    if (h%n_conv == 0){
      vect_conv.push_back (cum_mean_tmp/h);
    }
    
    for (i = bd_paticle; i < N-bd_paticle; i++) {
      n = k*i;
      r1 = (X[n+k] - X[n] - a);
      cum_mean[i-bd_paticle]+= (chi* (r1) +
                                Alpha* (r1*r1) +
                                bet* (r1*r1*r1))*X[n+dim]/m;
    }
    for (i = 1; i <= N; i++) {
      n = k*i;
      Temperature[i-1]+= X[n+dim]*X[n+dim]/m;
    }
  }
  

  std::ostringstream filename;

  filename << p.sparams["dir"] << "/ergodic_N_" << N << "_Tr_" << p.dparams["Tr"] << ".dat";
  // Open file and save ergodic_mean
  std::ofstream outfile (filename.str());

  if (outfile.is_open()) {
    for (i = 0; i < N-bd_paticle*2; i++) {
      outfile << "Ergodic mean: " << cum_mean[i]/no_step << std::endl;
    }
    outfile <<std::endl;
    outfile << "@@"<<std::endl;
    outfile <<std::endl;
    for (i = 0; i < N; i++) {
      outfile << "Temperature: " << Temperature[i]/no_step << std::endl;
    }
  } else {
    std::cerr << "Error opening file: " << filename.str() << std::endl;
  }

  if (outfile_conv.is_open ()){
    for (i = 0; i <vect_conv.size (); i++) {
      outfile_conv << n_conv*(i+1) << "  " << vect_conv[i] << std::endl;
    }
  } else {
    std::cerr << "Error opening file: " << convergence_mean.str() << std::endl;
  }

  outfile.close();
  std::cout << "Ergodic mean saved in file: " << filename.str() << std::endl;
  outfile_conv.close();
  std::cout << "Ergodic mean convergence saved in file: " << convergence_mean.str() << std::endl;

}