#include "ode_solvers.h"
#include "ode_func.h"
#include "utils.h"
#include "param.h"

#include <random>
#include <mpi.h>
#include <fstream>
#include <sstream>
#include <cstring> // Per strcmp
#include <iterator> 

int main(int argc, char** argv) {
  MPI_Init(&argc, &argv);

  int rank, size;
  MPI_Comm mpicomm = MPI_COMM_WORLD;
  MPI_Comm_rank(mpicomm, &rank);
  MPI_Comm_size(mpicomm, &size);

  double start_time = MPI_Wtime();  // Start timer for ETA calculation

  if (argc > 1) {
    //lettura parmetri
    p.read (argv[1]);
    if (argc > 2) {
      //scrittura parametri
      p.write (argv[2]);
    }
  }
  
  // varibili lette da file 
  int dim = p.iparams["dim"];
  int N = p.iparams["N"];
  int step = p.iparams["step"];
  int catene_CI = p.iparams["catene_CI"];  // numero di catene per generare CI
  int catene_scelte = p.iparams["catene_scelte"];  // numero di catene
  bool save_conditions = p.bparams["save_conditions"];
  bool time_average = p.bparams["time_average"];
  double dt = p.dparams["dt"];
  double ergodic_mean;
  ////////////////////////////////////
  int  neq = (N+2)*2*dim + 2;
  std::vector<double> X(neq);
  long int h;
  std::vector<double> X_tot;
  
  // genera le condizioni iniziali
  if (rank == 0 && save_conditions){
    std::cout << "\nSalvatggio condizioni iniziali su " << catene_CI << " catene";
  	save_condizioni_iniziali(catene_CI);
  }
  if (rank == 0 && time_average){
    std::cout << "\nTime average ";
  	ergodic_mean=compute_mean(catene_CI);
    std::ostringstream filename;
    
    filename << p.sparams["dir"] << "/ergodic_N_" << N << "_Tr_" << p.dparams["Tr"] << ".dat";
    // Open file and save ergodic_mean
    std::ofstream outfile(filename.str());
    
    if (outfile.is_open()) {
        outfile << "Ergodic mean: " << ergodic_mean << std::endl;
        outfile.close();
        std::cout << "Ergodic mean saved in file: " << filename.str() << std::endl;
    } else {
        std::cerr << "Error opening file: " << filename.str() << std::endl;
    }
  // Creazione del nome del file con N e Tr
    
     
  }

if (!time_average){
  std::ofstream fdata;
  // Solo il processo 0 legge il file binario
  std::ostringstream file_name;
  // Creazione del nome del file con N e Tr
  file_name << p.sparams["dir"] << "/ttcf_mil_N_" << N << "_Tr_" << p.dparams["Tr"] << ".dat";

  
  if (rank == 0) {
  	std::cout << "\nnumero di catene scelto: " << catene_scelte <<std::endl<<std::endl;
    read_conditions(X_tot, catene_scelte, neq);
    system("mkdir -p single_data");
    fdata.open(file_name.str(), std::ios::out | std::ios::trunc);
    fdata  << std::setiosflags(std::ios::scientific); 
    fdata << std::setprecision(4);
    if (!fdata.is_open()) {
        std::cerr << "Error: Could not open file " << file_name.str() << std::endl;
    }
  }
  MPI_Barrier (mpicomm);
  // Determinare quante condizioni iniziali deve gestire ogni processo
  int conditiozioni_per_processo = catene_scelte / size;
  int remainder = catene_scelte % size;
  std::vector<int> sendcounts(size, conditiozioni_per_processo * neq);
  std::vector<int> displs(size, 0);

  for (int i = 0; i < remainder; ++i) {
      sendcounts[i] += neq;
  }

  for (int i = 1; i < size; ++i) {
    displs[i] = displs[i - 1] + sendcounts[i - 1];
  }

  int local_conditions = sendcounts[rank] / neq;
  std::cout << "Rank " << rank << " ha ricevuto " << local_conditions << " catene\n";

  // Vettore locale piatto per ciascun processo
  std::vector<double> X_local_flat(sendcounts[rank]);

  // Distribuzione dei dati appiattiti
  MPI_Scatterv(
    rank == 0 ? X_tot.data() : nullptr, sendcounts.data(), displs.data(), MPI_DOUBLE,
    X_local_flat.data(), sendcounts[rank], MPI_DOUBLE, 0, MPI_COMM_WORLD
  );

  //pulizia del vettore X_tot
  std::vector<double>().swap(X_tot);

  // Riorganizzare i dati ricevuti in `std::vector<std::vector<double>>`
  std::vector<std::vector<double>> X_local(local_conditions, std::vector<double>(neq));

  for (int i = 0; i < local_conditions; ++i) {
    std::copy(X_local_flat.begin() + i * neq, X_local_flat.begin() + (i + 1) * neq, X_local[i].begin());
  }

  //pulizia del vettore X_local_flat
  std::vector<double>().swap(X_local_flat);


  std::vector<double> omega_vec (X_local.size());
  double ttcf_mean = 0;
  double ttcf_mean_prev = 0;
  double ttcf_mean_integral = 0;

  double obs_mean = 0;
  double obs_mean_prev = 0;
  double obs_mean_integral = 0;
  double T_init=1.0;
  double omega_mean = 0;

  for (int i = 0; i < X_local.size(); ++i) {
    omega_vec[i] = omega_0(X_local[i], T_init); 
    ttcf_mean_prev += TTCF(observable_bulk, omega_vec[i],X_local[i], T_init);
    obs_mean_prev += observable_bulk(X_local[i]);
    omega_mean += omega_vec[i];
  }


  MPI_Reduce(rank == 0 ? MPI_IN_PLACE : &ttcf_mean_prev, &ttcf_mean_prev, 1, MPI_DOUBLE, MPI_SUM, 0, mpicomm);
  MPI_Reduce(rank == 0 ? MPI_IN_PLACE : &omega_mean, &omega_mean, 1, MPI_DOUBLE, MPI_SUM, 0, mpicomm);
  MPI_Reduce(rank == 0 ? MPI_IN_PLACE : &obs_mean_prev, &obs_mean_prev, 1, MPI_DOUBLE, MPI_SUM, 0, mpicomm);

  if (rank == 0){
    ttcf_mean_prev = ttcf_mean_prev/catene_scelte;
    obs_mean_prev = obs_mean_prev/catene_scelte;
    omega_mean = omega_mean/catene_scelte;
    std::cout << "Media della omega: " << omega_mean << std::endl;
    std::cout << "Media osservabile: " << obs_mean_prev << std::endl;
  }


  // Evolvere le condizioni iniziali
  double t = 0.0;
  std::vector<double> t_vec (X_local.size(),0);


  // Imposta la dimensione del buffer
  double max_memory = 1.e8;  // 100 MB per vector buff
  int max_procs = 64; // # di processori su nodo a memoria condivisa
  size_t buffer_size = max_memory/(max_procs*10);
  if (rank==0){
    buffer_size = buffer_size < step ? buffer_size : step;
    std::cout <<"\nBuffer size: " << buffer_size << std::endl;
  }
  std::vector<std::string> buffer;  // Buffer in RAM per accumulare i dati
  buffer.reserve (buffer_size);
  
  for ( h = 1; h <= step; ++h) {
    ttcf_mean = 0;
    for (int i = 0; i < X_local.size(); ++i) {
      RK4Step(t_vec[i], X_local[i], betaFPUT, dt,neq);   // integration of the function
      // RK4Step(t, X_local[i], AlfaBetaFPUT, dt,neq);   // integration of the function
      t_vec[i] += dt;
      ttcf_mean += TTCF(observable_bulk, omega_vec[i],X_local[i], T_init);
    }

    MPI_Reduce(rank == 0 ? MPI_IN_PLACE : &ttcf_mean, &ttcf_mean, 1, MPI_DOUBLE, MPI_SUM, 0, mpicomm);

    if (rank == 0){
      ttcf_mean = ttcf_mean/catene_scelte;
      // // integrare sul tempo trapezi int_a^b (f(a)+ f(b))*(b-a)*0.5
      ttcf_mean_integral += (ttcf_mean + ttcf_mean_prev)*dt*0.5;
      std::ostringstream oss;
      oss << ttcf_mean << " " << ttcf_mean_integral << std::endl;  
      buffer.push_back(oss.str());
      // Quando il buffer è pieno, scrivi tutto sul file
      if (buffer.size() >= buffer_size) {
          std::copy(buffer.begin(), buffer.end(), std::ostream_iterator<std::string>(fdata, ""));
          buffer.clear();  // Svuota il buffer dopo la scrittura
      }
      ttcf_mean_prev = ttcf_mean;
    }
       
  }

  if (rank == 0){
    if (buffer.size() > 0){
      std::copy(buffer.begin(), buffer.end(), std::ostream_iterator<std::string>(fdata, ""));
      buffer.clear();  // Svuota il buffer dopo la scrittura
    }
    fdata.close();
  }

  double end_time = MPI_Wtime();  // end timer for ETA calculation
  if (rank == 0)
    std::cout <<"\nTotal simulation time: " << end_time - start_time << std::endl;
}
  MPI_Finalize();
	return 0;
}