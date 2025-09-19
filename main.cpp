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

#include <chrono>


int main (int argc, char** argv)
{
  MPI_Init (&argc, &argv);

  int rank, size;
  MPI_Comm mpicomm = MPI_COMM_WORLD;
  MPI_Comm_rank (mpicomm, &rank);
  MPI_Comm_size (mpicomm, &size);

  double start_time = MPI_Wtime(); // Start timer for ETA calculation

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
  int catene_CI = p.iparams["catene_CI"]; // numero di catene per generare CI
  int catene_scelte = p.iparams["catene_scelte"]; // numero di catene
  bool save_conditions = p.bparams["save_conditions"];
  bool time_average = p.bparams["time_average"];
  double dt = p.dparams["dt"];
  double ergodic_mean;
  ////////////////////////////////////
  int neq = (N+2)*2*dim + 2;
  std::vector<double> X (neq);
  long int h;
  std::vector<double> X_tot;

  // genera le condizioni iniziali
  if (rank == 0 && save_conditions) {
    std::cout << "\nSalvatggio condizioni iniziali su " << catene_CI << " catene";
    save_condizioni_iniziali (catene_CI);
  }

  if (rank == 0 && time_average) {
    std::cout << "\nTime average ";
    compute_mean ();
  }

  if (!time_average) {
    std::ofstream fdata;
    // Solo il processo 0 legge il file binario
    std::ostringstream file_name;
    // Creazione del nome del file con N e Tr
    file_name << p.sparams["dir"] << "/ttcf_mil_N_" << N << "_Tr_" << p.dparams["Tr"] << ".dat";

    // Numero massimo da leggere dal file
    const int MAX_FROM_FILE = 500000;
    // Quante condizioni effettive servono
    int total_conditions = catene_scelte;
    // Decidi quante leggere
    int read_conditions_num = std::min(total_conditions, MAX_FROM_FILE);

    if (rank == 0) {
      std::cout << "\nnumero di catene scelto: " << catene_scelte <<std::endl<<std::endl;
      if (read_conditions_num >= MAX_FROM_FILE)
        read_conditions_subset(X_tot, neq);   // legge 500k
      else
        read_conditions(X_tot, read_conditions_num, neq); // legge meno

      system ("mkdir -p single_data");
      fdata.open (file_name.str(), std::ios::out | std::ios::trunc);
      fdata << std::setiosflags (std::ios::scientific);
      fdata << std::setprecision (4);

      if (!fdata.is_open()) {
        std::cerr << "Error: Could not open file " << file_name.str() << std::endl;
      }
    }

    MPI_Barrier (mpicomm);
    // --- Distribuzione delle condizioni lette ---
    int base_conditions_per_proc = read_conditions_num / size;
    int remainder = read_conditions_num % size;
    
    std::vector<int> sendcounts(size, base_conditions_per_proc * neq);
    for (int i = 0; i < remainder; ++i) sendcounts[i] += neq;

    std::vector<int> displs(size, 0);
    for (int i = 1; i < size; ++i) displs[i] = displs[i-1] + sendcounts[i-1];

    int local_conditions_read = sendcounts[rank] / neq;
    std::cout << "Rank " << rank << " ha letto e ricevuto " << local_conditions_read << " catene\n";

    // Vettore locale piatto per ciascun processo
    std::vector<double> X_local_flat (sendcounts[rank]);

    // Distribuzione dei dati appiattiti
    MPI_Scatterv (
      rank == 0 ? X_tot.data() : nullptr, sendcounts.data(), displs.data(), MPI_DOUBLE,
      X_local_flat.data(), sendcounts[rank], MPI_DOUBLE, 0, MPI_COMM_WORLD
    );

    //pulizia del vettore X_tot
    std::vector<double>().swap (X_tot);

    // Riorganizzare i dati ricevuti in `std::vector<std::vector<double>>`
    std::vector<std::vector<double>> X_local (local_conditions_read, std::vector<double> (neq));

    for (int i = 0; i < local_conditions_read; ++i) 
      std::copy (X_local_flat.begin() + i * neq, X_local_flat.begin() + (i + 1) * neq, X_local[i].begin());
    
    // --- Ora generiamo condizioni extra se necessario ---
    int extra_needed = (total_conditions - read_conditions_num)/2;
    int extra_per_proc = extra_needed / size;
    if (rank < extra_needed % size) extra_per_proc++;

    // Generazione condizioni addizionali con simmetria
    for (int i = 0; i < extra_per_proc; ++i) {
      std::vector<double> cond(neq);
      generate_condition(cond);
      // salva condizione originale
      X_local.push_back(cond);
      // crea e salva la condizione simmetrica
      std::vector<double> cond_sym(neq);
      for (int j = 0; j < neq; ++j) {
        if (j < neq - 2) {
          cond_sym[j] = (j % 2 != 0 ? -cond[j] : cond[j]);
        } else {
          cond_sym[j] = -cond[j];
        }
      }
      X_local.push_back(cond_sym);
    }

    // Ora X_local contiene sia quelle lette che quelle generate
    std::cout << "Rank " << rank << " ha " << X_local.size() << " condizioni totali\n";
    
      //pulizia del vettore X_local_flat
    std::vector<double>().swap (X_local_flat);



    std::vector<double> omega_vec (X_local.size());
    double ttcf_mean = 0;
    double ttcf_mean_prev = 0;
    double ttcf_mean_integral = 0;

    double obs_mean = 0;
    double obs_mean_prev = 0;
    double obs_mean_integral = 0;
    double T_init=1.0;
    double omega_mean = 0;
    double loc = N*0.5;

    for (int i = 0; i < X_local.size(); ++i) {
      omega_vec[i] = omega_0 (X_local[i], T_init);
      ttcf_mean_prev += TTCF (observable_bulk, omega_vec[i],X_local[i], T_init);
      obs_mean_prev += observable_bulk (X_local[i]);
      omega_mean += omega_vec[i];
    }


    MPI_Reduce (rank == 0 ? MPI_IN_PLACE : &ttcf_mean_prev, &ttcf_mean_prev, 1, MPI_DOUBLE, MPI_SUM, 0, mpicomm);
    MPI_Reduce (rank == 0 ? MPI_IN_PLACE : &omega_mean, &omega_mean, 1, MPI_DOUBLE, MPI_SUM, 0, mpicomm);
    MPI_Reduce (rank == 0 ? MPI_IN_PLACE : &obs_mean_prev, &obs_mean_prev, 1, MPI_DOUBLE, MPI_SUM, 0, mpicomm);

    if (rank == 0) {
      ttcf_mean_prev = ttcf_mean_prev/catene_scelte;
      obs_mean_prev = obs_mean_prev/catene_scelte;
      omega_mean = omega_mean/catene_scelte;
      std::cout << "Media della omega: " << omega_mean << std::endl;
      std::cout << "Media osservabile: " << obs_mean_prev << std::endl;
    }


    // Evolvere le condizioni iniziali
  double t = 0.0;
  std::vector<double> t_vec(X_local.size(), 0.0);

  // Preallocate local time series of TTCF means
  std::vector<double> local_ttcf(step, 0.0);

  // Vector to store omega values
  for (int i = 0; i < X_local.size(); ++i) {
    omega_vec[i] = omega_0(X_local[i], T_init);
  }

  // Main time loop
  for (int h = 0; h < step; ++h) {
    double ttcf_mean = 0.0;

    for (int i = 0; i < X_local.size(); ++i) {
      RK4Step(t, X_local[i], AlfaBetaFPUT, dt, neq);
      t_vec[i] += dt;
      ttcf_mean += TTCF(observable_bulk, omega_vec[i], X_local[i], T_init);
    }

    // Store per-step average in the local buffer
    local_ttcf[h] = ttcf_mean / X_local.size();
    t += dt;
  }

  // Now reduce to global_ttcf vector
  std::vector<double> global_ttcf(step, 0.0);

  MPI_Reduce(
    local_ttcf.data(),
    global_ttcf.data(),
    step,
    MPI_DOUBLE,
    MPI_SUM,
    0,
    mpicomm
  );

  // Rank 0: finalize average and write to file
  std::ostringstream output_buffer;
  output_buffer << std::setprecision(8) << std::scientific;

  for (int h = 0; h < step; ++h) {
    global_ttcf[h] /= size;
    output_buffer << global_ttcf[h] << "\n";
  }

  fdata << output_buffer.str(); // One single write
  fdata.close();
  }
  MPI_Finalize();
  return 0;
}