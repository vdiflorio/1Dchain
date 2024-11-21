#include "ode_solvers.h"
#include "ode_func.h"
#include "utils.h"

#include <random>
#include <mpi.h>
#include <fstream>

int main(int argc, char **argv) {
	
	MPI_Init(&argc, &argv);

  int rank, size;
  MPI_Comm mpicomm = MPI_COMM_WORLD;
  MPI_Comm_rank(mpicomm, &rank);
  MPI_Comm_size(mpicomm, &size);



  int      neq = (N+2)*2*dim + 2;
  std::vector<double> X(neq);
  long int step = 50000;
  long int h;
  std::vector<double> X_tot;
  int num_catene = 1;  // numero di catene per generare CI
  int num_condizioni = 10;  // numero di catene

  // genera le condizioni iniziali
  if (rank == 0 && false){
    std::cout << "\nSalvatggio condizioni iniziali su " << num_catene << " catene";
  	save_condizioni_iniziali(num_catene);
  }
  std::ofstream fdata;
  // Solo il processo 0 legge il file binario
  if (rank == 0) {
  	std::cout << "\nnumero di catene scelto: " << num_condizioni <<std::endl<<std::endl;
    read_conditions(X_tot, num_condizioni, neq);
    fdata.open("ttcf.dat");
    fdata  << std::setiosflags(std::ios::scientific); 
    fdata << std::setprecision(4);
  }
  MPI_Barrier (mpicomm);
  // Determinare quante condizioni iniziali deve gestire ogni processo
  int conditiozioni_per_processo = num_condizioni / size;
  int remainder = num_condizioni % size;
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

  // // Debug: stampa dei dati ricevuti per ciascun rank
  // for (const auto& cond : X_local) {
  //   std::cout << "Rank " << rank << " :  ";
  //   for (double val : cond) {
  //     std::cout << std::setprecision(3)<< val << " ";
  //   }
  //   std::cout << "\n";
  // }

  std::vector<double> omega_vec (X_local.size());
  double ttcf_mean = 0;
  double ttcf_mean_prev = 0;
  double ttcf_mean_integral = 0;

  double obs_mean = 0;
  double obs_mean_prev = 0;
  double obs_mean_integral = 0;

  double omega_mean = 0;

  for (int i = 0; i < X_local.size(); ++i) {
    omega_vec[i] = omega_0(X_local[i], Tl); 
    ttcf_mean_prev += TTCF(observable, omega_vec[i],X_local[i], Tl);
    obs_mean_prev += observable (X_local[i]);
    omega_mean += omega_vec[i];
  }

  MPI_Barrier (mpicomm);
  if (rank == 0) {
    MPI_Reduce (MPI_IN_PLACE, &ttcf_mean_prev, 1, MPI_DOUBLE, MPI_SUM, 0,
                mpicomm);
    MPI_Reduce (MPI_IN_PLACE, &obs_mean_prev, 1, MPI_DOUBLE, MPI_SUM, 0,
                mpicomm);
    MPI_Reduce (MPI_IN_PLACE, &omega_mean, 1, MPI_DOUBLE, MPI_SUM, 0,
                mpicomm);
  } else {
    MPI_Reduce (&ttcf_mean_prev, &ttcf_mean_prev, 1, MPI_DOUBLE, MPI_SUM, 0,
                mpicomm);
    MPI_Reduce (&obs_mean_prev, &obs_mean_prev, 1, MPI_DOUBLE, MPI_SUM, 0,
                mpicomm);
    MPI_Reduce (&omega_mean, &omega_mean, 1, MPI_DOUBLE, MPI_SUM, 0,
                mpicomm);
  }
  if (rank == 0){
    ttcf_mean_prev = ttcf_mean_prev/num_condizioni;
    obs_mean_prev = obs_mean_prev/num_condizioni;
    omega_mean = omega_mean/num_condizioni;
    std::cout << "Media della omega: " << omega_mean << std::endl;
  }


  // Evolvere le condizioni iniziali
  double t = 0.0;
  std::vector<double> t_vec (X_local.size(),0);
  double dt = 1.e-3;
  for ( h = 1; h <= step; ++h) {
    ttcf_mean = 0;
    obs_mean = 0;
    for (int i = 0; i < X_local.size(); ++i) {
      RK4Step(t_vec[i], X_local[i], Chain1, dt,neq);   // integration of the function
      // RK4Step(t, X_local[i], AlfaBeta, dt,neq);   // integration of the function
      t_vec[i] += dt;
      ttcf_mean += TTCF(observable, omega_vec[i],X_local[i], Tl);
      obs_mean += observable(X_local[i]);
    }
    if (rank == 0) {
      MPI_Reduce (MPI_IN_PLACE, &ttcf_mean, 1, MPI_DOUBLE, MPI_SUM, 0,
                  mpicomm);
      MPI_Reduce (MPI_IN_PLACE, &obs_mean, 1, MPI_DOUBLE, MPI_SUM, 0,
                  mpicomm);

    } else {
      MPI_Reduce (&ttcf_mean, &ttcf_mean, 1, MPI_DOUBLE, MPI_SUM, 0,
                  mpicomm);
      MPI_Reduce (&obs_mean, &obs_mean, 1, MPI_DOUBLE, MPI_SUM, 0,
                  mpicomm);
    }
    if (rank == 0){
      ttcf_mean = ttcf_mean/num_condizioni;
      obs_mean = obs_mean/num_condizioni;
      // integrare sul tempo trapezi int_a^b (f(a)+ f(b))*(b-a)*0.5
      ttcf_mean_integral += (ttcf_mean + ttcf_mean_prev)*dt*0.5;
      obs_mean_integral += (obs_mean + obs_mean_prev)*dt*0.5;
      // salva su file ogni ttcf_mean_integral
      // salva su file ttcf_mean
      fdata << t_vec[0] << " "<<ttcf_mean 
                        << "  " << ttcf_mean_integral
                        << "  " << ttcf_mean_integral - obs_mean_integral*omega_mean << std::endl;
      ttcf_mean_prev = ttcf_mean;
      obs_mean_prev = obs_mean;
    }
    // stop 
    MPI_Barrier (mpicomm);
  }

  if (rank == 0)
    fdata.close ();

  MPI_Finalize();
	return 0;
}