#include "ode_solvers.h"
#include "ode_func.h"
#include "utils.h"

#include <random>
#include <mpi.h>

int main(int argc, char **argv) {
	
	MPI_Init(&argc, &argv);

  int rank, size;
  MPI_Comm mpicomm = MPI_COMM_WORLD;
  MPI_Comm_rank(mpicomm, &rank);
  MPI_Comm_size(mpicomm, &size);



  int      neq = (N+2)*2*dim + 2;
  std::vector<double> X(neq);
  long int step = 8000000;
  long int h;
  std::vector<double> X_tot;
  int num_catene = 1;  // numero di catene per generare CI
  int num_condizioni = 9;  // numero di catene

  // genera le condizioni iniziali
  if (rank == 0 && true){
    std::cout << "\nSalvatggio condizioni iniziali su " << num_catene << "catene";
  	save_condizioni_iniziali(num_catene);
  }

  // Solo il processo 0 legge il file binario
  if (rank == 0) {
  	std::cout << "\nnumero di catene scelto: " << num_condizioni <<std::endl<<std::endl;
    read_conditions(X_tot, num_condizioni, neq);
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


  // Evolvere le condizioni iniziali
  double t = 0.0;
  double dt = 1.e-3;
  for ( h = 1; h <= step; ++h) {
    for (int i = 0; i < X_local.size(); ++i) {
      RK4Step(t, X_local[i], Chain1, dt,neq);   // integration of the function
      // RK4Step(t, X_local[i], AlfaBeta, dt,neq);   // integration of the function
      t += dt;
    }
  }


  MPI_Finalize();
	return 0;
}