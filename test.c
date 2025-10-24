#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

// Hamiltonian RHS for one trajectory
void hamiltonian_rhs(double* q, double* p, double* dq, double* dp, int d) {
    for (int i = 0; i < d; i++) {
        int left  = (i == 0) ? d-1 : i-1;
        int right = (i == d-1) ? 0 : i+1;

        dq[i] = p[i];
        dp[i] = -cos(q[i] + (q[left] + q[right])/2) * 2.0; // chain rule
    }
}

// RK4 step for one trajectory
void rk4_step(double* q, double* p, int d, double dt,
              double* k1_q, double* k1_p,
              double* k2_q, double* k2_p,
              double* k3_q, double* k3_p,
              double* k4_q, double* k4_p,
              double* temp_q, double* temp_p) {

    hamiltonian_rhs(q, p, k1_q, k1_p, d);

    for (int i = 0; i < d; i++) { temp_q[i] = q[i] + dt*k1_q[i]/2; temp_p[i] = p[i] + dt*k1_p[i]/2; }
    hamiltonian_rhs(temp_q, temp_p, k2_q, k2_p, d);

    for (int i = 0; i < d; i++) { temp_q[i] = q[i] + dt*k2_q[i]/2; temp_p[i] = p[i] + dt*k2_p[i]/2; }
    hamiltonian_rhs(temp_q, temp_p, k3_q, k3_p, d);

    for (int i = 0; i < d; i++) { temp_q[i] = q[i] + dt*k3_q[i]; temp_p[i] = p[i] + dt*k3_p[i]; }
    hamiltonian_rhs(temp_q, temp_p, k4_q, k4_p, d);

    for (int i = 0; i < d; i++) {
        q[i] += dt*(k1_q[i] + 2*k2_q[i] + 2*k3_q[i] + k4_q[i])/6.0;
        p[i] += dt*(k1_p[i] + 2*k2_p[i] + 2*k3_p[i] + k4_p[i])/6.0;
    }
}

int main() {
    int d = 200;
    int n_traj = 10; // adjust based on RAM
    int t_steps = 1000;
    double dt = 0.01;

    // Allocate main arrays
    double* q = malloc(n_traj * d * sizeof(double));
    double* p = malloc(n_traj * d * sizeof(double));

    // Allocate temporary arrays for RK4 (one trajectory at a time)
    double *k1_q = malloc(d*sizeof(double)), *k1_p = malloc(d*sizeof(double));
    double *k2_q = malloc(d*sizeof(double)), *k2_p = malloc(d*sizeof(double));
    double *k3_q = malloc(d*sizeof(double)), *k3_p = malloc(d*sizeof(double));
    double *k4_q = malloc(d*sizeof(double)), *k4_p = malloc(d*sizeof(double));
    double *temp_q = malloc(d*sizeof(double)), *temp_p = malloc(d*sizeof(double));

    // Initialize random initial conditions
    srand(time(NULL));
    for (int n = 0; n < n_traj; n++)
        for (int i = 0; i < d; i++) {
            q[n*d + i] = 0.1 * ((double)rand()/RAND_MAX - 0.5);
            p[n*d + i] = 0.1 * ((double)rand()/RAND_MAX - 0.5);
        }

    double observable = 0.0;

    clock_t start = clock();

    for (int step = 0; step < t_steps; step++) {
        for (int n = 0; n < n_traj; n++) {
            rk4_step(q + n*d, p + n*d, d, dt,
                     k1_q, k1_p, k2_q, k2_p, k3_q, k3_p, k4_q, k4_p, temp_q, temp_p);

            // // Example observable: mean q for this trajectory
            // double mean_q = 0.0;
            // for (int i = 0; i < d; i++) mean_q += q[n*d + i];
            // mean_q /= d;

            // observable += mean_q;
        }

        // if (step % 100 == 0)
        //     printf("Step %d done, current observable: %.5f\n", step, observable / ((step+1)*n_traj));
    }

    clock_t end = clock();
    double elapsed = (double)(end - start) / CLOCKS_PER_SEC;
    printf("Time elapsed for evolution: %.3f seconds\n", elapsed);

    observable /= (t_steps * n_traj);
    printf("Final observable: %.5f\n", observable);

    free(q); free(p);
    free(k1_q); free(k1_p); free(k2_q); free(k2_p);
    free(k3_q); free(k3_p); free(k4_q); free(k4_p);
    free(temp_q); free(temp_p);

    return 0;
}