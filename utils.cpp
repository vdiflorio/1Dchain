// File for auxiliary functions (ex. Omega)
#include "ode_func.h"

double omega_0(double *Y, double T){
    // Function to compute the dissipation function wrt the initial density
    int     j,k;   
    k = 2*dim;
    double  omega=0.0;
    
    for(j =0; j<dim; j++){
        omega+=Y[k*(N+2)+1]*Y[k*N+dim+j]*Y[k*N+dim+j]/m*(1/Tr-1/T)+Y[k*(N+2)]*Y[k+dim+j]*Y[k+dim+j]/m*(1/Tl-1/T);


    } 

    return omega;
}


double TTCF(double (*obs) (double *), double (*omega)( double *, double )){

}