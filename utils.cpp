// File for auxiliary functions (ex. Omega)
#include "ode_func.h"
#include "utils.h"

double omega_0(double *Y, double T){
    int     j,k;   
    k = 2*dim;
    double  omega=0.0;
    
    for(j =0; j<dim; j++){
        omega+=Y[k*(N+2)+1]*Y[k*N+dim+j]*Y[k*N+dim+j]*(1.0/Tr-1.0/T)/m +
               Y[k*(N+2)]*Y[k+dim+j]*Y[k+dim+j]*(1.0/Tl-1.0/T)/m;
    } 

    return omega;
}
