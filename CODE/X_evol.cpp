/********** EVOLVE X IN TIME ********/

#include "definitions.h"

using namespace std;

/**** First step: X(t) = X(t - h_N), i.e., X(h_N) = X(0) ****/
double evolve_X(gsl_vector *c, gsl_matrix *S, R *R_B, double lambda, double dE0_dX, double X, double X_old){
    gsl_vector *c_tmp = gsl_vector_alloc(2*N);
    gsl_vector_memcpy(c_tmp, c);

    /**** The function "normalization" also changes the c, that is the reason of c_tmp ****/
    double norm = normalization(c_tmp, S);

    /**** This is a mistery ****/
    double lambda_tmp = lambda*0.25*(m + h*gamma_el);

    /**** Equation of motion for X ****/
    double X_new = 2.*X*M_N - X_old*(M_N - 0.5*gamma_N*h_N);

    /* mysterious factor 2 */
    X_new -= 2.*h_N*h_N*(dE0_dX + lambda_tmp*norm);
    X_new = X_new/(M_N + 0.5*gamma_N*h_N);
    
    /**** Update the nuclear position B ****/
    R_B->x = X_new;  
    return X_new;
}

/**** The matrices must be dQ/dX and dH/dX when called ****/
double compute_dE0_dX(gsl_matrix *F, gsl_matrix *H, gsl_vector *c, double X){
    double dE0_dX=0.;
    dE0_dX = compute_E0(F, H, c) - 1./(X*X);
    return dE0_dX;
}


