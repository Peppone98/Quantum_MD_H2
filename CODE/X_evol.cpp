/********** EVOLVE X IN TIME ********/

#include "definitions.h"

using namespace std;

/**** First step: X(t) = X(t - h_N), i.e., X(h_N) = X(0) ****/
double evolve_X(double Q[2*N][2*N][2*N][2*N], gsl_vector *c, gsl_matrix *H, gsl_matrix *S, R *R_B, double lambda, double dE0_dX, double X, double X_old){
    gsl_vector *c_tmp = gsl_vector_alloc(2*N);
    gsl_vector_memcpy(c_tmp, c);

    /**** "normalization" also changes the c, that is the reason of c_tmp ****/
    double norm = normalization(c_tmp, S);
    double X_new = X*M_N - 0.5*X_old*(M_N - gamma_N*h_N);
    X_new -= h_N*h_N*(dE0_dX - lambda*norm);
    X_new = 2.*X_new/(M_N + gamma_N*h_N);

    /**** Update the nuclear position B ****/
    R_B->x = X_new;  
    return X_new;
}

/**** The matrices must be dQ/dX and dH/dX when called ****/
double compute_dE0_dX(double Q[2*N][2*N][2*N][2*N], gsl_vector *c, gsl_matrix *H, R R_A, R R_B){
    double X = scalar_prod(R_A, R_B);
    double dE0_dX=0.;
    dE0_dX = compute_E0(Q, c, H, R_A, R_B) - 1/X;
    dE0_dX -= 1/(X*X);
    return dE0_dX;
}


