/********** EVOLVE X IN TIME ********/

#include "definitions.h"

using namespace std;


double evolve_X(gsl_vector *c, gsl_matrix *S, R *R_B, double lambda, double dE0_dX, double X, double X_old, string s){
    gsl_vector *c_tmp = gsl_vector_alloc(2*N);
    gsl_vector_memcpy(c_tmp, c);
    double norm, lambda_tmp;

    /**** The function "normalization" also changes the c, that is the reason of c_tmp ****/
    norm = normalization(c_tmp, S);
    lambda_tmp = lambda;

    if(s == "CP"){
        /**** Readjust the lambda with the numerical factor (see "Evolution of C" algorithm in the text) ****/
        lambda_tmp = -lambda*0.25*(m + h*gamma_el);
    }

    if(s == "BO_CG"){
        /**** Eigenvalue adjustement ****/
        lambda_tmp = 2.0*lambda;
    }

    /**** Equation of motion for X. The 0.5 factor is due to the derivative (X_new - X_old)/2h_N in damped term ****/
    double X_new = 2.*X*M_N - X_old*(M_N - 0.5*gamma_N*h_N);

    /**** Factor 2: remember that the reduced mass is M_N/2 ****/
    X_new -= 2.*h_N*h_N*(dE0_dX - lambda_tmp*norm);
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

