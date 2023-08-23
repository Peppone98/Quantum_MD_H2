
/********** SHAKE ALGORITHM IMPLEMENTATION ********/

#include "definitions.h"

using namespace std;


double dsigma_dX(gsl_vector *c, R R_A, R R_B){
    gsl_matrix *dS_dX = gsl_matrix_alloc(2*N, 2*N);
    gsl_vector *tmp = gsl_vector_alloc(2*N);
    create_dS_dX(dS_dX, R_A, R_B);
    double result = 0.0;
    
    /**** c^T * dS_dX * c ****/
    gsl_blas_dgemv(CblasNoTrans, 1., dS_dX, c, 0., tmp);
    gsl_blas_ddot(c, tmp, &result);
    gsl_matrix_free(dS_dX);
    gsl_vector_free(tmp);
    return result;
}


/**** sigma is the constraint c^T * dS_dX * c - 1 = 0 ****/
double sigma(gsl_vector *c, R R_A, R R_B){
    gsl_matrix *S = gsl_matrix_alloc(2*N, 2*N);
    gsl_vector *tmp = gsl_vector_alloc(2*N);
    create_S(S, R_A, R_B);
    double result = 0.0;
    
    /**** c^T * dS_dX * c - 1 ****/
    gsl_blas_dgemv(CblasNoTrans, 1., S, c, 0., tmp);
    gsl_blas_ddot(c, tmp, &result);
    return result - 1.0;
}



double Get_X_shake(double X, double X_old, gsl_vector *c, R R_A, R R_B, double dE0_dX, double lambda_guess, double eps){
    double denominator, numerator = 1.0, X_new;

    /**** Compute the (partial) X_new with the first guess of lambda (provided by the CP equilibration run) ****/
    double dsig_dX = dsigma_dX(c, R_A, R_B);
    X_new = 2.0*X - X_old - h_N*h_N/M_N*(2.0*dE0_dX + lambda_guess*dsig_dX);
    
    /**** Shake iterative cycle ****/
    while(fabs(numerator) > 1E-3){
        /**** Update the R_B ****/
        R_B.x = X_new;

        /**** (complete) denominator: C^T * dS(X[n])/dX * C * C^T * dS(X_new)/dX * C ****/
        denominator = dsig_dX*dsigma_dX(c, R_A, R_B);

        /**** numerator: C^T * S(X_new) * C - 1 ****/
        numerator = sigma(c, R_A, R_B);

        /**** Update the X and R_B ****/
        X_new = X_new - h_N*h_N/M_N*dsig_dX*numerator/denominator;
    }
    return X_new;
}


