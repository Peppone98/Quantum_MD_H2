
/********** SHAKE ALGORITHM IMPLEMENTATION ********/

#include "definitions.h"

using namespace std;

void partial_update_shake(double X, double X_old, R *R_B, double dE0_dX, double lambda_guess, double dsig_dX){
    
    /**** Update of the X with the guess of the Lagrange multiplier lambda ****/
    double X_new = 2.0*X - X_old - h_N*h_N/M_N*(2.0*dE0_dX + lambda_guess*dsig_dX);
    R_B->x = X_new;
}


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



double sigma(gsl_vector *c, R R_A, R R_B){
    gsl_matrix *S = gsl_matrix_alloc(2*N, 2*N);
    gsl_vector *tmp = gsl_vector_alloc(2*N);
    create_S(S, R_A, R_B);
    double result = 0.0;
    
    /**** c^T * dS_dX * c ****/
    gsl_blas_dgemv(CblasNoTrans, 1., S, c, 0., tmp);
    gsl_blas_ddot(c, tmp, &result);
    return result;
}



double Get_X_shake(double X, double X_old, gsl_vector *c, R R_A, R R_B, double dE0_dX, double lambda_guess, double eps){
    double denominator, numerator, new_lambda = 0.0, old_lambda = lambda_guess;

    /**** This is needed to complete the evolution ****/
    double dsig_dX = dsigma_dX(c, R_A, R_B);

    /**** First iteration ****/
    denominator = dsigma_dX(c, R_A, R_B);
    std::cout << "Before: " << R_B.x << endl;
    partial_update_shake(X, X_old, &R_B, dE0_dX, lambda_guess, denominator);
    std::cout << "After: " << R_B.x << endl;
    denominator = denominator*dsigma_dX(c, R_A, R_B);
    numerator = sigma(c, R_A, R_B) - 1.0;
    new_lambda = lambda_guess - numerator/denominator;

    /**** Iterative shake cycle ****/
    while(fabs(new_lambda - old_lambda) > eps){
        old_lambda = new_lambda;
        denominator = dsigma_dX(c, R_A, R_B);
        std::cout << "Before: " << R_B.x << endl;
        partial_update_shake(X, X_old, &R_B, dE0_dX, new_lambda, denominator);
        std::cout << "After: " << R_B.x << endl;
        denominator = denominator*dsigma_dX(c, R_A, R_B);
        numerator = sigma(c, R_A, R_B) - 1.0;
        new_lambda = old_lambda - numerator/denominator;
    }
    return 2.0*X - X_old - h_N*h_N/M_N*(2.0*dE0_dX - new_lambda*dsig_dX);
}


