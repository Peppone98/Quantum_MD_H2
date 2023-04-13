#ifndef DEFINITIONS_H
#define DEFINITIONS_H

#include <iostream>
#include <iomanip>
#include <fstream>
#include <cstdlib>
#include <string>
#include <cmath>
#include <gsl/gsl_math.h>
#include <gsl/gsl_eigen.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>

struct R {
    double x;
    double y;
    double z;
};

const double M_N = 2000.; /* nuclear mass */
const double gamma_N = 15.; /* nuclear damping */
const int CP_iter = 50; /* Car-Parrinello iterations */
const double m = 2.; /* fictitious mass */
const double gamma_el= 1.; /* electronic damping */
const double h = 0.1; /* electronic time scale */
const double h_N = 43*h; /* nuclear time scale*/
const int N = 4;
const double a[N] = {13.00773, 1.962079, 0.444529, 0.1219492};
const double pi = 3.1415926; 


/******** ROUTINES FOR MATRIX ELEMENTS *****/
double scalar_prod(R R_A, R R_B);
double K(double alpha, double beta, R R_A, R R_B);
R R_weighted(double alpha, double beta, R R_A, R R_B);
double overlap(double alpha, double beta, R R_A, R R_B);
double laplacian(double alpha, double beta, R R_A, R R_B);
double el_nucl(double alpha, double beta, R R_A, R R_B, R R_C);
double F0(double x);
double direct_term(double alpha, double beta, R R_A, R R_B, double alpha_prime, double beta_prime, R R_A_prime, R R_B_prime);


/******* ROUTINES FOR SOLVING ROOTHAN EQUATION ****/
void diag_S(gsl_matrix *S, gsl_matrix *U);
void create_eval_evec(gsl_matrix *A, gsl_matrix *B, gsl_vector *C);
void create_V(gsl_matrix *U, gsl_vector *L);
void multiply(gsl_matrix *A, gsl_matrix *B, gsl_matrix *C);
double solve_FC_eSC(gsl_matrix *F, gsl_matrix *V, gsl_matrix *U);


/******** ROUTINES FOR FILLING MATRICES ******/
void create_S(gsl_matrix *S, R R_A, R R_B);
void one_body_H(gsl_matrix *H, R R_A, R R_B);
void build_Q(double Q[2*N][2*N][2*N][2*N], R R_A, R R_B);
void two_body_F(double Q[2*N][2*N][2*N][2*N], gsl_vector *c, gsl_matrix *F, R R_A, R R_B);


/******** COMPUTE TOTAL ENERGY AND PRINT SOLUTION ****/
double compute_E0(double Q[2*N][2*N][2*N][2*N], gsl_vector *c, gsl_matrix *H, R R_A, R R_B);
void print_orbital(gsl_vector *c_new, R R_A, R R_B);


/******** EVOLVE THE COEFFICIENTS IN TIME ********/
double update_c(gsl_matrix *F, gsl_matrix *S, gsl_vector *c, gsl_vector *c_old);
void partial_evolution(gsl_matrix *F, gsl_vector *c, gsl_vector *c_old);
void complete_evolution(gsl_matrix *S, gsl_vector *c, gsl_vector *c_old, double lambda);
double solve_eq2degree(double a, double b, double c, int solution);
double lowest_positive_root(double a, double b, double c);
double normalization(gsl_vector *c, gsl_matrix *S);
double get_C(gsl_vector *c, gsl_matrix *S);
double get_B(gsl_vector *c, gsl_vector *c_old, gsl_matrix *S);
double get_A(gsl_vector *c_old, gsl_matrix *S);


/******** NUCLEAR MOTION MATRIX ELEMENTS ********/
double doverlap_dX(double alpha, double beta, R R_A, R R_B);
double dlaplacian_dX(double alpha, double beta, R R_A, R R_B);
double del_nucl_dX(double alpha, double beta, double X, R R_A, R R_B);
double dF0_dt(double t);
double ddirect_term_dX(double alpha, double beta, R R_A, R R_B, double alpha_prime, double beta_prime, R R_A_prime, R R_B_prime, double X);

/******** ROUTINES FOR FILLING NUCLEAR MATRICES ******/
void create_dS_dX(gsl_matrix *S, R R_A, R R_B);
void one_body_dH_dX(gsl_matrix *H, R R_A, R R_B, double X);
void build_dQ_dX(double Q[2*N][2*N][2*N][2*N], R R_A, R R_B, double X);

/********** EVOLVE X IN TIME ********/
double evolve_X(double Q[2*N][2*N][2*N][2*N], gsl_vector *c, gsl_matrix *H, gsl_matrix *S, R *R_B, double lambda, double dE0_dX, double X, double X_old);
double compute_dE0_dX(double Q[2*N][2*N][2*N][2*N], gsl_vector *c, gsl_matrix *H, R R_A, R R_B);


#endif