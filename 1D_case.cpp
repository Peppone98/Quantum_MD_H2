
/**** Oscillatore armonico in 1 direzione ****/

/* -1/2(y_i+1 - 2y_i + y_i-1) + 1/2 x_i^2 y_i = ey_i */

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

using namespace std;

const double L = 8;
const int n = 50;
const double h = L/n;

void fill_pot(double pot[n]);
void fill_matrix(gsl_matrix *A, double pot[n]);
void create_eval_evec(gsl_matrix *A, gsl_matrix *B, gsl_vector *C);


int main(){
    gsl_matrix *A = gsl_matrix_alloc(n, n);
    gsl_matrix *evec = gsl_matrix_alloc(n, n);
    gsl_vector *eval = gsl_vector_alloc(n);
    double pot[n];
    fill_pot(pot);
    fill_matrix(A, pot);
    create_eval_evec(A, evec, eval);
    cout << "Eigenvalues: " << endl;
    cout << gsl_vector_get(eval, 0)/(h*h) << endl;
    cout << gsl_vector_get(eval, 1)/(h*h) << endl;
    cout << gsl_vector_get(eval, 2)/(h*h) << endl;

    ofstream my_eigenvec;
	my_eigenvec.open("eigenvectors.txt", ios :: out | ios :: trunc);
    for(int i=0; i<n; i++){
        my_eigenvec << gsl_matrix_get(evec, i, 0) << "   " << gsl_matrix_get(evec, i, 1) << "   " << gsl_matrix_get(evec, i, 2) << endl;
    }
    my_eigenvec.close();
}

void fill_matrix(gsl_matrix *A, double pot[n]){
    int i;
    gsl_matrix_set(A, 0, 0, 1.0 + pot[0]*h*h);
    gsl_matrix_set(A, 0, 1, -0.5);

    /* construct rows [1:n-2] */
    for (i = 1; i < n - 1; ++i){
        gsl_matrix_set(A, i, i + 1, -0.5);
        gsl_matrix_set(A, i, i, 1.0 + pot[i]*h*h);
        gsl_matrix_set(A, i, i - 1, -0.5);
    }

    /* construct last row */
    gsl_matrix_set(A, n - 1, n - 1, 1.0 + pot[n-1]*h*h);
    gsl_matrix_set(A, n - 1, n - 2, -0.5);
}


void fill_pot(double pot[n]){
    double x;
    int i;
    for(i=0; i<n; i++){
        x = -L/2. + i*h;
        pot[i] = 0.5*x*x;
    }
}


void create_eval_evec(gsl_matrix *A, gsl_matrix *B, gsl_vector *C){
    gsl_eigen_symmv_workspace *w = gsl_eigen_symmv_alloc(n); 
	gsl_eigen_symmv(A, C, B, w); 
	gsl_eigen_symmv_free(w); 
	gsl_eigen_symmv_sort(C, B, GSL_EIGEN_SORT_VAL_ASC);
}


void kronecker(gsl_matrix *A, gsl_matrix *B, gsl_matrix *C, int dim_A, int dim_B){
    int i, j, k, l;
    double a_ij, b_kl;
    for(i=0; i<dim_A; i++){
        for(j=0; j<dim_A; j++){
            a_ij = gsl_matrix_get(A, i, j);

            for(k=0; k<dim_B; k++){
                for(l=0; l<dim_B; l++){
                    b_kl = gsl_matrix_get(B, k, l);
                    gsl_matrix_set(C, i*dim_B + k, j*dim_B + l, a_ij*b_kl);
                }
            }
        }
    }
}