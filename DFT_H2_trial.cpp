
/**** Costruisco una mesh, calcolo la matrice, aggiungo la diagonale ****/


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

const int len = 1000;
const double L = 5;
const double h = L/len;

void kronecker(gsl_matrix *A, gsl_matrix *B, gsl_matrix *C, int dim_A, int dim_B);


int main(){
    int dim_A = 2, dim_B = 5, i, j;
    gsl_matrix *A = gsl_matrix_alloc(dim_A, dim_A);
    gsl_matrix *B = gsl_matrix_alloc(dim_B, dim_B);
    cout << "Matrix A: " << endl;
    for(i=0; i<dim_A; i++){
        for(j=0; j<dim_A; j++){
            gsl_matrix_set(A, i, j, i + j);
            cout << gsl_matrix_get(A, i, j) << "     ";
        }
        cout << "     " << endl;
    }
    
    cout << "Matrix B: " << endl;
    for(i=0; i<dim_B; i++){
        for(j=0; j<dim_B; j++){
            gsl_matrix_set(B, i, j, i - j);
            cout << gsl_matrix_get(B, i, j) << "     ";
        }
        cout << "     " << endl;
    }


    double x[len], y[len], z[len];
    for(i=0; i<len; i++){
        x[i] = -L/2. + i*h;
        y[i] = -L/2. + i*h;
        z[i] = -L/2. + i*h;
    }



    double dim_C = dim_A*dim_B;
    gsl_matrix *C = gsl_matrix_alloc(dim_C, dim_C);
    kronecker(A, B, C, dim_A, dim_B);

    cout << "Matrix C: " << endl;
    for(i=0; i<dim_C; i++){
        for(j=0; j<dim_C; j++){
            cout << gsl_matrix_get(C, i, j) << "     ";
        }
        cout << "   " << endl;
    }
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

/* non so esattamente come mettere il potenziale nelle matrici */

void harm_pot(double pot[len][len][len]){
    int i, j, k;
    double pot=0.;
    for(i=0; i<len; i++){
        for(j=0; j<len; j++){
            for(k=0; k<len; k++){
                pot[i][j][k] = 0.5*(x[i]*x[i] + y[j]*y[j] + z[k]*z[k]);
            }
        }
    }
}