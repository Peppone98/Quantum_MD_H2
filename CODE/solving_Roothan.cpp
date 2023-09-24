/******* ROUTINES FOR SOLVING ROOTHAN EQUATION ****/

#include "definitions.h"

using namespace std;

void diag_S(gsl_matrix *S, gsl_matrix *U){
	gsl_vector *L = gsl_vector_alloc(2*N);
	create_eval_evec(S, U, L); 
	create_V(U, L);
	gsl_vector_free(L);
}

void create_V(gsl_matrix *U, gsl_vector *L){
	int i;
	double x=0.;
	for(i=0; i<2*N; i++){
		x = gsl_vector_get(L, i);
		gsl_vector_set(L, i, 1./sqrt(x)); 
	}
	gsl_matrix_scale_columns(U, L); 
}

/*The diagonal and lower triangular part of A are 
destroyed during the computation, but the strict upper 
triangular part is not referenced*/
void create_eval_evec(gsl_matrix *A, gsl_matrix *B, gsl_vector *C){
    gsl_eigen_symmv_workspace *w = gsl_eigen_symmv_alloc(2*N); 
	gsl_eigen_symmv(A, C, B, w); 
	gsl_eigen_symmv_free(w); 
	gsl_eigen_symmv_sort(C, B, GSL_EIGEN_SORT_VAL_ASC);
}
/**** eigenvectors are normalised to unit magnitude ****/


void multiply(gsl_matrix *A, gsl_matrix *B, gsl_matrix *C){
	int i, j, k;
	double val, f1, f2; 
	for(i=0; i<2*N; i++){
		for(j=0; j<2*N; j++){
			for(k=0; k<2*N; k++){
				f1 = gsl_matrix_get(A, i, k);
				f2 = gsl_matrix_get(B, k, j);
				val += f1*f2;
			}				
			gsl_matrix_set(C, i, j, val);
			val = 0.;
		}
	}
}


double solve_FC_eSC(gsl_matrix *F, gsl_matrix *V, gsl_matrix *U){
	double val=0.; 
	gsl_matrix_memcpy(U, V);
	gsl_matrix *tmp = gsl_matrix_alloc(2*N, 2*N);
	gsl_vector *L = gsl_vector_alloc(2*N);

	/**** Compute V^T * F * V, where F is the Fock matrix ****/
	multiply(F, U, tmp); 
	gsl_matrix_transpose(U); 
	multiply(U, tmp, F); 

	/**** Diagonalize Fock matrix ****/
	create_eval_evec(F, tmp, L); 
	gsl_matrix_transpose(U);

	/**** Now, tmp contains the eigenvectors of F (see above lines), while U contains V ****/ 
	multiply(U, tmp, U); 
	val = gsl_vector_get(L, 0);
	gsl_matrix_free(tmp);
	gsl_vector_free(L);
	return val;
}
