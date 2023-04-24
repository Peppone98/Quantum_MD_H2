/******** COMPUTE TOTAL ENERGY AND PRINT SOLUTION *******/

#include "definitions.h"

using namespace std;

void print_orbital(gsl_vector *c_new, R R_A, R R_B){
	double *f, h=0.003;
	R r;
	int mesh = 1000, i, n; 
	f = new double[mesh];
	ofstream myfile;
	myfile.open("Psi_0.txt", ios :: out | ios :: trunc);
    myfile.precision(6);
	for(i=0; i<mesh; i++){
		r.x = -1. + R_A.x + i*h;
		for(n=0; n<N; n++){
			f[i] += gsl_vector_get(c_new, n)*exp(-(r.x - R_A.x)*(r.x - R_A.x)*a[n]);
			f[i] += gsl_vector_get(c_new, n + N)*exp(-(r.x - R_B.x)*(r.x - R_B.x)*a[n]);
		}
		myfile << r.x << "       " << f[i] << endl; 	
	}
	myfile.close();
}



double compute_E0(gsl_matrix *F, gsl_matrix *H, gsl_vector *c){
	/**** Create copies of F and c to avoid unwanted changes ****/
    gsl_vector *tmp_c = gsl_vector_alloc(2*N);
	gsl_matrix *tmp_F = gsl_matrix_alloc(2*N, 2*N);
	gsl_matrix_memcpy(tmp_F, F);

    /**** Add H to the Fock matrix ****/
    gsl_matrix_add(tmp_F, H);
    gsl_blas_dgemv(CblasNoTrans, 1., tmp_F, c, 0., tmp_c);

    double result = 0.;
    for(int i=0; i<2*N; i++){
        result += gsl_vector_get(c, i)*gsl_vector_get(tmp_c, i);
    }
    return result;
}