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


double Nuclear_kinetic_en(double X_new, double X_old){
	double v = (X_new - X_old)/h_N;

	/**** The reduced mass is 0.5*M_N, this justifies the 0.25 factor ****/
	return 0.25*M_N*v*v;
}


double Electron_electron_en(gsl_vector *c, gsl_matrix *F){
	double result = 0.0;
	gsl_vector *tmp = gsl_vector_alloc(2*N);

    /**** c^T * F * c ****/
    gsl_blas_dgemv(CblasNoTrans, 1., F, c, 0., tmp);
    gsl_blas_ddot(c, tmp, &result);
    return result;
}


double One_body(gsl_vector *c, gsl_matrix *H){
	double result = 0.0;
	gsl_vector *tmp = gsl_vector_alloc(2*N);

    /**** c^T * H * c ****/
    gsl_blas_dgemv(CblasNoTrans, 1., H, c, 0., tmp);
    gsl_blas_ddot(c, tmp, &result);
    return 2.0*result;
}


double XC_energy(gsl_vector *c, gsl_matrix *V_xc){
	double result = 0.0;
	gsl_vector *tmp = gsl_vector_alloc(2*N);

    /**** c^T * H * c ****/
    gsl_blas_dgemv(CblasNoTrans, 1., V_xc, c, 0., tmp);
    gsl_blas_ddot(c, tmp, &result);
    return 2.0*result;
}


double Fictitious_kin_energy(gsl_vector *c, gsl_vector *c_old, gsl_matrix *S, gsl_matrix *dS_dX, double X, double X_old, R R_A, R R_B){
    int p, q;
    double c_p, c_q, c_old_p, c_old_q, c_p_dot, c_q_dot, sum_1=0.0, sum_2=0.0, sum_3=0.0, X_dot;

    /**** Get the velocity ****/
    X_dot = (X - X_old)/h_N;

	/**** Cycle for the three terms of the braket product <psi_dot|psi_dot> ****/
    for(p=0; p<N; p++){
        c_p = gsl_vector_get(c, p);
        c_old_p = gsl_vector_get(c_old, p);
        for(q=0; q<N; q++){
            c_q = gsl_vector_get(c, q);
            c_old_q = gsl_vector_get(c_old, q);
            c_p_dot = (c_p - c_old_p)/h_N;
            c_q_dot = (c_q - c_old_q)/h_N;
            sum_1 += c_p_dot*c_q_dot*gsl_matrix_get(S, p, q);
			sum_2 += c_p_dot*c_q*gsl_matrix_get(dS_dX, p, q)*X_dot;
			sum_3 += c_p*c_q*laplacian(a[p], a[q], R_A, R_B)*X_dot*X_dot;
        }
    }

	/**** Factor 0.5*m because it's a kinetic energy ****/
	return 0.5*m*(sum_1 + sum_2 + sum_3);

}