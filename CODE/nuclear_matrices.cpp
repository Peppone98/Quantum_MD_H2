
/******** ROUTINES FOR FILLING NUCLEAR MATRICES ******/

#include "definitions.h"

using namespace std;

void create_dS_dX(gsl_matrix *S, R R_A, R R_B){
	int p, q;
	double val1=0., val2=0.;
	for(p=0; p<N; p++){
		for(q=0; q<=p; q++){

			val1 = doverlap_dX(a[p], a[q], R_A, R_A);
	        gsl_matrix_set(S, p, q, val1);
            gsl_matrix_set(S, q, p, val1);
            gsl_matrix_set(S, p + N, q + N, val1);
            gsl_matrix_set(S, q + N, p + N, val1);

            val2 = doverlap_dX(a[p], a[q], R_A, R_B);
	        gsl_matrix_set(S, p, q + N, val2);
            gsl_matrix_set(S, q, p + N, val2);
            gsl_matrix_set(S, p + N, q, val2);
            gsl_matrix_set(S, q + N, p, val2);

		}
	}
}



void one_body_dH_dX(gsl_matrix *H, R R_A, R R_B, double X){
    int p, q;
    double val1=0., val2=0.;
    if(sqrt(scalar_prod(R_A, R_B)) != X){
        cout << "Warning: give correct distances" << endl;
    }
    cout << "Quick check: "<< sqrt(scalar_prod(R_A, R_B)) << "  " << X << endl;
    for(p = 0; p < N; p++){
        for(q = 0; q <= p; q++){

            val1 = dlaplacian_dX(a[p], a[q], R_A, R_A);
            val1 -= del_nucl_dX(a[p], a[q], X, R_A, R_A);
            gsl_matrix_set(H, p, q, val1);  
            gsl_matrix_set(H, q, p, val1); 
            gsl_matrix_set(H, p + N, q + N, val1);  
            gsl_matrix_set(H, q + N, p + N, val1); 

            val2 = dlaplacian_dX(a[p], a[q], R_A, R_B);
            val2 -= del_nucl_dX(a[p], a[q], X, R_A, R_B);
            gsl_matrix_set(H, p, q + N, val2);  
            gsl_matrix_set(H, q, p + N, val2); 
            gsl_matrix_set(H, p + N, q, val2);  
            gsl_matrix_set(H, q + N, p, val2);
            
            }
        }
}



void build_dQ_dX(double Q[2*N][2*N][2*N][2*N], R R_A, R R_B, double X){
    int p, q, t, s;
    for(p=0; p<N; p++){
        for(q=0; q<N; q++){
            for(t=0; t<N; t++){
                for(s=0; s<N; s++){
                    Q[p][q][t][s] = ddirect_term_dX(a[p], a[q], R_A, R_A, a[t], a[s], R_A, R_A, X);
                    Q[p][q][t][s + N] = ddirect_term_dX(a[p], a[q], R_A, R_A, a[t], a[s], R_A, R_B, X);
                    Q[p][q][t + N][s] = ddirect_term_dX(a[p], a[q], R_A, R_A, a[t], a[s], R_B, R_A, X);
                    Q[p][q][t + N][s + N] = ddirect_term_dX(a[p], a[q], R_A, R_A, a[t], a[s], R_B, R_B, X);
                    Q[p][q + N][t][s] = ddirect_term_dX(a[p], a[q], R_A, R_B, a[t], a[s], R_A, R_A, X);
                    Q[p][q + N][t][s + N] = ddirect_term_dX(a[p], a[q], R_A, R_B, a[t], a[s], R_A, R_B, X);
                    Q[p][q + N][t + N][s] = ddirect_term_dX(a[p], a[q], R_A, R_B, a[t], a[s], R_B, R_A, X);
                    Q[p][q + N][t + N][s + N] = ddirect_term_dX(a[p], a[q], R_A, R_B, a[t], a[s], R_B, R_B, X);
                    Q[p + N][q][t][s] = ddirect_term_dX(a[p], a[q], R_B, R_A, a[t], a[s], R_A, R_A, X);
                    Q[p + N][q + N][t][s] = ddirect_term_dX(a[p], a[q], R_B, R_B, a[t], a[s], R_A, R_A, X);
                    Q[p + N][q][t + N][s] = ddirect_term_dX(a[p], a[q], R_B, R_A, a[t], a[s], R_B, R_A, X);
                    Q[p + N][q][t][s + N] = ddirect_term_dX(a[p], a[q], R_B, R_A, a[t], a[s], R_A, R_B, X);
                    Q[p + N][q + N][t + N][s] = ddirect_term_dX(a[p], a[q], R_B, R_B, a[t], a[s], R_B, R_A, X);
                    Q[p + N][q + N][t][s + N] = ddirect_term_dX(a[p], a[q], R_B, R_B, a[t], a[s], R_A, R_B, X);
                    Q[p + N][q][t + N][s + N] = ddirect_term_dX(a[p], a[q], R_B, R_A, a[t], a[s], R_B, R_B, X);
                    Q[p + N][q + N][t + N][s + N] = ddirect_term_dX(a[p], a[q], R_B, R_B, a[t], a[s], R_B, R_B, X);
                }
            }
        }
    }
}
