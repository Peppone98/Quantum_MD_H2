/******** ROUTINES FOR FILLING MATRICES ******/

#include "definitions.h"

void create_S(gsl_matrix *S, R R_A, R R_B){
	int p, q;
	double val1=0., val2=0.;
	for(p=0; p<N; p++){
		for(q=0; q<=p; q++){

			val1 = overlap(a[p], a[q], R_A, R_A);
	        gsl_matrix_set(S, p, q, val1);
            gsl_matrix_set(S, q, p, val1);
            gsl_matrix_set(S, p + N, q + N, val1);
            gsl_matrix_set(S, q + N, p + N, val1);

            val2 = overlap(a[p], a[q], R_A, R_B);
	        gsl_matrix_set(S, p, q + N, val2);
            gsl_matrix_set(S, q, p + N, val2);
            gsl_matrix_set(S, p + N, q, val2);
            gsl_matrix_set(S, q + N, p, val2);

		}
	}
}


void one_body_H(gsl_matrix *H, R R_A, R R_B){
    int p, q;
    double val1=0., val2=0.;
    for(p = 0; p < N; p++){
        for(q = 0; q <= p; q++){

            val1 = laplacian(a[p], a[q], R_A, R_A);
            val1 += el_nucl(a[p], a[q], R_A, R_A, R_A);
            val1 += el_nucl(a[p], a[q], R_A, R_A, R_B);
            gsl_matrix_set(H, p, q, val1);  
            gsl_matrix_set(H, q, p, val1); 
            gsl_matrix_set(H, p + N, q + N, val1);  
            gsl_matrix_set(H, q + N, p + N, val1); 

            val2 = laplacian(a[p], a[q], R_A, R_B);
            val2 += el_nucl(a[p], a[q], R_A, R_B, R_A);
            val2 += el_nucl(a[p], a[q], R_A, R_B, R_B);
            gsl_matrix_set(H, p, q + N, val2);  
            gsl_matrix_set(H, q, p + N, val2); 
            gsl_matrix_set(H, p + N, q, val2);  
            gsl_matrix_set(H, q + N, p, val2);
            
            }
        }
}

void build_Q(double Q[2*N][2*N][2*N][2*N], R R_A, R R_B){
    int p, q, t, s;
    for(p=0; p<N; p++){
        for(q=0; q<N; q++){
            for(t=0; t<N; t++){
                for(s=0; s<N; s++){
                    Q[p][q][t][s] = direct_term(a[p], a[q], R_A, R_A, a[t], a[s], R_A, R_A);
                    Q[p][q][t][s + N] = direct_term(a[p], a[q], R_A, R_A, a[t], a[s], R_A, R_B);
                    Q[p][q][t + N][s] = direct_term(a[p], a[q], R_A, R_A, a[t], a[s], R_B, R_A);
                    Q[p][q][t + N][s + N] = direct_term(a[p], a[q], R_A, R_A, a[t], a[s], R_B, R_B);
                    Q[p][q + N][t][s] = direct_term(a[p], a[q], R_A, R_B, a[t], a[s], R_A, R_A);
                    Q[p][q + N][t][s + N] = direct_term(a[p], a[q], R_A, R_B, a[t], a[s], R_A, R_B);
                    Q[p][q + N][t + N][s] = direct_term(a[p], a[q], R_A, R_B, a[t], a[s], R_B, R_A);
                    Q[p][q + N][t + N][s + N] = direct_term(a[p], a[q], R_A, R_B, a[t], a[s], R_B, R_B);
                    Q[p + N][q][t][s] = direct_term(a[p], a[q], R_B, R_A, a[t], a[s], R_A, R_A);
                    Q[p + N][q + N][t][s] = direct_term(a[p], a[q], R_B, R_B, a[t], a[s], R_A, R_A);
                    Q[p + N][q][t + N][s] = direct_term(a[p], a[q], R_B, R_A, a[t], a[s], R_B, R_A);
                    Q[p + N][q][t][s + N] = direct_term(a[p], a[q], R_B, R_A, a[t], a[s], R_A, R_B);
                    Q[p + N][q + N][t + N][s] = direct_term(a[p], a[q], R_B, R_B, a[t], a[s], R_B, R_A);
                    Q[p + N][q + N][t][s + N] = direct_term(a[p], a[q], R_B, R_B, a[t], a[s], R_A, R_B);
                    Q[p + N][q][t + N][s + N] = direct_term(a[p], a[q], R_B, R_A, a[t], a[s], R_B, R_B);
                    Q[p + N][q + N][t + N][s + N] = direct_term(a[p], a[q], R_B, R_B, a[t], a[s], R_B, R_B);
                }
            }
        }
    }
}

void two_body_F(double Q[2*N][2*N][2*N][2*N], gsl_vector *c, gsl_matrix *U, gsl_matrix *F, R R_A, R R_B){
    int p, q, t, s;
    double val = 0., c1, c2;
    for(p=0; p<2*N; p++){
        for(q=0; q<2*N; q++){
            val = 0.;
            for(t=0; t<2*N; t++){
                for(s=0; s<2*N; s++){
                    c1 = gsl_vector_get(c, t);
                    c2 = gsl_vector_get(c, s);
                    val += c1*c2*Q[p][t][q][s];
                    gsl_matrix_set(F, p, q, val);
                }
            }
        }
    }
}
