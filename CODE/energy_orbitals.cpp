/******** COMPUTE TOTAL ENERGY AND PRINT SOLUTION ****/

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

double compute_E0(double Q[2*N][2*N][2*N][2*N], gsl_vector *c, gsl_matrix *H, R R_A, R R_B){
    int p, q, t, s;
    double E0 = 0., c_p=0., c_q=0., c_t=0., c_s=0., H_pq=0.;
    for(p=0; p<2*N; p++){
        for(q=0; q<2*N; q++){
            c_p = gsl_vector_get(c, p);
            c_q = gsl_vector_get(c, q);
            H_pq = gsl_matrix_get(H, p, q);
            E0 += 2.*c_p*c_q*H_pq;
            for(t=0; t<2*N; t++){
                for(s=0; s<2*N; s++){
                    c_t = gsl_vector_get(c, t);
                    c_s = gsl_vector_get(c, s);
                    E0 += c_p*c_q*c_t*c_s*Q[p][t][q][s];
                }
            }
        }
    }
    return E0 + 1/sqrt(scalar_prod(R_A, R_B));
}