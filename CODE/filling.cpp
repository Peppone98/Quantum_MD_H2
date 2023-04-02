

#include <iostream>
#include <iomanip>
#include <fstream>
#include <cmath>
#include <cstdlib>
#include <string>
#include <gsl/gsl_math.h>
#include <gsl/gsl_eigen.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>

using namespace std;

struct R {
    double x;
    double y;
    double z;
};

const int N = 4;
const double pi = 3.1415926;
const double a[N] = {13.00773, 1.962079, 0.444529, 0.1219492};

void create_S(gsl_matrix *S, R R_A, R R_B);
void one_body_H(gsl_matrix *H, R R_A, R R_B);
void two_body_F(gsl_vector *c, gsl_matrix *F, R R_A, R R_B);

double scalar_prod(R R_A, R R_B);
double K(double alpha, double beta, R R_A, R R_B);
R R_weighted(double alpha, double beta, R R_A, R R_B);
double overlap(double alpha, double beta, R R_A, R R_B);
double laplacian(double alpha, double beta, R R_A, R R_B);
double el_nucl(double alpha, double beta, R R_A, R R_B, R R_C);
double F0(double x);
double direct_term(double alpha, double beta, R R_A, R R_B, double alpha_prime, double beta_prime, R R_A_prime, R R_B_prime);



int main(){
    R R_A, R_B;
    R_A.x = 0., R_A.y = 0., R_A.z = 0.;
    R_B.x = 1., R_B.y = 0., R_B.z = 0.;
    gsl_matrix *S = gsl_matrix_alloc(2*N, 2*N);
    gsl_matrix *H = gsl_matrix_alloc(2*N, 2*N);
    gsl_matrix *F = gsl_matrix_alloc(2*N, 2*N);
    gsl_vector *c = gsl_vector_alloc(2*N);
	gsl_vector_set_all(c, 1.);

    create_S(S, R_A, R_B);
    one_body_H(H, R_A, R_B);
    two_body_F(c, F, R_A, R_B);
    gsl_matrix_add(F, H);

    for(int i=0; i<2*N; i++){
        for(int j=0; j<2*N; j++){
            cout << gsl_matrix_get(H, i, j) << "    ";
        }
        cout << endl;
    }
    cout << direct_term(a[0], a[0], R_A, R_A, a[0], a[0], R_A, R_A) << endl;
    cout << direct_term(a[1], a[2], R_A, R_A, a[1], a[1], R_B, R_B) << endl;

}

/************ FUNCTIONS TO CREATE MATRIX ELEMETS ***************/

double scalar_prod(R R_A, R R_B){
    double x = R_A.x - R_B.x;
    double y = R_A.y - R_B.y;
    double z = R_A.z - R_B.z;
    return x*x + y*y + z*z;
}

double K(double alpha, double beta, R R_A, R R_B){
    return exp(-alpha*beta*scalar_prod(R_A, R_B)/(alpha + beta));
}

R R_weighted(double alpha, double beta, R R_A, R R_B){
    R R_w;
    R_A.x = alpha*R_A.x, R_A.y = alpha*R_A.y, R_A.z = alpha*R_A.z;
    R_B.x = beta*R_B.x, R_B.y = beta*R_B.y, R_B.z = beta*R_B.z;
    R_w.x = (R_A.x + R_B.x)/(alpha + beta);
    R_w.y = (R_A.y + R_B.y)/(alpha + beta);
    R_w.z = (R_A.z + R_B.z)/(alpha + beta);
    return R_w;
}

double overlap(double alpha, double beta, R R_A, R R_B){
	return pow(pi/(alpha + beta), 1.5)*K(alpha, beta, R_A, R_B);
}


double laplacian(double alpha, double beta, R R_A, R R_B){
    double tmp = 3. + 2.*log(K(alpha, beta, R_A, R_B));
    return alpha*beta/(alpha + beta)*tmp*overlap(alpha, beta, R_A, R_B);
}

double F0(double x){
    if(x == 0.){
        return 1;
    }else{
        return sqrt(pi/(4.*x))*erf(sqrt(x));
    }
}

double el_nucl(double alpha, double beta, R R_A, R R_B, R R_C){
    double tmp = -2*pi/(alpha + beta)*K(alpha, beta, R_A, R_B);
    R R_P = R_weighted(alpha, beta, R_A, R_B);
    return tmp*F0((alpha + beta)*scalar_prod(R_P, R_C));
}


double direct_term(double alpha, double beta, R R_A, R R_B, double alpha_prime, double beta_prime, R R_A_prime, R R_B_prime){
    double tmp = 2.*pow(pi, 2.5)/((alpha + alpha_prime)*(beta + beta_prime)*sqrt(alpha + alpha_prime + beta + beta_prime));
    tmp = tmp*K(alpha, alpha_prime, R_A, R_A_prime)*K(beta, beta_prime, R_B, R_B_prime);
    R R_P, R_Q;
    R_P = R_weighted(alpha, alpha_prime, R_A, R_A_prime);
    R_Q = R_weighted(beta, beta_prime, R_B, R_B_prime);
    return tmp*F0(-log(K(alpha + alpha_prime, beta + beta_prime, R_P, R_Q)));
}


/************** FILLING THE MATRICES **************/

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

void two_body_F(gsl_vector *c, gsl_matrix *F, R R_A, R R_B){
	double c1=0., c2=0., k=0., val1=0., val2=.0, val3=0., val4=0.;
    int p, q, t, s;
    for(p = 0; p < N; p++){
        for(q = 0; q <= p; q++){
            val1=0., val2=.0, val3=0., val4=0.;
            for(t = 0; t < N; t++){
                for(s = 0; s < N; s++){
                    /* first submatrix */
                    c1 = gsl_vector_get(c, t);
                    c2 = gsl_vector_get(c, s);
                    val1 += c1*c2*direct_term(a[p], a[q], R_A, R_A, a[t], a[s], R_A, R_A);
                    val1 += c1*c2*direct_term(a[p], a[q], R_A, R_A, a[t], a[s], R_A, R_B);
                    val1 += c1*c2*direct_term(a[p], a[q], R_A, R_A, a[t], a[s], R_B, R_A);
                    val1 += c1*c2*direct_term(a[p], a[q], R_A, R_A, a[t], a[s], R_B, R_B);

                    /* second submatrix */
                    c1 = gsl_vector_get(c, t + N);
                    c2 = gsl_vector_get(c, s);
                    val2 += c1*c2*direct_term(a[p], a[q], R_A, R_B, a[t], a[s], R_A, R_A);
                    val2 += c1*c2*direct_term(a[p], a[q], R_A, R_B, a[t], a[s], R_A, R_B);
                    val2 += c1*c2*direct_term(a[p], a[q], R_A, R_B, a[t], a[s], R_B, R_A);
                    val2 += c1*c2*direct_term(a[p], a[q], R_A, R_B, a[t], a[s], R_B, R_B);

                    /* third submatrix */
                    c1 = gsl_vector_get(c, t);
                    c2 = gsl_vector_get(c, s + N);
                    val3 += c1*c2*direct_term(a[p], a[q], R_B, R_A, a[t], a[s], R_A, R_A);
                    val3 += c1*c2*direct_term(a[p], a[q], R_B, R_A, a[t], a[s], R_A, R_B);
                    val3 += c1*c2*direct_term(a[p], a[q], R_B, R_A, a[t], a[s], R_B, R_A);
                    val3 += c1*c2*direct_term(a[p], a[q], R_B, R_A, a[t], a[s], R_B, R_B);

                    /* fourth submatrix */
                    c1 = gsl_vector_get(c, t + N);
                    c2 = gsl_vector_get(c, s + N);
                    val4 += c1*c2*direct_term(a[p], a[q], R_B, R_B, a[t], a[s], R_A, R_A);
                    val4 += c1*c2*direct_term(a[p], a[q], R_B, R_B, a[t], a[s], R_A, R_B);
                    val4 += c1*c2*direct_term(a[p], a[q], R_B, R_B, a[t], a[s], R_B, R_A);
                    val4 += c1*c2*direct_term(a[p], a[q], R_B, R_B, a[t], a[s], R_B, R_B);
                }
            }

            /* first submatrix */
            gsl_matrix_set(F, p, q, val1);  
            gsl_matrix_set(F, q, p, val1);

            /* second submatrix */
            gsl_matrix_set(F, p, q + N, val2);  
            gsl_matrix_set(F, q, p + N, val2);

            /* third submatrix */
            gsl_matrix_set(F, p + N, q, val3);  
            gsl_matrix_set(F, q + N, p, val3);

            /* fourth submatrix */
            gsl_matrix_set(F, p + N, q + N, val4);  
            gsl_matrix_set(F, q + N, p + N, val4);

        }
    }
}